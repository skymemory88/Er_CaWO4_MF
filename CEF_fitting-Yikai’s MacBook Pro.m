function [B_fit_cm, g_eff_fit, stats] = CEF_fitting(eigenE, fieldData)
%% CEF g-tensor fitting for Er:CaWO4
% Fits the crystal-field Stevens coefficients so the lowest Kramers doublet
% reproduces a target effective-spin g tensor, while also matching the
% provided zero-field 16-level hyperfine spectrum (ground-state doublet
% coupled to I=7/2). The spectrum is used as an additional residual term
% and is treated up to a global energy scale via an optimal scalar A0.
% Optional field-dependent spectra can be supplied to further constrain the
% fit; these should be provided as 16-level ladders evaluated at known
% magnetic-field vectors (matching the layout returned by MF_Er_CaWO4_v1b).
%
% Input
%   eigenE : 1x16 or 16x1 vector of zero-field energies (meV or MHz OK).
%            The residual uses a linear scale fit, so absolute units cancel.
%            Values are sorted internally and mean-centered for matching.
%   fieldData (optional) : struct or struct array with fields
%       .Bvec      : 3xN array of field vectors (Tesla) for each spectrum
%       .energies  : 16xN array of corresponding energy levels (same units
%                    as eigenE)
%       .weight    : scalar residual weight (default 1)
%       .gN        : 1x3 nuclear g tensor (default [0.1618 0.1618 0.1618])
%       .Q         : 1x3 quadrupole coefficients in energy units (default 0)
%
% Output
%   B_fit_cm   : best-fit Stevens coefficients (cm^-1)
%   g_eff_fit  : resulting effective g tensor [gx gy gz]
%   stats      : struct with fields resnorm, exitflag, output, A0_best, funcCount
%
% Notes
% - No optimizer changes beyond extending the residual; existing BayesOpt/
%   lsqnonlin/fminsearch flows are preserved.
% - The hyperfine residual is built from H_hf = A0 * (Jproj.Jx⊗Ix + Jy⊗Iy + Jz⊗Iz),
%   where Jproj.* are the projected J operators of the lowest doublet from the
%   current CEF trial, and I=7/2 for 167Er. A0 is eliminated per-evaluation via
%   a scalar least-squares fit to Eobs16 so energy units need not be specified.

if nargin < 2
    fieldData = [];
end

%% Target and search configuration -----------------------------------------
targetG = [8.3, 8.3, 1.26]; % target pseudospin g tensor
targetAeff_MHz = [-871.1, -871.1, -130.3]; % optional: target effective hyperfine (MHz)
reportHyperfine = false; % if true, print inferred A0 and A_eff
baseB_cm = [-753.279 796.633 -376.577 -133.463 -3.30013 -84.5397 -9.7149];
fitMask = [true, true, true, true, true, true, true]; % true -> parameter is varied

L = 6; % Er
S = 3/2; % Er

option.Bayes = true; % enable Bayesian optimisation instead of LM (if available)
option.bounds = true; % toggle simple box bounds in cm^-1
option.polish = true; % run a local lsqnonlin polish after BayesOpt if available
lowerB_cm = baseB_cm - 400; % crude bounds to keep the search stable
upperB_cm = baseB_cm + 400;
lowerB_cm(~fitMask) = baseB_cm(~fitMask);
upperB_cm(~fitMask) = baseB_cm(~fitMask);

lsqAvailable = exist('lsqnonlin', 'file') == 2;
if lsqAvailable
    lsqOpts = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
        'StepTolerance', 1e-3, 'FunctionTolerance', 1e-3, 'MaxIterations', 1e4);
    if option.bounds
        freeLower = lowerB_cm(fitMask);
        freeUpper = upperB_cm(fitMask);
    else
        freeLower = [];
        freeUpper = [];
    end
else
    lsqOpts = [];
    freeLower = [];
    freeUpper = [];
end

bayesSummary = [];

fieldData = normalizeFieldData(fieldData);

%% Optimiser call -----------------------------------------------------------
free0 = baseB_cm(fitMask);

% Prepare observed 16-level spectrum: sort and mean-center (units cancel)
eigenE = eigenE(:).';
assert(numel(eigenE) == 16, 'eigen-energy must have 16 elements.');
[eigE, ~] = sort(eigenE);
eigE_Norm = eigE - min(eigE); % normalize the eigen-energies

% Precompute nuclear spin operators once (I = 7/2)
I = 7/2;
[~, ~, ~, Ix, Iy, Iz, ~, ~, ~, ~, ~, ~] = spin_operators(1/2, I);

% Residual combines g-tensor misfit and zero-field 16-level spectrum
% Weighting kept simple and dimensionless
wG = 1.0;  % g residual weight
wE = 1.0;  % energy residual weight
objective = @(freeB) fitResidual(freeB, baseB_cm, fitMask, targetG, L, S, ...
                                 Ix, Iy, Iz, eigE_Norm, wG, wE, fieldData);

if option.Bayes && exist('bayesopt', 'file') == 2
    fprintf('Running Bayesian optimisation over the selected CEF coefficients...\n');
    optVars = optimizableVariable.empty;
    idx = 0;
    for jj = 1:numel(baseB_cm)
        if fitMask(jj)
            idx = idx + 1;
            name = sprintf('B%d', jj);
            if option.bounds
                bounds = [lowerB_cm(jj), upperB_cm(jj)];
            else
                span = max(800, abs(baseB_cm(jj)) + 800);
                bounds = [baseB_cm(jj) - span, baseB_cm(jj) + span];
            end
            optVars(idx) = optimizableVariable(name, bounds); 
        end
    end

    meritFcn = @(tbl) bayesObjective(tbl, baseB_cm, fitMask, targetG, L, S, Ix, Iy, Iz, eigE_Norm, wG, wE, fieldData);
    results = bayesopt(meritFcn, optVars, ...
        'MaxObjectiveEvaluations', 1e4, ...
        'IsObjectiveDeterministic', true, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 0, ...
        'PlotFcn', []);

    bestPoint = results.XAtMinObjective;
    freeOpt = zeros(sum(fitMask), 1);
    names = bestPoint.Properties.VariableNames;
    for kk = 1:numel(names)
        freeOpt(kk) = bestPoint{1, names{kk}};
    end
    resnorm = results.MinObjective;
    exitflag = 0;
    output = struct('funcCount', results.NumObjectiveEvaluations, ...
                    'iterations', results.NumObjectiveEvaluations, ...
                    'phase', 'bayesopt');
    bayesSummary = struct('MinObjective', results.MinObjective, ...
                          'NumObjectiveEvaluations', results.NumObjectiveEvaluations);
    residual = sqrt(resnorm) * sign(resnorm); 
    if option.polish && lsqAvailable
        fprintf('Polishing BayesOpt solution with lsqnonlin...\n');
        [freePolish, resnormPolish, residualPolish, exitflagPolish, outputPolish] = ...
            lsqnonlin(objective, freeOpt, freeLower, freeUpper, lsqOpts);
        if exitflagPolish > 0 && resnormPolish <= resnorm
            freeOpt = freePolish;
            resnorm = resnormPolish;
            residual = residualPolish;
            exitflag = exitflagPolish;
            output = outputPolish;
            output.phase = 'lsqnonlin';
            output.bayesEvaluations = results.NumObjectiveEvaluations;
        else
            output.polish = outputPolish;
            output.polishExitflag = exitflagPolish;
            output.polishResnorm = resnormPolish;
        end
    end
elseif lsqAvailable
    [freeOpt, resnorm, residual, exitflag, output] = lsqnonlin(objective, free0, freeLower, freeUpper, lsqOpts);
    output.phase = 'lsqnonlin';
else
    fprintf('lsqnonlin not available; falling back to fminsearch.\n');
    simplexObj = @(freeB) sum(objective(freeB).^2);
    [freeOpt, resnorm, exitflag, output] = fminsearch(simplexObj, free0, optimset('Display', 'iter'));
    residual = objective(freeOpt);
    if isstruct(output)
        output.phase = 'fminsearch';
    end
end

%% Results ------------------------------------------------------------------
B_fit_cm = baseB_cm;
B_fit_cm(fitMask) = freeOpt;

[g_eff_fit, Jproj] = computeGeff(B_fit_cm, L, S);
misfit = (g_eff_fit - targetG) ./ targetG;

fprintf('\nOptimisation exit flag: %g\n', exitflag);
fprintf('Final resnorm: %.6g\n', resnorm);
if isstruct(output) && isfield(output, 'funcCount')
    fprintf('Function evaluations: %d\n', output.funcCount);
end
if isstruct(output) && isfield(output, 'iterations')
    fprintf('Iterations: %d\n', output.iterations);
end
fprintf('\nBest-fit CEF parameters (cm^-1):\n    %s\n', mat2str(B_fit_cm, 6));
fprintf('Resulting g tensor: [%.4f %.4f %.4f]\n', g_eff_fit);
fprintf('Relative g error:  [%.2e %.2e %.2e]\n', misfit);

fprintf('\nProjected J operators in the doublet (dimensionless):\n');
fprintf('Jx -> off-diagonal magnitude %.4f\n', abs(Jproj.Jx(1,2)));
fprintf('Jy -> off-diagonal magnitude %.4f\n', abs(Jproj.Jy(1,2)));
fprintf('Jz -> diagonal values [%.4f, %.4f]\n', Jproj.Jz(1,1), Jproj.Jz(2,2));

if reportHyperfine
    gJ = gLande(L, S);
    [A0_axes_MHz, A0_rec_MHz, Aeff_fromA0_MHz] = inferA0AndAeff(g_eff_fit, gJ, targetAeff_MHz);
    fprintf('\nHyperfine mapping via isotropic A0 (J dot I):\n');
    fprintf('A0 estimates per axis (MHz): [%.2f %.2f %.2f]\n', A0_axes_MHz);
    fprintf('Recommended A0 (MHz): %.2f\n', A0_rec_MHz);
    fprintf('Implied A_eff from A0 and fitted g (MHz): [%.1f %.1f %.1f]\n', Aeff_fromA0_MHz);
end

% Report the best-fit scalar A0 (in same units as Eobs16 input)
Hhf_unit = kron(Jproj.Jx, Ix) + kron(Jproj.Jy, Iy) + kron(Jproj.Jz, Iz);
evals = eig((Hhf_unit + Hhf_unit')/2);
evals = sort(real(evals));
evals = evals - mean(evals);
A0_best = (sum(evals .* eigE_Norm(:))) / (sum(evals .* evals) + eps);
fprintf('\nZero-field hyperfine scale A0 (units of Eobs16): %.6g\n', A0_best);

% Collect stats
stats = struct();
stats.resnorm = resnorm;
stats.exitflag = exitflag;
stats.output = output;
if ~isempty(bayesSummary)
    stats.bayes = bayesSummary;
end
if isstruct(output) && isfield(output, 'funcCount')
    stats.funcCount = output.funcCount;
elseif exist('results', 'var')
    stats.funcCount = results.NumObjectiveEvaluations; %#ok<*NASGU>
end
stats.A0_best = A0_best;

%% Helper functions ---------------------------------------------------------
function res = fitResidual(freeB, baseB, mask, targetG, L, S, Ix, Iy, Iz, Eobs16c, wG, wE, fieldData)
    % Assemble trial CF set
    Btrial = baseB;
    Btrial(mask) = freeB;

    % Compute effective g and projected J for the lowest doublet
    [g_eff, Jproj] = computeGeff(Btrial, L, S);

    % g residual (dimensionless, compare by components)
    gres = (g_eff(:) - targetG(:)) ./ targetG(:);

    % Zero-field hyperfine spectrum in doublet ⊗ I space (16x16)
    % Build dimensionless hyperfine Hamiltonian for A0 = 1
    Hhf_unit = kron(Jproj.Jx, Ix) + kron(Jproj.Jy, Iy) + kron(Jproj.Jz, Iz);
    evals = eig((Hhf_unit + Hhf_unit')/2); % enforce Hermiticity numerically
    evals = sort(real(evals));
    evals = evals - mean(evals);

    % Optimal scalar A0 (same units as Eobs16c) by linear least squares
    num = sum(evals .* Eobs16c(:));
    den = sum(evals .* evals) + eps;
    A0_opt = num / den;

    % Energy residual (dimensionless via RMS scaling of observed)
    epred = A0_opt * evals;
    eobs = Eobs16c(:);
    scale = sqrt(mean(eobs.^2)) + eps;
    eres = (epred - eobs) / scale;

    % Field-dependent spectra, if supplied
    fieldRes = computeFieldResidual(fieldData, Jproj, Hhf_unit, A0_opt, Ix, Iy, Iz, L, S);

    % Stack with weights
    res = [wG * gres; wE * eres; fieldRes];
end

function fieldRes = computeFieldResidual(fieldData, Jproj, Hhf_unit, A0_opt, Ix, Iy, Iz, L, S)
    if isempty(fieldData)
        fieldRes = [];
        return;
    end

    muB_meV_T = 9.2740100783e-24 * 6.241509074e21;
    muN_meV_T = 5.0507837461e-27 * 6.241509074e21;
    gJ = gLande(L, S);

    IdJ = eye(size(Jproj.Jx));
    IxOp = kron(IdJ, Ix);
    IyOp = kron(IdJ, Iy);
    IzOp = kron(IdJ, Iz);
    Ix2 = Ix * Ix;
    Iy2 = Iy * Iy;
    Iz2 = Iz * Iz;

    fieldRes = [];
    for dd = 1:numel(fieldData)
        fd = fieldData(dd);
        Bvec = fd.Bvec;
        Eobs = fd.energies;
        if size(Bvec,1) ~= 3
            error('fieldData(%d).Bvec must be a 3xN array of field components.', dd);
        end
        if size(Eobs,1) ~= 16
            error('fieldData(%d).energies must be a 16xN array.', dd);
        end

        gN = getStructField(fd, 'gN', [0.1618, 0.1618, 0.1618]);
        Q = getStructField(fd, 'Q', [0, 0, 0]);
        weight = getStructField(fd, 'weight', 1);

        for kk = 1:size(Bvec,2)
            Bk = Bvec(:, kk);

            He = -muB_meV_T * gJ * (Bk(1)*Jproj.Jx + Bk(2)*Jproj.Jy + Bk(3)*Jproj.Jz);
            He = (He + He') / 2;
            Hze = kron(He, eye(size(Ix,1)));

            Hn = -muN_meV_T * (Bk(1) * gN(1) * IxOp + ...
                               Bk(2) * gN(2) * IyOp + ...
                               Bk(3) * gN(3) * IzOp);
            Hq = Q(1) * kron(IdJ, Ix2) + ...
                 Q(2) * kron(IdJ, Iy2) + ...
                 Q(3) * kron(IdJ, Iz2);

            Htot = Hze + A0_opt * Hhf_unit + Hn + Hq;
            Htot = (Htot + Htot') / 2;

            evals = eig(Htot);
            evals = sort(real(evals));
            evals = evals - mean(evals);

            obs = Eobs(:, kk);
            obs = obs(:) - mean(obs);
            scale = sqrt(mean(obs.^2)) + eps;

            fieldRes = [fieldRes; weight * ((evals - obs) / scale)];
        end
    end
end

function val = getStructField(s, name, defaultVal)
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = defaultVal;
    end
end

function fd = normalizeFieldData(fd)
    if isempty(fd)
        return;
    end
    if ~isstruct(fd)
        error('fieldData must be provided as a struct or struct array.');
    end
    for jj = 1:numel(fd)
        if isfield(fd(jj), 'B') && ~isfield(fd(jj), 'Bvec')
            fd(jj).Bvec = fd(jj).B;
        end
        if ~isfield(fd(jj), 'Bvec')
            error('fieldData(%d) is missing the Bvec field.', jj);
        end
        if ~isfield(fd(jj), 'energies')
            error('fieldData(%d) is missing the energies field.', jj);
        end
        fd(jj).Bvec = double(fd(jj).Bvec);
        fd(jj).energies = double(fd(jj).energies);
    end
end

function merit = bayesObjective(tbl, baseB, mask, targetG, L, S, Ix, Iy, Iz, Eobs16c, wG, wE, fieldData)
    freeB = zeros(sum(mask), 1);
    names = tbl.Properties.VariableNames;
    for kk = 1:numel(names)
        freeB(kk) = tbl{1, names{kk}};
    end
    res = fitResidual(freeB, baseB, mask, targetG, L, S, Ix, Iy, Iz, Eobs16c, wG, wE, fieldData);
    merit = sum(res.^2);
end

function [g_eff, Jproj] = computeGeff(B_cm, L, S)
    conv_cm_to_meV = 0.123983; % 1 cm^-1 to meV
    J = L + S;
    gJ = gLande(L, S);

    B_meV = B_cm * conv_cm_to_meV;
    Hcf = cf(J, B_meV, 0);
    [~, ~, ~, ~, ~, ~, Jx, Jy, Jz, ~, ~, ~] = spin_operators(J, 0);

    [~, Jproj] = projectDoublet(Hcf, Jx, Jy, Jz);
    g_eff = zeros(1, 3);
    g_eff(1) = 2 * gJ * abs(Jproj.Jx(1, 2));
    g_eff(2) = 2 * gJ * abs(Jproj.Jy(1, 2));
    g_eff(3) = 2 * gJ * abs(Jproj.Jz(1, 1));
end
end
