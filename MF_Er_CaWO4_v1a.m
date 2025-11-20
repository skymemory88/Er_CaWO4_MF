%% Main Script
clearvars

% Setup options
Options.CEF = false; % use crystal field parameters
Options.nEn = 16; % number of eigenstates to include
Options.ndE = [1:15]; % specify which transition orders to calculate and plot
Options.dT_dB = false; % First derivatives of transitions
Options.zefoz = [-14E-6, 14E-6, 1E-5]; % ZEFOZ criteria (field range, criteria)
Options.d2T_dB2 = false; % Second derivatives of transitions
Options.chi = false; % calculate magnetic susceptibility
Options.parallel = false; % use parallelization for susceptibility calculation
Options.freq = linspace(2.6, 3.8, 601); % frequency range for susceptibility [GHz]
Options.temp = 0.1; % temperature for susceptibility calculation [K]
Options.gamma = 1e-5; % spin linewidth [meV]

% Perturbation theory options
Options.perturbation = true; % enable second-order perturbation theory analysis
Options.baseField = [0, 0, 0]; % base magnetic field for perturbation [T]
Options.Bfluctuation = 14e-6; % perturbation amplitude [T]

Options.linSty = {'--', '-', '-.', ':'}; % line styles for different transitions
Options.colors = {'b', 'r', 'g', 'm', 'c', 'k', 'y'}; % colors for different transitions

B0 = linspace(-0.2, 0.2, 401); % [T]
theta = 0; % [rad] deviation angle from c-axis
phi = 0; % [rad] rotation angle within ab-plane, relative to a-axis
Bfield = [B0*sin(theta)*cos(phi); B0*sin(theta)*sin(phi); B0*cos(theta)];

% Initialize
[const, params, hamCF, spinOps, Options] = initialization(Options);

% Hamiltonian diagnolization
[eigE, eigW, dEn] = eigEnergy(Bfield, Options, const, params, hamCF, spinOps);
[~, ~, zefozModes] = deltaE(Options, params.J, params.I, dEn, Bfield); % Transition modes

% Susceptibility calculation
if Options.chi
    fprintf('\n--- Computing Susceptibility ---\n');
    results = MF_susceptibility(Options, const, params, eigE, eigW, Bfield, spinOps);
end

% Perturbation theory analysis using standalone pert_2nd.m
if Options.perturbation
    % Calculate Hamiltonian at base field
    ham0 = buildHamiltonian(Options.baseField, hamCF, spinOps, params, const);

    % Extract unperturbed eigenvalues and eigenvectors (first nEn states)
    [wv0, ee0] = eig(ham0, 'vector');
    ee0 = real(ee0);
    [ee0, sortIdx] = sort(ee0);
    wv0 = wv0(:, sortIdx);
    ee0 = ee0 - ee0(1); % Set ground state to zero energy

    eigE0 = ee0(1:Options.nEn);
    eigW0 = wv0(:, 1:Options.nEn);

    % Electronic and Nuclear Zeeman perturbations
    muBconv = const.muB * const.J2meV;
    muNconv = const.muN * const.J2meV;
    H_pert_x = muBconv * params.gE(1) * spinOps.Jxh + muNconv * params.gN(1) * spinOps.Ixh;
    H_pert_y = muBconv * params.gE(2) * spinOps.Jyh + muNconv * params.gN(2) * spinOps.Iyh;
    H_pert_z = muBconv * params.gE(3) * spinOps.Jzh + muNconv * params.gN(3) * spinOps.Izh;

    % Adaptive degeneracy tolerance (relative to spectral width)
    Espan = max(eig(ham0)) - min(eig(ham0));
    tolE  = max(1e-12, 1e-8 * max(1, Espan));

    % Full-tensor attempt (useful away from degeneracy). If degeneracy exists,
    % per-state K_ij is not unique -> tensors appear in results.blockK.
    pertResults = pert_2nd(ham0, {H_pert_x, H_pert_y, H_pert_z}, ...
        'threshold', tolE, 'verbose', true, 'calcWaveFunc', true, 'calcTensor', true);
    % cellfun(@numel, pertResults.blocks) % debugging

    % ALSO compute along the crystal axes to get well-defined per-state slopes
    % and curvatures even at B0 = 0 (Kramers degeneracy).
    resX = pert_2nd(ham0, {H_pert_x, H_pert_y, H_pert_z}, ...
        'threshold', tolE, 'verbose', false, 'calcWaveFunc', true, 'calcTensor', false, ...
        'direction', [1;0;0]);
    resY = pert_2nd(ham0, {H_pert_x, H_pert_y, H_pert_z}, ...
        'threshold', tolE, 'verbose', false, 'calcWaveFunc', true, 'calcTensor', false, ...
        'direction', [0;1;0]);
    resZ = pert_2nd(ham0, {H_pert_x, H_pert_y, H_pert_z}, ...
        'threshold', tolE, 'verbose', false, 'calcWaveFunc', true, 'calcTensor', false, ...
        'direction', [0;0;1]);

    % Attach axis-wise results to the main struct for plotting convenience
    pertResults.axis.X = resX;  % fields: dE1 (N×1), dE2 (N×1), etc.
    pertResults.axis.Y = resY;
    pertResults.axis.Z = resZ;

    % Optional: noise covariance; either scalar RMS (isotropic) or 3×3 matrix
    sigmaB = 10e-6; % 10 µT RMS
    NoiseCov = sigmaB; % isotropic; or build a 3×3 covariance

    % Visualize only the top ZEFOZ candidates (top 8 for consistent 2×2 layouts):
    [T_top, idxTop] = plot_S1S2_top_zefoz(pertResults, const, Options, ...
        'TopN', 8, 'NoiseCov', NoiseCov);

    % % compute S1 & S2 for selected modes
    % % Options = struct('pairs', nchoosek(1:N,2));  % all m>n
    % [T_all, Trans] = compute_all_S1S2(pertResults, const, Options);

    plotPertb(pertResults, Options, const);

end

fprintf('All calculations completed.\n');

%% Initialization
function [const, params, hamCF, spinOps, Opts] = initialization(Opts)
% Constants
const.hbar = 1.05457E-34;
const.muB = 9.274e-24;
const.kB = 1.3806e-23;
const.muN = 5.05078e-27;
const.mu0 = 4e-7 * pi;
const.kB_meV = 8.61733e-2;
const.J2meV = 6.24151e+21;
const.Gh2mV = const.hbar * 2*pi * 10^9 * const.J2meV;

I = 7/2;
if Opts.CEF
    J = 15/2;
    Bxx = [567.02 -945.07 1008.10 0 -22.16 936.29 0.854] * 0.123983;
    hamCF = cf(J, Bxx, 0);
    hamCF = kron(hamCF, eye(2*I+1));
else
    J = 1/2;
    hamCF = 0;
end

% Parameter set 2
A = [-871.1 -871.1 -130.3] / 1000 * const.Gh2mV; % [meV] hyperfine interaction
gE = [8.3 8.3 1.26]; % electronic g-factor
gN = [-0.1618 -0.1618 -0.1618]; % nuclear g-factor
Q = [1.67 1.67 -3.34] / 1000 * const.Gh2mV; % nuclear quadrupler interaction strength

params.J = J;
params.I = I;
params.A = A;
params.gE = gE;
params.gN = gN;
params.Q = Q;

EnMax = (2*J+1)*(2*I+1);
if Opts.nEn > EnMax; Opts.nEn = EnMax; end

% Pre-calculate spin operators once
[~,~,~,~,~,~,spinOps.Jxh,spinOps.Jyh,spinOps.Jzh,...
    spinOps.Ixh,spinOps.Iyh,spinOps.Izh] = spin_operators(J, I);
end

%% Helper function: Building total Hamiltonian
function ham = buildHamiltonian(Bfield, hamCF, spinOps, params, const)
muBconv = const.muB * const.J2meV;
muNconv = const.muN * const.J2meV;  % CRITICAL FIX: Add proper unit conversion for nuclear term

% Electronic Zeeman term
ham_Ez = muBconv * (params.gE(1) * spinOps.Jxh * Bfield(1) + ...
    params.gE(2) * spinOps.Jyh * Bfield(2) + ...
    params.gE(3) * spinOps.Jzh * Bfield(3));

% Nuclear terms
% Hyperfine interaction
ham_Nz = muNconv * (params.gN(1) * spinOps.Ixh * Bfield(1) + ...
    params.gN(2) * spinOps.Iyh * Bfield(2) + ...
    params.gN(3) * spinOps.Izh * Bfield(3));
% nuclear Zeeman
ham_hyp = params.A(1)*spinOps.Ixh*spinOps.Jxh + params.A(2)*spinOps.Iyh*spinOps.Jyh + params.A(3)*spinOps.Izh*spinOps.Jzh;

% Quadrupolar interaction
ham_Qd = params.Q(1)*spinOps.Ixh*spinOps.Ixh + params.Q(2)*spinOps.Iyh*spinOps.Iyh + params.Q(3)*spinOps.Izh*spinOps.Izh;

ham = hamCF + ham_Ez + ham_Nz + ham_hyp + ham_Qd;
end

%% Eigen-energies and eigen-functions
function [eigE, eigW, dEn] = eigEnergy(Bf, Opts, const, params, hamCF, spinOps)

nb = size(Bf, 2);
% Spin operators now passed in as parameter

eigE = zeros(Opts.nEn, nb);
eigW = zeros(Opts.nEn, Opts.nEn, nb);
dEn = zeros(Opts.nEn, Opts.nEn, nb);

if nb > 50 && license('test', 'Parallel_Computing_Toolbox')
    parfor ii = 1:nb
        [eigE(:,ii), eigW(:,:,ii), dEn(:,:,ii)] = ...
            eigEi(Bf(:,ii), hamCF, spinOps, params, const, Opts);
    end
else
    for ii = 1:nb
        [eigE(:,ii), eigW(:,:,ii), dEn(:,:,ii)] = ...
            eigEi(Bf(:,ii), hamCF, spinOps, params, const, Opts);
    end
end
end

%% Helper function: Hamiltonian diagnolization
function [eigEfield, eigWfield, dEnField] = eigEi(Bfield, hamCF, ...
    spinOps, params, const, Opts)

ham = buildHamiltonian(Bfield, hamCF, spinOps, params, const);

[wv, ee] = eig(ham, 'vector');
ee = real(ee);

[ee, sortIdx] = sort(ee);
wv = wv(:, sortIdx);

ee = ee - ee(1);

eigEfield = ee(1:Opts.nEn);
eigWfield = wv(1:Opts.nEn, 1:Opts.nEn);

enGhz = ee(1:Opts.nEn) / const.Gh2mV;
[Ex, Ey] = meshgrid(enGhz, enGhz);
dEnField = Ex - Ey;
end

%% First and second derivatives
function zefozModes = dT_dB(Opts, deltaEn, Bf, EnMax, modeInfo)
fprintf('\nCalculating transition derivatives (optimized)...\n');

Bz = Bf(3,:);

% Check grid uniformity
dBgrid = diff(Bz);
isUniform = std(dBgrid) < 1e-12 * mean(abs(dBgrid));

if ~isUniform
    warning('Non-uniform grid detected. Using adaptive derivatives.');
end

dT = cell(length(Opts.ndE), 1);
d2T = cell(length(Opts.ndE), 1);

% Process each transition order
for tIdx = 1:length(Opts.ndE)
    transOrd = Opts.ndE(tIdx);
    if transOrd >= EnMax || isempty(deltaEn{tIdx})
        continue;
    end

    dataMat = deltaEn{tIdx};

    % Calculate derivatives based on grid type
    if isUniform
        % Uniform grid derivatives
        [nt, nb] = size(dataMat);

        if nb < 3
            warning('Insufficient points for accurate derivatives');
            dT{tIdx} = NaN(size(dataMat));
            d2T{tIdx} = NaN(size(dataMat));
            continue;
        end

        dB = (Bz(2) - Bz(1)) * 1000; % [mT]

        dT{tIdx} = zeros(nt, nb);

        % Vectorized first derivatives
        if nb >= 3
            % Forward diff for first point
            dT{tIdx}(:, 1) = (-3*dataMat(:,1) + 4*dataMat(:,2) - dataMat(:,3)) / (2*dB);

            % Central diff for interior
            if nb > 2
                dT{tIdx}(:, 2:end-1) = (dataMat(:, 3:end) - dataMat(:, 1:end-2)) / (2*dB);
            end

            % Backward diff for last point
            dT{tIdx}(:, end) = (dataMat(:,end-2) - 4*dataMat(:,end-1) + 3*dataMat(:,end)) / (2*dB);
        end

        % Second derivatives
        if Opts.d2T_dB2 && nb >= 4
            d2T{tIdx} = zeros(nt, nb);
            dB2 = dB^2;

            d2T{tIdx}(:, 1) = (2*dataMat(:,1) - 5*dataMat(:,2) + 4*dataMat(:,3) - dataMat(:,4)) / dB2;

            if nb > 3
                d2T{tIdx}(:, 2:end-1) = (dataMat(:, 3:end) - 2*dataMat(:, 2:end-1) + dataMat(:, 1:end-2)) / dB2;
            end

            d2T{tIdx}(:, end) = (-dataMat(:,end-3) + 4*dataMat(:,end-2) - 5*dataMat(:,end-1) + 2*dataMat(:,end)) / dB2;
        else
            d2T{tIdx} = [];
        end

    else
        % Non-uniform grid derivatives
        [nt, nb] = size(dataMat);

        if nb < 2
            dT{tIdx} = NaN(size(dataMat));
            d2T{tIdx} = NaN(size(dataMat));
            continue;
        end

        dT{tIdx} = zeros(nt, nb);
        BzmT = Bz * 1000;

        for i = 1:nb
            if i == 1
                dB = BzmT(2) - BzmT(1);
                dT{tIdx}(:, i) = (dataMat(:, 2) - dataMat(:, 1)) / dB;
            elseif i == nb
                dB = BzmT(end) - BzmT(end-1);
                dT{tIdx}(:, i) = (dataMat(:, end) - dataMat(:, end-1)) / dB;
            else
                dBback = BzmT(i) - BzmT(i-1);
                dBforw = BzmT(i+1) - BzmT(i);

                w1 = -dBforw / (dBback * (dBback + dBforw));
                w2 = (dBforw - dBback) / (dBback * dBforw);
                w3 = dBback / (dBforw * (dBback + dBforw));

                dT{tIdx}(:, i) = w1 * dataMat(:, i-1) + w2 * dataMat(:, i) + w3 * dataMat(:, i+1);
            end
        end

        if Opts.d2T_dB2 && nb >= 3
            d2T{tIdx} = zeros(nt, nb);

            for i = 2:nb-1
                h1 = BzmT(i) - BzmT(i-1);
                h2 = BzmT(i+1) - BzmT(i);

                w1 = 2 / (h1 * (h1 + h2));
                w2 = -2 / (h1 * h2);
                w3 = 2 / (h2 * (h1 + h2));

                d2T{tIdx}(:, i) = w1 * dataMat(:, i-1) + w2 * dataMat(:, i) + w3 * dataMat(:, i+1);
            end

            d2T{tIdx}(:, 1) = d2T{tIdx}(:, 2);
            d2T{tIdx}(:, end) = d2T{tIdx}(:, end-1);
        else
            d2T{tIdx} = [];
        end
    end

    fprintf('Processed %d transitions for order %d\n', size(dataMat,1), transOrd);
end

% Display derivative statistics
fprintf('\n=== Derivative Statistics ===\n');

for tIdx = 1:length(Opts.ndE)
    transOrd = Opts.ndE(tIdx);
    if transOrd >= EnMax || isempty(dT{tIdx})
        continue;
    end

    derivData = dT{tIdx};
    fprintf('\nTransition order %d:\n', transOrd);

    validMask = ~isnan(derivData);

    for transIdx = 1:size(derivData, 1)
        iSt = modeInfo{tIdx}.iState(transIdx);
        fSt = modeInfo{tIdx}.fStates(transIdx);

        validData = derivData(transIdx, validMask(transIdx, :));
        if ~isempty(validData)
            maxAbs = max(abs(validData));
            meanVal = mean(validData);
            stdVal = std(validData);

            fprintf('  |%d> → |%d>: Max |d/dB| = %.3e, Mean = %.3e ± %.3e GHz/mT\n', ...
                iSt, fSt, maxAbs, meanVal, stdVal);

            if ~isempty(d2T{tIdx})
                secData = d2T{tIdx}(transIdx, validMask(transIdx, :));
                if ~isempty(secData)
                    maxAbs2 = max(abs(secData));
                    meanVal2 = mean(secData);
                    stdVal2 = std(secData);
                    fprintf('    Max |d²/dB²| = %.3e, Mean = %.3e ± %.3e GHz/mT²\n', ...
                        maxAbs2, meanVal2, stdVal2);
                end
            end
        end
    end
end
fprintf('=============================\n\n');

% Plot derivatives
if Opts.dT_dB
    BzPlot = Bf(3,:) * 1000;

    if Opts.d2T_dB2 && ~isempty(d2T)
        figure('Position', [100, 100, 1200, 800]);
        subplotCfg = [2, 1];
    else
        figure('Position', [100, 100, 1200, 400]);
        subplotCfg = [1, 1];
    end

    maxTrans = sum(cellfun(@(x) size(x,1), dT));
    colors = lines(maxTrans);
    colorIdx = 1;

    % Plot first derivatives
    subplot(subplotCfg(1), subplotCfg(2), 1);
    hold on;

    plotHandles = [];
    legendLabels = {};

    for tIdx = 1:length(Opts.ndE)
        transOrd = Opts.ndE(tIdx);
        if transOrd >= EnMax || isempty(dT{tIdx}); continue; end

        dataPlot = dT{tIdx};
        nt = size(dataPlot, 1);

        for transIdx = 1:nt
            plotColor = colors(colorIdx, :);
            colorIdx = colorIdx + 1;

            p = plot(BzPlot, dataPlot(transIdx, :), 'Color', plotColor, 'LineWidth', 1.5);

            iSt = modeInfo{tIdx}.iState(transIdx);
            fSt = modeInfo{tIdx}.fStates(transIdx);
            label = sprintf('d/dB |%d⟩→|%d⟩', iSt, fSt);

            plotHandles(end+1) = p;
            legendLabels{end+1} = label;
        end
    end

    xlabel('Magnetic Field (mT)');
    ylabel('dE/dB (GHz/mT)');
    title('First Derivatives (Optimized)');
    grid on;

    if length(plotHandles) <= 20
        legend(plotHandles, legendLabels, 'Location', 'bestoutside', 'FontSize', 8);
    end

    % Plot second derivatives
    if subplotCfg(1) == 2
        subplot(2, 1, 2);
        hold on;

        plotHandles = [];
        legendLabels = {};
        colorIdx = 1;

        for tIdx = 1:length(Opts.ndE)
            transOrd = Opts.ndE(tIdx);
            if transOrd >= EnMax || isempty(d2T{tIdx}); continue; end

            dataPlot = d2T{tIdx};
            nt = size(dataPlot, 1);

            for transIdx = 1:nt
                plotColor = colors(colorIdx, :);
                colorIdx = colorIdx + 1;

                p = plot(BzPlot, dataPlot(transIdx, :), 'Color', plotColor, 'LineWidth', 1.5);

                iSt = modeInfo{tIdx}.iState(transIdx);
                fSt = modeInfo{tIdx}.fStates(transIdx);
                label = sprintf('d²/dB² |%d⟩→|%d⟩', iSt, fSt);

                plotHandles(end+1) = p;
                legendLabels{end+1} = label;
            end
        end

        xlabel('Magnetic Field (mT)');
        ylabel('d²E/dB² (GHz/mT²)');
        title('Second Derivatives (Optimized)');
        grid on;

        if length(plotHandles) <= 20
            legend(plotHandles, legendLabels, 'Location', 'bestoutside', 'FontSize', 8);
        end
    end

    sgtitle('Transition Derivatives - Vectorized', 'FontSize', 12);
end

% ZEFOZ filtering
zefozModes = perfZEFOZfilt(Opts, deltaEn, dT, d2T, Bf, modeInfo, EnMax);
end


%% Helper function: Initialize ZEFOZ structure
function zefozModes = initZEFOZstruct()
zefozModes = struct();
zefozModes.modeOrd = [];
zefozModes.modeIdx = [];
zefozModes.iSt = [];
zefozModes.fSt = [];
zefozModes.freq = {};
zefozModes.dTdB = {};
zefozModes.d2TdB2 = {};
end

%% Helper function: Vectorized transition filtering
function filtIdx = filtTransVect(derivData, fieldIdx, thresh)
[nt, ~] = size(derivData);
filtIdx = [];

for transIdx = 1:nt
    derivWin = derivData(transIdx, fieldIdx);
    absDerivs = abs(derivWin);

    if all(absDerivs <= thresh)
        filtIdx(end+1) = transIdx;
    end
end
end

%% Helper function: ZEFOZ filtering
function zefozModes = perfZEFOZfilt(Opts, deltaEn, dT, d2T, Bf, modeInfo, EnMax)
zefozModes = [];

if ~isfield(Opts, 'zefoz') || length(Opts.zefoz) < 3
    return;
end

fprintf('\n=== ZEFOZ Filtering ===\n');

Brange = sort(Opts.zefoz(1:2));
Blow = Brange(1);
Bhi = Brange(2);
zefozThresh = abs(Opts.zefoz(3));

fprintf('Field window: [%.6f, %.6f] T ([%.3f, %.3f] mT)\n', ...
    Blow, Bhi, Blow*1000, Bhi*1000);
fprintf('Derivative threshold: ±%.2e GHz/mT\n', zefozThresh);

% Find field indices
B0 = vecnorm(Bf, 2, 1) .* sign(sum(Bf, 1));
fieldMask = (B0 >= Blow) & (B0 <= Bhi);
fieldIdx = find(fieldMask);

if isempty(fieldIdx)
    warning('No field points in ZEFOZ window');
    return;
end

fprintf('Analyzing %d field points...\n', length(fieldIdx));

% Initialize results
zefozModes = initZEFOZstruct();

% Process each transition order
for tIdx = 1:length(Opts.ndE)
    transOrd = Opts.ndE(tIdx);
    if transOrd >= EnMax || isempty(dT{tIdx}); continue; end

    derivData = dT{tIdx};
    freqData = deltaEn{tIdx};

    secDerivData = [];
    if ~isempty(d2T{tIdx})
        secDerivData = d2T{tIdx};
    end

    filtIdx = filtTransVect(derivData, fieldIdx, zefozThresh);

    for locIdx = 1:length(filtIdx)
        transIdx = filtIdx(locIdx);

        zefozModes.modeOrd(end+1) = transOrd;
        zefozModes.modeIdx(end+1) = transIdx;
        zefozModes.iSt(end+1) = modeInfo{tIdx}.iState(transIdx);
        zefozModes.fSt(end+1) = modeInfo{tIdx}.fStates(transIdx);

        zefozModes.freq{end+1} = freqData(transIdx, :);
        zefozModes.dTdB{end+1} = derivData(transIdx, :);

        if ~isempty(secDerivData)
            zefozModes.d2TdB2{end+1} = secDerivData(transIdx, :);
        else
            % Compute simple second derivative using finite differences
            freq_row = freqData(transIdx, :);
            Bz_row = Bf(3,:);
            if length(freq_row) >= 3 && length(Bz_row) >= 3
                % Use central difference for second derivative
                dB = diff(Bz_row);
                if std(dB) < 1e-12 * mean(abs(dB))  % uniform grid
                    dB_uniform = dB(1) * 1000; % convert to mT
                    d2freq = zeros(size(freq_row));
                    d2freq(2:end-1) = (freq_row(3:end) - 2*freq_row(2:end-1) + freq_row(1:end-2)) / dB_uniform^2;
                    d2freq(1) = d2freq(2);  % extrapolate edges
                    d2freq(end) = d2freq(end-1);
                    zefozModes.d2TdB2{end+1} = d2freq;
                else
                    zefozModes.d2TdB2{end+1} = zeros(size(freq_row));  % fallback for non-uniform grid
                end
            else
                zefozModes.d2TdB2{end+1} = zeros(size(freq_row));  % fallback for insufficient points
            end
        end

        fprintf('Found: |%d⟩ → |%d⟩ (Order %d)\n', ...
            zefozModes.iSt(end), zefozModes.fSt(end), transOrd);
    end
end

if ~isempty(zefozModes.modeOrd)
    fprintf('Total ZEFOZ transitions: %d\n', length(zefozModes.modeOrd));
    plotZEFOZres(zefozModes, Bf, fieldIdx, Opts);
else
    fprintf('\nNo ZEFOZ transitions found.\n');
end
end

% Plot ZEFOZ results
function plotZEFOZres(zefozModes, Bf, fieldIdx, Opts)
figure('Position', [100, 50, 1400, 1000]);

Brange = sort(Opts.zefoz(1:2));
Blow = Brange(1) * 1000;
Bhi = Brange(2) * 1000;

BzPlot = Bf(3,:) * 1000;
nModes = length(zefozModes.modeOrd);
colors = lines(nModes);

% Frequencies
subplot(3,1,1); hold on; box on;

yTemp = ylim;
patch([Blow Bhi Bhi Blow], [yTemp(1) yTemp(1) yTemp(2) yTemp(2)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

plotHandles = [];
for i = 1:nModes
    plot(BzPlot, zefozModes.freq{i}, '-', 'Color', [colors(i,:) 0.3], 'LineWidth', 1);
    h = plot(BzPlot(fieldIdx), zefozModes.freq{i}(fieldIdx), ...
        '-', 'Color', colors(i,:), 'LineWidth', 2);
    plotHandles(end+1) = h;
end

ylim auto; yLims = ylim;
delete(findobj(gca, 'Type', 'patch'));
patch([Blow Bhi Bhi Blow], [yLims(1) yLims(1) yLims(2) yLims(2)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
uistack(findobj(gca, 'Type', 'line'), 'top');

xlabel('Magnetic Field (mT)'); ylabel('Frequency (GHz)');
title(sprintf('ZEFOZ Transitions: %d modes', nModes)); grid on;

legendLabels = cell(1, nModes);
for i = 1:nModes
    legendLabels{i} = sprintf('|%d⟩→|%d⟩', zefozModes.iSt(i), zefozModes.fSt(i));
end

if nModes <= 15
    legend(plotHandles, legendLabels, 'Location', 'best', 'FontSize', 9);
end

% First derivatives
subplot(3,1,2); hold on; box on;

thresh = Opts.zefoz(3);
patch([Blow Bhi Bhi Blow], [-abs(thresh) -abs(thresh) abs(thresh) abs(thresh)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'k', 'FaceAlpha', 0.3, 'LineStyle', '--');

for i = 1:nModes
    plot(BzPlot, zefozModes.dTdB{i}, '-', 'Color', [colors(i,:) 0.3], 'LineWidth', 1);
    plot(BzPlot(fieldIdx), zefozModes.dTdB{i}(fieldIdx), '-', 'Color', colors(i,:), 'LineWidth', 2);
end

plot(xlim, [0 0], 'k-', 'LineWidth', 0.5);
xlabel('Magnetic Field (mT)'); ylabel('dE/dB (GHz/mT)');
title('First Derivatives'); grid on;

% Second derivatives
subplot(3,1,3); hold on; box on;

yTemp = [-1 1] * max(cellfun(@(x) max(abs(x)), zefozModes.d2TdB2));
patch([Blow Bhi Bhi Blow], [yTemp(1) yTemp(1) yTemp(2) yTemp(2)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

for i = 1:nModes
    plot(BzPlot, zefozModes.d2TdB2{i}, '-', 'Color', [colors(i,:) 0.3], 'LineWidth', 1);
    plot(BzPlot(fieldIdx), zefozModes.d2TdB2{i}(fieldIdx), '-', 'Color', colors(i,:), 'LineWidth', 2);
end

plot(xlim, [0 0], 'k-', 'LineWidth', 0.5);
xlabel('Magnetic Field (mT)'); ylabel('d²E/dB² (GHz/mT²)');
title('Second Derivatives'); grid on;

sgtitle(sprintf('ZEFOZ Analysis - %d Transitions', nModes), 'FontSize', 14, 'FontWeight', 'bold');
end

%% Transition modes
function [deltaEn, modeInfo, zefozModes] = deltaE(Opts, J, I, dEn, Bf)
EnMax = (2*J+1)*(2*I+1);
if Opts.nEn > EnMax; Opts.nEn = EnMax; end

% Calculate transitions
deltaEn = cell(length(Opts.ndE), 1);
modeInfo = cell(length(Opts.ndE), 1);

for tIdx = 1:length(Opts.ndE)
    transOrd = Opts.ndE(tIdx);
    if transOrd < Opts.nEn
        nt = Opts.nEn - transOrd;
        deltaEn{tIdx} = zeros(nt, size(dEn,3));

        for ii = 1:size(dEn,3)
            deltaEn{tIdx}(:, ii) = diag(squeeze(dEn(:,:,ii)), transOrd);
        end

        modeInfo{tIdx}.iState = (0:nt-1)';
        modeInfo{tIdx}.fStates = (transOrd:Opts.nEn-1)';
    end
end

fprintf('Calculated transitions for orders: %s\n', mat2str(Opts.ndE));

% Plot basic transitions
plotBasicTrans(Opts, deltaEn, Bf, modeInfo);

% Calculate derivatives
zefozModes = [];
if Opts.dT_dB
    zefozModes = dT_dB(Opts, deltaEn, Bf, EnMax, modeInfo);
end
end

%% Plot basic transitions
function plotBasicTrans(Opts, deltaEn, Bf, modeInfo)
figure; hold on; box on;

Bz = Bf(3,:) * 1000;
plotHandles = [];
legendLabels = {};
colors = lines(20);
colorIdx = 1;

for tIdx = 1:length(Opts.ndE)
    if isempty(deltaEn{tIdx}); continue; end

    dataPlot = deltaEn{tIdx};
    nt = size(dataPlot, 1);

    for transIdx = 1:nt
        color = colors(mod(colorIdx-1, size(colors,1))+1, :);
        linSty = mod(transIdx,length(Opts.linSty))+1;
        p = plot(Bz, dataPlot(transIdx, :), Opts.linSty{linSty}, 'Color', color, 'LineWidth', 1.5);

        iSt = modeInfo{tIdx}.iState(transIdx);
        fSt = modeInfo{tIdx}.fStates(transIdx);
        label = sprintf('|%d⟩→|%d⟩', iSt, fSt);

        plotHandles(end+1) = p;
        legendLabels{end+1} = label;
        colorIdx = colorIdx + 1;
    end
end

xlabel('Magnetic Field (mT)'); ylabel('Energy Gap (GHz)');
title('Eigenstate Transitions (Optimized)'); grid on;

if length(plotHandles) <= 25
    legend(plotHandles, legendLabels, 'Location', 'bestoutside', 'FontSize', 9);
end
end


%% CaWO4-specific perturbation analysis and plotting
function plotPertb(pertResults, Options, const)
% Plot CaWO4:Er³⁺ specific perturbation theory results
%
% Inputs:
%   pertResults - Output from pert_2nd() function
%   Options - Options structure with ndE field
%   const - Constants structure
%   params - Parameters structure

fprintf('\n=== CaWO4:Er³⁺ Perturbation Analysis ===\n');

% Extract results from pert_2nd() output
eigE0 = pertResults.eigE0;
eigW0 = pertResults.eigW0;
dE1_tensor = pertResults.dE1;
dE2_tensor = pertResults.dE2;

ns = length(eigE0);

% Prefer axis-wise results if they exist (robust at B0=0)
useAxis = isfield(pertResults, 'axis');
if useAxis
    dE1_axis = [pertResults.axis.X.dE1, pertResults.axis.Y.dE1, pertResults.axis.Z.dE1];  % [ns x 3]
    dE2_axis = [pertResults.axis.X.dE2, pertResults.axis.Y.dE2, pertResults.axis.Z.dE2];  % [ns x 3]
    dE1_GHz  = dE1_axis / const.Gh2mV;   % GHz/T
    dE2_GHz  = dE2_axis / const.Gh2mV;   % GHz/T^2  (these are *coefficients*, curvature = 2*dE2)
else
    dE1_GHz  = dE1_tensor / const.Gh2mV;  % [ns x 3]
    dE2_GHz  = dE2_tensor / const.Gh2mV;  % [ns x 3 x 3]
end

% Create main analysis figure
figure('Position', [100, 100, 1400, 1000]);

% Plot 1: Energy level corrections
subplot(2,2,1);
states = 0:ns-1;
bar(states, dE1_GHz, 'grouped');
xlabel('State Index');
ylabel('Energy Shift (GHz/T)');
title('First-Order Energy Corrections');
legend('dE/dB_x', 'dE/dB_y', 'dE/dB_z', 'Location', 'best');
grid on;

% Plot 2: Second-order diagonal corrections
subplot(2,2,2);
if ndims(dE2_GHz) == 3
    % Full tensor case: extract diagonal elements
    dE2_diag = [squeeze(dE2_GHz(:,1,1)), squeeze(dE2_GHz(:,2,2)), squeeze(dE2_GHz(:,3,3))];
    bar(states, abs(dE2_diag), 'grouped');
    legend('|K_{xx}|', '|K_{yy}|', '|K_{zz}|', 'Location', 'best');
    if all(abs(dE2_diag(:)) < 1e-15)
        text(0.5, 0.5, {'Tensor elements are zero', '(degenerate system)'}, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
else
    % Diagonal-only case: plot per-component diagonal corrections
    bar(states, abs(dE2_GHz), 'grouped');
    legend('|E2_x|', '|E2_y|', '|E2_z|', 'Location', 'best');
end
xlabel('State Index');
ylabel('|E2| (GHz/T²)');
title('Second-Order Diagonal Corrections');
grid on;
if any(abs(dE2_GHz(:)) > 1e-15)
    set(gca, 'YScale', 'log');
else
    ylim([0, 1]);  % Set reasonable scale when all values are zero
end

% Plot 3: Wavefunction probability amplitude changes
subplot(2,2,3);
% Visualize wavefunction probability amplitude changes due to perturbation
% Show probability redistribution among basis states for ground state
ground_state_idx = 1; % Ground state
unperturbed_probs = abs(eigW0(:, ground_state_idx)).^2;

% Define threshold for significant probability states
non_zero_threshold = 1e-6;

% Calculate perturbed probabilities (first-order correction)
if isfield(pertResults, 'dPsi1') && ~isempty(pertResults.dPsi1)
    % Use first-order wavefunction correction for B_z direction (component 3)
    dPsi_z = pertResults.dPsi1{3}; % Z-direction perturbation
    perturbed_wf = eigW0(:, ground_state_idx) + dPsi_z(:, ground_state_idx);
    perturbed_probs = abs(perturbed_wf).^2;
    perturbed_probs = perturbed_probs / sum(perturbed_probs); % Normalize

    % Plot probability changes
    basis_states = 1:ns; % Use 1-based indexing for clarity
    bar(basis_states, [unperturbed_probs, perturbed_probs], 'grouped');
    xlabel('Basis State |J_z, I_z⟩');
    ylabel('Probability Amplitude²');
    title('Ground State Wavefunction Changes (B_z Perturbation)');
    legend({'Unperturbed |ψ₀⟩²', 'Perturbed |ψ₀⁽¹⁾⟩²'}, 'Location', 'best');
    
    % Add basis state labels only for states with non-vanishing amplitudes
    if ns == 16 % J=1/2, I=7/2 case
        % Generate |J_z, I_z⟩ labels for all states
        Jz_vals = [1/2, -1/2]; % J_z = ±1/2
        Iz_vals = 7/2:-1:-7/2;  % I_z = 7/2, 5/2, ..., -7/2
        state_labels = cell(1, ns);
        for j = 1:length(Jz_vals)
            for i = 1:length(Iz_vals)
                idx = (j-1)*length(Iz_vals) + i;
                state_labels{idx} = sprintf('|%+g,%+g⟩', Jz_vals(j), Iz_vals(i));
            end
        end
        
        % Only show labels for states with significant probability
        significant_states = find(max([unperturbed_probs, perturbed_probs], [], 2) > non_zero_threshold);
        if ~isempty(significant_states)
            set(gca, 'XTick', significant_states, 'XTickLabel', state_labels(significant_states));
            xtickangle(45);
        end
    end
else
    % Fallback: show unperturbed probabilities only
    basis_states = 1:ns; % Use 1-based indexing for clarity
    bar(basis_states, unperturbed_probs);
    xlabel('Basis State |J_z, I_z⟩');
    ylabel('Probability Amplitude²');
    title('Ground State Wavefunction (Unperturbed)');
    text(0.5, 0.8, 'Wavefunction corrections not calculated', 'Units', 'normalized');
    
    % Add basis state labels only for states with non-vanishing amplitudes
    if ns == 16 % J=1/2, I=7/2 case
        % Generate |J_z, I_z⟩ labels for all states
        Jz_vals = [1/2, -1/2]; % J_z = ±1/2
        Iz_vals = 7/2:-1:-7/2;  % I_z = 7/2, 5/2, ..., -7/2
        state_labels = cell(1, ns);
        for j = 1:length(Jz_vals)
            for i = 1:length(Iz_vals)
                idx = (j-1)*length(Iz_vals) + i;
                state_labels{idx} = sprintf('|%+g,%+g⟩', Jz_vals(j), Iz_vals(i));
            end
        end
        
        % Only show labels for states with significant probability (unperturbed only in this case)
        significant_states = find(unperturbed_probs > non_zero_threshold);
        if ~isempty(significant_states)
            set(gca, 'XTick', significant_states, 'XTickLabel', state_labels(significant_states));
            xtickangle(45);
        end
    end
end
grid on;

% Plot 4: Total field sensitivity
subplot(2,2,4);
total_linear = sqrt(sum(dE1_GHz.^2, 2));

if ndims(dE2_GHz) == 3
    % Full tensor case: sum over all tensor components
    total_quadratic = sqrt(sum(sum(dE2_GHz.^2, 3), 2));
else
    % Diagonal-only case: sum over components
    total_quadratic = sqrt(sum(dE2_GHz.^2, 2));
end

yyaxis left;
bar(states - 0.2, total_linear, 0.4, 'FaceAlpha', 0.7);
ylabel('Total Linear Response (GHz/T)');
xlabel('State Index');

yyaxis right;
bar(states + 0.2, total_quadratic, 0.4, 'FaceAlpha', 0.7);
ylabel('Total Quadratic Response (GHz/T²)');
title('Combined Field Sensitivity');
legend('Linear', 'Quadratic', 'Location', 'best');
grid on;

sgtitle('CaWO4:Er³⁺ Second-Order Perturbation Analysis', 'FontSize', 14, 'FontWeight', 'bold');

% Calculate and plot transition responses
fprintf('Computing CaWO4-specific transition responses...\n');
plotCaWO4TransitionResponses(pertResults, Options, const);

fprintf('CaWO4:Er³⁺ perturbation analysis completed.\n');
end

function plotCaWO4TransitionResponses(pertResults, Options, const)
% Plot transition frequency responses for CaWO4:Er³⁺ - all orders with multi-figure support

eigE0 = pertResults.eigE0;
ns = length(eigE0);

dE1_tensor = pertResults.dE1;
dE2_tensor = pertResults.dE2;

useAxis = isfield(pertResults, 'axis');
if useAxis
    dE1_axis = [pertResults.axis.X.dE1, pertResults.axis.Y.dE1, pertResults.axis.Z.dE1];  % [ns x 3]
    dE2_axis = [pertResults.axis.X.dE2, pertResults.axis.Y.dE2, pertResults.axis.Z.dE2]; % [ns x 3]
end

% Process all orders in Options.ndE, with 6 subplots per figure
maxSubplots = 6;
validOrders = Options.ndE(Options.ndE < ns); % Only valid transition orders
nValidOrders = length(validOrders);

if nValidOrders == 0
    fprintf('  No valid transition orders to plot.\n');
    return;
end

nFigures = ceil(nValidOrders / maxSubplots);
% Generate multiple figures as needed
for figIdx = 1:nFigures
    % Calculate which orders go in this figure
    startIdx = (figIdx - 1) * maxSubplots + 1;
    endIdx = min(figIdx * maxSubplots, nValidOrders);
    ordersInFig = endIdx - startIdx + 1;

    % Create new figure with offset
    figPos = [150 + (figIdx-1)*50, 150 + (figIdx-1)*50, 1400, 800];
    figure('Position', figPos);

    % Determine subplot layout
    if ordersInFig <= 2
        nRows = 1; nCols = ordersInFig;
    elseif ordersInFig <= 4
        nRows = 2; nCols = 2;
    else % ordersInFig <= 6
        nRows = 2; nCols = 3;
    end

    % Plot transitions for this figure
    for localIdx = 1:ordersInFig
        globalIdx = startIdx + localIdx - 1;
        transitionOrder = validOrders(globalIdx);

        nt = ns - transitionOrder;

        % Calculate linear responses (GHz/T)
        dtrans_x = zeros(nt, 1);
        dtrans_y = zeros(nt, 1);
        dtrans_z = zeros(nt, 1);

        % Calculate quadratic responses (GHz/T²)
        d2trans_x = zeros(nt, 1);
        d2trans_y = zeros(nt, 1);
        d2trans_z = zeros(nt, 1);

        for i = 1:nt
            initial_state = i;
            final_state = i + transitionOrder;

            % Linear responses
            if useAxis
                dtrans_x(i) = (dE1_axis(final_state,1) - dE1_axis(initial_state,1)) / const.Gh2mV;
                dtrans_y(i) = (dE1_axis(final_state,2) - dE1_axis(initial_state,2)) / const.Gh2mV;
                dtrans_z(i) = (dE1_axis(final_state,3) - dE1_axis(initial_state,3)) / const.Gh2mV;
            else
                dtrans_x(i) = (dE1_tensor(final_state,1) - dE1_tensor(initial_state,1)) / const.Gh2mV;
                dtrans_y(i) = (dE1_tensor(final_state,2) - dE1_tensor(initial_state,2)) / const.Gh2mV;
                dtrans_z(i) = (dE1_tensor(final_state,3) - dE1_tensor(initial_state,3)) / const.Gh2mV;
            end


            % Quadratic responses - handle both tensor and diagonal cases
            if useAxis
                d2trans_x(i) = (dE2_axis(final_state,1) - dE2_axis(initial_state,1)) / const.Gh2mV;
                d2trans_y(i) = (dE2_axis(final_state,2) - dE2_axis(initial_state,2)) / const.Gh2mV;
                d2trans_z(i) = (dE2_axis(final_state,3) - dE2_axis(initial_state,3)) / const.Gh2mV;
            else
                if ndims(dE2_tensor) == 3
                    d2trans_x(i) = (dE2_tensor(final_state,1,1) - dE2_tensor(initial_state,1,1)) / const.Gh2mV;
                    d2trans_y(i) = (dE2_tensor(final_state,2,2) - dE2_tensor(initial_state,2,2)) / const.Gh2mV;
                    d2trans_z(i) = (dE2_tensor(final_state,3,3) - dE2_tensor(initial_state,3,3)) / const.Gh2mV;
                else
                    d2trans_x(i) = (dE2_tensor(final_state,1) - dE2_tensor(initial_state,1)) / const.Gh2mV;
                    d2trans_y(i) = (dE2_tensor(final_state,2) - dE2_tensor(initial_state,2)) / const.Gh2mV;
                    d2trans_z(i) = (dE2_tensor(final_state,3) - dE2_tensor(initial_state,3)) / const.Gh2mV;
                end
            end

        end

        subplot(nRows, nCols, localIdx);

        % Plot linear responses
        trans_idx = 0:nt-1;
        plot(trans_idx, dtrans_x, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'dν/dB_x');
        hold on;
        plot(trans_idx, dtrans_y, 'gs-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'dν/dB_y');
        plot(trans_idx, dtrans_z, 'b^-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'dν/dB_z');

        xlabel('Transition Index');
        ylabel('Linear Response (GHz/T)');
        title(sprintf('Order %d Transitions', transitionOrder));
        legend('show', 'Location', 'best', 'FontSize', 8);
        grid on;
    end

    % Set figure title
    if nFigures == 1
        sgtitle('CaWO4:Er³⁺ Transition Frequency Responses', 'FontSize', 14, 'FontWeight', 'bold');
    else
        sgtitle(sprintf('CaWO4:Er³⁺ Transition Responses - Figure %d/%d', figIdx, nFigures), ...
            'FontSize', 14, 'FontWeight', 'bold');
    end
end

% Also create a quadratic response figure (optional - can be enabled)
if nValidOrders <= 6
    d2TdB2_plot(pertResults, Options, const, validOrders);
end
end

function d2TdB2_plot(pertResults, ~, const, validOrders)
% Plot quadratic (second-order) transition responses for CaWO4:Er³⁺

eigE0 = pertResults.eigE0;
dE2_tensor = pertResults.dE2;
ns = length(eigE0);
nValidOrders = length(validOrders);

% Create quadratic response figure
figure('Position', [200, 200, 1400, 800]);

% Determine subplot layout
if nValidOrders <= 2
    nRows = 1; nCols = nValidOrders;
elseif nValidOrders <= 4
    nRows = 2; nCols = 2;
else % nValidOrders <= 6
    nRows = 2; nCols = 3;
end

for orderIdx = 1:nValidOrders
    transitionOrder = validOrders(orderIdx);
    nt = ns - transitionOrder;

    % Calculate quadratic responses (GHz/T²)
    d2trans_x = zeros(nt, 1);
    d2trans_y = zeros(nt, 1);
    d2trans_z = zeros(nt, 1);

    for i = 1:nt
        initial_state = i;
        final_state = i + transitionOrder;

        % Quadratic responses - handle both tensor and diagonal cases
        if ndims(dE2_tensor) == 3
            d2trans_x(i) = (dE2_tensor(final_state, 1, 1) - dE2_tensor(initial_state, 1, 1)) / const.Gh2mV;
            d2trans_y(i) = (dE2_tensor(final_state, 2, 2) - dE2_tensor(initial_state, 2, 2)) / const.Gh2mV;
            d2trans_z(i) = (dE2_tensor(final_state, 3, 3) - dE2_tensor(initial_state, 3, 3)) / const.Gh2mV;
        else
            % For diagonal-only case, use per-component corrections
            d2trans_x(i) = (dE2_tensor(final_state, 1) - dE2_tensor(initial_state, 1)) / const.Gh2mV;
            d2trans_y(i) = (dE2_tensor(final_state, 2) - dE2_tensor(initial_state, 2)) / const.Gh2mV;
            d2trans_z(i) = (dE2_tensor(final_state, 3) - dE2_tensor(initial_state, 3)) / const.Gh2mV;
        end
    end

    subplot(nRows, nCols, orderIdx);

    % Plot quadratic responses on log scale
    trans_idx = 0:nt-1;
    semilogy(trans_idx, abs(d2trans_x), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '|d²ν/dB²_x|');
    hold on;
    semilogy(trans_idx, abs(d2trans_y), 'gs-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '|d²ν/dB²_y|');
    semilogy(trans_idx, abs(d2trans_z), 'b^-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '|d²ν/dB²_z|');

    xlabel('Transition Index');
    ylabel('|Quadratic Response| (GHz/T²)');
    title(sprintf('Order %d Quadratic Transitions', transitionOrder));
    legend('show', 'Location', 'best', 'FontSize', 8);
    grid on;
end

sgtitle('CaWO4:Er³⁺ Quadratic Transition Frequency Responses', 'FontSize', 14, 'FontWeight', 'bold');
end