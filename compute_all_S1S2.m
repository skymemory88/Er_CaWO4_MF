function [T, Trans] = compute_all_S1S2(pertResults, const, Options)
% compute_all_S1S2
%   Build S1 (GHz/T) and S2 (GHz/T^2) for a set of transitions, including
%   higher-order transitions |n> -> |n + x| with x in Options.ndE.
%
%   INPUTS
%     pertResults : struct returned by your perturb2nd calls, with fields:
%                   - dE1 (N×3) and dE2 (N×3×3)  [optional, full tensor]
%                   - axis.X/Y/Z.dE1 (N×1), axis.X/Y/Z.dE2 (N×1)  [recommended]
%                   - eigE0 (N×1) energies at B0, in meV  [optional]
%     const       : constants, must include const.Gh2mV (meV per GHz)
%     Options     : struct with fields (all optional)
%                   - ndE : scalar or vector of positive integers (orders).
%                           Default 1 (adjacent transitions only).
%                   - pairs : custom [K×2] list of (m,n) with m>n.
%                             If provided, overrides ndE.
%                   - maxTransitions : limit total transitions calculated (default: unlimited)
%                   - onlyZefoz : if true, pre-filter for low S1 values (default: false)
%                   - zefozThreshold : S1 threshold for pre-filtering (default: 1e-3 GHz/T)
%
%   OUTPUTS
%     T     : table summarizing each transition (sorted for ZEFOZ search)
%     Trans : struct array with detailed S1/S2 per transition
%
%   NOTES
%   * S2 uses the TRUE curvature (Hessian), i.e., 2× the stored coefficient.
%   * Axis-wise S1/S2 are robust at degeneracy. Full tensor S2 is reported
%     only if pertResults.dE2 exists and represents the coefficient tensor.
%   * Memory optimization: Use maxTransitions to limit memory usage for large systems

    if nargin < 3, Options = struct; end
    if ~isfield(Options,'ndE') || isempty(Options.ndE), Options.ndE = 1; end
    if ~isfield(Options,'maxTransitions'), Options.maxTransitions = inf; end
    if ~isfield(Options,'onlyZefoz'), Options.onlyZefoz = false; end
    if ~isfield(Options,'zefozThreshold'), Options.zefozThreshold = 1e-3; end

    % ---------- detect availability ----------
    hasTensor = isfield(pertResults,'dE1') && isfield(pertResults,'dE2') ...
             && ~isempty(pertResults.dE1) && ~isempty(pertResults.dE2) ...
             && ndims(pertResults.dE2)==3 && size(pertResults.dE2,2)==3;

    hasAxis = isfield(pertResults,'axis') && isfield(pertResults.axis,'X') ...
           && isfield(pertResults.axis,'Y') && isfield(pertResults.axis,'Z') ...
           && isfield(pertResults.axis.X,'dE1') && isfield(pertResults.axis.X,'dE2');

    % ---------- number of levels ----------
    if hasTensor
        N = size(pertResults.dE1,1);
    elseif hasAxis
        N = numel(pertResults.axis.Z.dE1);
    else
        error('compute_all_S1S2:NoData', ...
              'Need either full tensor (dE1,dE2) or axis.X/Y/Z results.');
    end

    % ---------- build transition list ----------
    if isfield(Options,'pairs') && ~isempty(Options.pairs)
        pairs = Options.pairs;
        % Validate
        if size(pairs,2) ~= 2, error('pairs must be K×2'); end
        if any(pairs(:,1) <= pairs(:,2))
            error('pairs must have m>n for each row.');
        end
    else
        ndE = unique(Options.ndE(:).');  % row vector
        if any(ndE < 1 | ndE ~= round(ndE))
            error('Options.ndE must be positive integer(s).');
        end
        pairs = [];
        for x = ndE
            if N - x >= 1
                pairs = [pairs; [(1+x: N)' (1: N-x)']];
            end
        end
    end

    Kpairs = size(pairs,1);
    GHz2meV = get_GHz2meV(const);
    haveE = isfield(pertResults,'eigE0') && ~isempty(pertResults.eigE0);

    % Memory optimization: limit transitions if requested
    if Kpairs > Options.maxTransitions
        fprintf('Warning: Limiting transitions from %d to %d for memory efficiency\n', ...
                Kpairs, Options.maxTransitions);
        pairs = pairs(1:Options.maxTransitions, :);
        Kpairs = Options.maxTransitions;
    end

    % ---------- preallocate ----------
    % Pre-allocate arrays for better memory management
    mvec   = zeros(Kpairs,1);
    nvec   = zeros(Kpairs,1);
    order  = zeros(Kpairs,1);
    freqGHZ= nan(Kpairs,1);
    S1norm = nan(Kpairs,1);
    S2fro  = nan(Kpairs,1);
    s1X = nan(Kpairs,1); s1Y = nan(Kpairs,1); s1Z = nan(Kpairs,1);
    s2X = nan(Kpairs,1); s2Y = nan(Kpairs,1); s2Z = nan(Kpairs,1);
    
    % Pre-filter for ZEFOZ if requested (memory optimization)
    if Options.onlyZefoz && hasAxis
        fprintf('Pre-filtering for ZEFOZ transitions (S1 < %.2e GHz/T)...\n', Options.zefozThreshold);
        validPairs = [];
        for k = 1:Kpairs
            m = pairs(k,1); n = pairs(k,2);
            s1_test = sqrt((pertResults.axis.X.dE1(m) - pertResults.axis.X.dE1(n))^2 + ...
                          (pertResults.axis.Y.dE1(m) - pertResults.axis.Y.dE1(n))^2 + ...
                          (pertResults.axis.Z.dE1(m) - pertResults.axis.Z.dE1(n))^2) / GHz2meV;
            if s1_test < Options.zefozThreshold
                validPairs = [validPairs; k];
            end
        end
        pairs = pairs(validPairs, :);
        Kpairs = length(validPairs);
        fprintf('  Retained %d/%d transitions after ZEFOZ pre-filtering\n', Kpairs, size(pairs,1));
        
        % Resize arrays after filtering
        mvec = mvec(1:Kpairs); nvec = nvec(1:Kpairs); order = order(1:Kpairs);
        freqGHZ = freqGHZ(1:Kpairs); S1norm = S1norm(1:Kpairs); S2fro = S2fro(1:Kpairs);
        s1X = s1X(1:Kpairs); s1Y = s1Y(1:Kpairs); s1Z = s1Z(1:Kpairs);
        s2X = s2X(1:Kpairs); s2Y = s2Y(1:Kpairs); s2Z = s2Z(1:Kpairs);
    end

    % Allocate Trans structure only for final number of transitions
    Trans(Kpairs) = struct('m',[],'n',[],'order',[], ...
                           'S1_vec_GHzT',[],'S2_ten_GHzT2',[], ...
                           's1X',[],'s1Y',[],'s1Z',[], ...
                           's2X',[],'s2Y',[],'s2Z',[]);

    % ---------- loop over transitions ----------
    for k = 1:Kpairs
        m = pairs(k,1); n = pairs(k,2); x = m - n;
        mvec(k) = m; nvec(k) = n; order(k) = x;

        % Frequency (optional)
        if haveE
            freqGHZ(k) = (pertResults.eigE0(m) - pertResults.eigE0(n)) / GHz2meV;
        end

        % ---- Full vector/tensor if available (non-degenerate safe region) ----
        if hasTensor
            S1_vec = (squeeze(pertResults.dE1(m,:)) - squeeze(pertResults.dE1(n,:))) / GHz2meV;  % 1×3
            Cm = squeeze(pertResults.dE2(m,:,:));   % meV/T^2, coefficient C
            Cn = squeeze(pertResults.dE2(n,:,:));
            Kmn = 2 * real(Cm - Cn) / GHz2meV;  % 3×3, curvature (GHz/T^2)

            Trans(k).S1_vec_GHzT  = S1_vec(:).';
            Trans(k).S2_ten_GHzT2 = 0.5*(Kmn + Kmn.');   % symmetrize
            S1norm(k) = norm(S1_vec);
            S2fro(k)  = norm(Trans(k).S2_ten_GHzT2,'fro');
        end

        % ---- Axis-wise (always compute; robust at degeneracy) ----
        if hasAxis
            d1x = (pertResults.axis.X.dE1(m) - pertResults.axis.X.dE1(n)) / GHz2meV;
            d1y = (pertResults.axis.Y.dE1(m) - pertResults.axis.Y.dE1(n)) / GHz2meV;
            d1z = (pertResults.axis.Z.dE1(m) - pertResults.axis.Z.dE1(n)) / GHz2meV;
            s1X(k) = d1x; s1Y(k) = d1y; s1Z(k) = d1z;

            d2x = 2*(pertResults.axis.X.dE2(m) - pertResults.axis.X.dE2(n)) / GHz2meV;
            d2y = 2*(pertResults.axis.Y.dE2(m) - pertResults.axis.Y.dE2(n)) / GHz2meV;
            d2z = 2*(pertResults.axis.Z.dE2(m) - pertResults.axis.Z.dE2(n)) / GHz2meV;
            s2X(k) = d2x; s2Y(k) = d2y; s2Z(k) = d2z;

            % If full vector missing, reconstruct S1norm from axis components
            if isnan(S1norm(k))
                S1norm(k) = norm([d1x d1y d1z]);
            end
        end

        % Store into Trans
        Trans(k).m = m; Trans(k).n = n; Trans(k).order = x;
        Trans(k).s1X = s1X(k); Trans(k).s1Y = s1Y(k); Trans(k).s1Z = s1Z(k);
        Trans(k).s2X = s2X(k); Trans(k).s2Y = s2Y(k); Trans(k).s2Z = s2Z(k);
    end

    % ---------- assemble table ----------
    T = table(mvec, nvec, order, freqGHZ, S1norm, S2fro, ...
              abs(s1X), abs(s1Y), abs(s1Z), abs(s2X), abs(s2Y), abs(s2Z), ...
        'VariableNames', {'m','n','order','freq_GHz', ...
                          'S1norm_GHzT','S2fro_GHzT2', ...
                          'abs_s1X','abs_s1Y','abs_s1Z', ...
                          'abs_s2X','abs_s2Y','abs_s2Z'});

    % Sort rows (ZEFOZ-first). If S2fro = NaN (no tensor), sorting ignores it.
    T = sortrows(T, {'S1norm_GHzT','S2fro_GHzT2','abs_s1Z','abs_s2Z'}, ...
                    {'ascend','ascend','ascend','ascend'});
end

% --------- helpers ---------
function val = get_GHz2meV(const)
% get_GHz2meV - Robust unit conversion with consistency checking
%   Returns conversion factor: GHz → meV (i.e., const.Gh2mV)
%   Priority: const.Gh2mV > calculated from fundamentals > CODATA fallback
    
    if isfield(const,'Gh2mV') && ~isempty(const.Gh2mV)
        val = const.Gh2mV;
        validateConversion(val, 'const.Gh2mV');
    else
        % Calculate from fundamental constants if available
        if isfield(const,'hbar') && isfield(const,'J2meV')
            val = const.hbar * 2*pi * 1e9 * const.J2meV;
            fprintf('Using calculated GHz→meV from fundamental constants: %.6e\n', val);
        else
            % Fallback with warning
            val = 4.135667696e-3; % h*1GHz in meV (CODATA 2018)
            fprintf('Warning: Using fallback GHz→meV conversion: %.6e\n', val);
        end
        validateConversion(val, 'calculated/fallback');
    end
end

function validateConversion(val, source)
% Validate unit conversion factor is reasonable
    expected = 4.135667696e-3; % CODATA 2018 value for h*1GHz in meV
    tolerance = 0.01; % 1% tolerance
    
    if abs(val - expected) / expected > tolerance
        fprintf('Warning: %s conversion factor %.6e differs significantly from expected %.6e\n', ...
                source, val, expected);
    end
end
