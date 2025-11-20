function [Tdir, s1u_all, s2u_all] = compute_all_S1S2_along_direction(pertResults, const, uhat, varargin)
% compute_all_S1S2_along_direction
%   For each transition, compute S1 and S2 along uhat (GHz/T and GHz/T^2).
%
% Usage:
%   [Tdir, s1u, s2u] = compute_all_S1S2_along_direction(pertResults, const, uhat);           % adjacent
%   [Tdir, s1u, s2u] = compute_all_S1S2_along_direction(pertResults, const, uhat, 'pairs','all');

    p = inputParser;
    addParameter(p,'pairs','adjacent');  % 'adjacent' | 'all' | [K×2]
    parse(p,varargin{:});
    opt = p.Results;

    % Number of levels and pairs
    if isfield(pertResults,'dE1') && ~isempty(pertResults.dE1)
        N = size(pertResults.dE1,1);
    elseif isfield(pertResults,'axis') && isfield(pertResults.axis,'Z') && ~isempty(pertResults.axis.Z.dE1)
        N = numel(pertResults.axis.Z.dE1);
    else
        error('Cannot infer number of levels.');
    end

    if ischar(opt.pairs) || isstring(opt.pairs)
        if strcmpi(opt.pairs,'adjacent')
            pairs = [(2:N)' (1:N-1)'];
        elseif strcmpi(opt.pairs,'all')
            [mm, nn] = find(triu(true(N),1));
            pairs = [mm nn];
        else
            error('pairs must be ''adjacent'', ''all'', or a [K×2] array.');
        end
    else
        pairs = opt.pairs;
    end

    K = size(pairs,1);
    s1u_all = zeros(K,1);
    s2u_all = zeros(K,1);

    for k = 1:K
        m = pairs(k,1); n = pairs(k,2);
        [s1u_all(k), s2u_all(k)] = transition_S1S2_direction(pertResults, m, n, const, uhat);
    end

    Tdir = table(pairs(:,1), pairs(:,2), abs(s1u_all), abs(s2u_all), ...
                 'VariableNames', {'m','n','abs_s1u_GHzT','abs_s2u_GHzT2'});
    Tdir = sortrows(Tdir, {'abs_s1u_GHzT','abs_s2u_GHzT2'}, {'ascend','ascend'});
end

% helper function
function [s1_u_GHzT, s2_u_GHzT2] = transition_S1S2_direction(pertResults, m, n, const, uhat)
% transition_S1S2_direction
%   Returns S1 and S2 along an arbitrary unit vector uhat (3×1 or 1×3).
%   Prefers full tensor when available; otherwise uses axis.* only for
%   cardinal directions.

    u = uhat(:);  u = u / norm(u);
    GHz2meV = get_GHz2meV(const);

    hasTensor = isfield(pertResults,'dE1') && isfield(pertResults,'dE2') ...
                && ndims(pertResults.dE2) == 3 && size(pertResults.dE2,2) == 3;

    if hasTensor
        dm1 = squeeze(pertResults.dE1(m,:)).';   % 3×1
        dn1 = squeeze(pertResults.dE1(n,:)).';
        S1_vec = (dm1 - dn1) / GHz2meV;          % GHz/T
        s1_u_GHzT = dot(S1_vec, u);

        Cm = squeeze(pertResults.dE2(m,:,:));    % meV/T^2, coefficient
        Cn = squeeze(pertResults.dE2(n,:,:));
        Kmn = 2 * real(Cm - Cn);                 % Hessian, meV/T^2
        S2_ten = Kmn / GHz2meV;                  % GHz/T^2
        s2_u_GHzT2 = u.' * S2_ten * u;

        return
    end

    % No tensor: allow X/Y/Z only
    if all(abs(u-[1;0;0])<1e-12) || all(abs(u-[-1;0;0])<1e-12)
        [s1_u_GHzT, s2_u_GHzT2] = transition_S1S2_axis(pertResults, m, n, const, 'X'); return
    elseif all(abs(u-[0;1;0])<1e-12) || all(abs(u-[0;-1;0])<1e-12)
        [s1_u_GHzT, s2_u_GHzT2] = transition_S1S2_axis(pertResults, m, n, const, 'Y'); return
    elseif all(abs(u-[0;0;1])<1e-12) || all(abs(u-[0;0;-1])<1e-12)
        [s1_u_GHzT, s2_u_GHzT2] = transition_S1S2_axis(pertResults, m, n, const, 'Z'); return
    else
        error(['No full tensor available and direction is not a cardinal axis. ', ...
               'Call perturb2nd(...,''direction'',uhat) for this uhat and attach it as pertResults.axis.U.']);
    end
end

function [s1_GHzT, s2_GHzT2] = transition_S1S2_axis(pertResults, m, n, const, axisChar)
% transition_S1S2_axis
%   Calculate S1 and S2 for a specific axis using axis.X/Y/Z data
    
    GHz2meV = get_GHz2meV(const);
    
    switch upper(axisChar)
        case 'X'
            if ~isfield(pertResults, 'axis') || ~isfield(pertResults.axis, 'X')
                error('pertResults.axis.X not available for axis calculation');
            end
            dE1_axis = pertResults.axis.X.dE1;
            dE2_axis = pertResults.axis.X.dE2;
        case 'Y' 
            if ~isfield(pertResults, 'axis') || ~isfield(pertResults.axis, 'Y')
                error('pertResults.axis.Y not available for axis calculation');
            end
            dE1_axis = pertResults.axis.Y.dE1;
            dE2_axis = pertResults.axis.Y.dE2;
        case 'Z'
            if ~isfield(pertResults, 'axis') || ~isfield(pertResults.axis, 'Z')
                error('pertResults.axis.Z not available for axis calculation');
            end
            dE1_axis = pertResults.axis.Z.dE1;
            dE2_axis = pertResults.axis.Z.dE2;
        otherwise
            error('axisChar must be ''X'', ''Y'', or ''Z''');
    end
    
    % Calculate transition sensitivities
    s1_GHzT = (dE1_axis(m) - dE1_axis(n)) / GHz2meV;
    s2_GHzT2 = 2 * (dE2_axis(m) - dE2_axis(n)) / GHz2meV;  % Factor of 2 for true curvature
end

function val = get_GHz2meV(const)
% get_GHz2meV - Unit conversion with validation
%   Returns const.Gh2mV (GHz → meV conversion factor)
    if isfield(const,'Gh2mV') && ~isempty(const.Gh2mV)
        val = const.Gh2mV;
    elseif isfield(const,'hbar') && isfield(const,'J2meV')
        val = const.hbar * 2*pi * 1e9 * const.J2meV;
    else
        val = 4.135667696e-3; % CODATA 2018 fallback: h*1GHz in meV
    end
    
    % Validate conversion factor
    expected = 4.135667696e-3;
    if abs(val - expected) / expected > 0.01
        warning('GHz→meV conversion factor %.6e differs from expected %.6e', val, expected);
    end
end