function A0_MHz = calculate_A0_from_Aeff(A_eff_MHz, g_eff, gJ)
%CALCULATE_A0_FROM_AEFF Convert effective anisotropic A to isotropic A0
%   A0_MHz = calculate_A0_from_Aeff(A_eff_MHz, g_eff, gJ)
%
%   The relationship is: A_eff_i = A0 * g_eff_i / (2*gJ)
%   Therefore: A0 = A_eff_i * (2*gJ) / g_eff_i
%
%   Inputs:
%     A_eff_MHz : 1x3 vector of effective hyperfine constants [Ax, Ay, Az] in MHz
%     g_eff     : 1x3 vector of effective g-factors [gx, gy, gz]
%     gJ        : Lande g-factor for the ion
%
%   Output:
%     A0_MHz    : Isotropic hyperfine constant in MHz
%
%   This function calculates A0 for each component and returns the mean.
%   A warning is issued if the components differ by more than 5%.

    if nargin < 3
        error('calculate_A0_from_Aeff requires 3 inputs: A_eff_MHz, g_eff, gJ');
    end

    % Ensure inputs are the correct size
    if numel(A_eff_MHz) ~= 3 || numel(g_eff) ~= 3
        error('A_eff_MHz and g_eff must be 1x3 or 3x1 vectors');
    end
    if numel(gJ) ~= 1
        error('gJ must be a scalar');
    end

    % Calculate A0 from each component
    A0_components = A_eff_MHz(:)' .* (2*gJ) ./ g_eff(:)';

    % Return the mean
    A0_MHz = mean(A0_components);

    % Check consistency
    rel_std = std(A0_components) / abs(mean(A0_components));
    if rel_std > 0.05
        warning('A0 values from different components vary by more than 5%% (std/mean = %.1f%%)', ...
            rel_std*100);
        fprintf('  A0 from each component: [%.2f, %.2f, %.2f] MHz\n', A0_components);
        fprintf('  Mean A0: %.2f MHz\n', A0_MHz);
    end
end
