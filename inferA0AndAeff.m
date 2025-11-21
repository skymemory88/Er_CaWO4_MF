function [A0_axes_MHz, A0_rec_MHz, Aeff_fromA0_MHz] = inferA0AndAeff(g_eff, gJ, targetAeff_MHz)
%INFERA0ANDAEFF Infer isotropic A0 from effective g-tensor and target A_eff
%   [A0_axes_MHz, A0_rec_MHz, Aeff_fromA0_MHz] = inferA0AndAeff(g_eff, gJ, targetAeff_MHz)
%
%   Calculates the isotropic hyperfine constant A0 from the effective
%   anisotropic hyperfine constants and the fitted g-tensor.
%
%   Inputs:
%     g_eff          : 1x3 vector of fitted effective g-factors [gx, gy, gz]
%     gJ             : Lande g-factor (scalar)
%     targetAeff_MHz : 1x3 vector of target effective A values [Ax, Ay, Az] in MHz
%
%   Outputs:
%     A0_axes_MHz    : 1x3 vector of A0 inferred from each axis
%     A0_rec_MHz     : Recommended A0 (mean of the three components)
%     Aeff_fromA0_MHz: 1x3 vector of A_eff recalculated from A0_rec and g_eff
%
%   The relationship is: A_eff_i = A0 * g_eff_i / (2*gJ)
%   Therefore: A0 = A_eff_i * (2*gJ) / g_eff_i

    % Calculate A0 from each axis
    A0_axes_MHz = targetAeff_MHz .* (2*gJ) ./ g_eff;

    % Recommended A0 is the mean
    A0_rec_MHz = mean(A0_axes_MHz);

    % Calculate what A_eff would be using the recommended A0 and fitted g
    Aeff_fromA0_MHz = A0_rec_MHz * g_eff / (2*gJ);
end
