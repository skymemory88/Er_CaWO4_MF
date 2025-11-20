%% Diagnostic script to identify discrepancies between CEF and effective models
% This script compares the eigenstates and excitation spectrum from:
% 1. CEF model with fitted parameters
% 2. Effective spin-1/2 model
clearvars

fprintf('===================================================================\n');
fprintf('CEF vs Effective Model Diagnostic\n');
fprintf('===================================================================\n\n');

%% Constants
const.hbar = 1.05457E-34;
const.muB = 9.274e-24;
const.muN = 5.05078e-27;
const.kB = 1.3806e-23;
const.mu0 = 4e-7 * pi;
const.kB_meV = 8.61733e-2;
const.J2meV = 6.24151e+21;
const.Gh2mV = const.hbar * 2*pi * 10^9 * const.J2meV;

%% Ion parameters
L = 6; S = 3/2; % Er3+
J = L + S;
I = 7/2; % Nuclear spin

gJ = gLande(L, S);
fprintf('Lande g-factor: gJ = %.6f\n', gJ);

% Target effective parameters (from effective model)
g_eff_target = [8.3, 8.3, 1.26];
A_eff_target_MHz = [-871.1, -871.1, -130.3];

fprintf('\nTarget effective parameters:\n');
fprintf('  g_eff = [%.3f, %.3f, %.3f]\n', g_eff_target);
fprintf('  A_eff = [%.1f, %.1f, %.1f] MHz\n', A_eff_target_MHz);

%% CEF parameters (from MF_Er_CaWO4_v1b.m line 116)
B_cm = [-330.74 2212.59 196.501 5.57413 -8.07197 -217.553 81.2414];
B_meV = B_cm * 0.123983; % cm^-1 to meV

fprintf('\nCEF parameters (cm^-1):\n');
fprintf('  B = [');
for i = 1:length(B_cm)
    fprintf('%.3f', B_cm(i));
    if i < length(B_cm), fprintf(', '); end
end
fprintf(']\n');

%% Build CEF Hamiltonian
Hcf = cf(J, B_meV, 0);

% Get J operators
[~,~,~,~,~,~,Jx,Jy,Jz,~,~,~] = spin_operators(J, 0);

%% Method 1: projectDoublet (used in CEF_fitting.m)
fprintf('\n--- Method 1: projectDoublet (zero-field) ---\n');
[proj1, Jproj1] = projectDoublet(Hcf, Jx, Jy, Jz);

fprintf('Ground doublet energies: [%.6f, %.6f] meV\n', proj1.energies(1), proj1.energies(2));
fprintf('Energy gap: %.6f meV\n', proj1.energies(2) - proj1.energies(1));

fprintf('\nProjected J operators (Method 1):\n');
fprintf('Jx:\n'); disp(Jproj1.Jx);
fprintf('Jy:\n'); disp(Jproj1.Jy);
fprintf('Jz:\n'); disp(Jproj1.Jz);

% Calculate effective g-factors
g_eff_1 = zeros(1,3);
g_eff_1(1) = 2 * gJ * abs(Jproj1.Jx(1,2));
g_eff_1(2) = 2 * gJ * abs(Jproj1.Jy(1,2));
g_eff_1(3) = 2 * gJ * abs(Jproj1.Jz(1,1));

fprintf('\nEffective g-factors (Method 1):\n');
fprintf('  g_eff = [%.6f, %.6f, %.6f]\n', g_eff_1);
fprintf('  Relative error: [%.2e, %.2e, %.2e]\n', ...
    (g_eff_1 - g_eff_target)./g_eff_target);

%% Method 2: SW_proj (used in MF_Er_CaWO4_v1b.m)
fprintf('\n--- Method 2: SW_proj (Schrieffer-Wolff) ---\n');

% Set up ion structure for SW_proj
ion_sw.J = [J];
ion_sw.Jx = Jx;
ion_sw.Jy = Jy;
ion_sw.Jz = Jz;
ion_sw.Hcf = Hcf;
ion_sw.h4 = 0;
ion_sw.gLande = [gJ];
ion_sw.idx = 1;

% Test at zero field
params_sw.temp = 0.1;
params_sw.field = [0; 0; 0];

[en_sw, ham_E_sw, basis_sw, ~, ~, ~] = SW_proj(const, ion_sw, params_sw);

fprintf('Ground doublet energies: [%.6f, %.6f] meV\n', en_sw(1,1,1), en_sw(2,1,1));
fprintf('Energy gap: %.6f meV\n', en_sw(2,1,1) - en_sw(1,1,1));

% Project J operators using SW basis
wav0_sw = basis_sw(:,:,1,1,1);
jx_sw = real(wav0_sw' * Jx * wav0_sw);
jy_sw = real(wav0_sw' * Jy * wav0_sw);
jz_sw = real(wav0_sw' * Jz * wav0_sw);

fprintf('\nProjected J operators (Method 2):\n');
fprintf('Jx:\n'); disp(jx_sw);
fprintf('Jy:\n'); disp(jy_sw);
fprintf('Jz:\n'); disp(jz_sw);

% Calculate effective g-factors
g_eff_2 = zeros(1,3);
g_eff_2(1) = 2 * gJ * abs(jx_sw(1,2));
g_eff_2(2) = 2 * gJ * abs(jy_sw(1,2));
g_eff_2(3) = 2 * gJ * abs(jz_sw(1,1));

fprintf('\nEffective g-factors (Method 2):\n');
fprintf('  g_eff = [%.6f, %.6f, %.6f]\n', g_eff_2);
fprintf('  Relative error: [%.2e, %.2e, %.2e]\n', ...
    (g_eff_2 - g_eff_target)./g_eff_target);

fprintf('\nDifference between methods:\n');
fprintf('  Δg_eff = [%.2e, %.2e, %.2e]\n', g_eff_2 - g_eff_1);

%% Calculate required isotropic A0 for effective hyperfine
fprintf('\n--- Hyperfine Analysis ---\n');

% From effective model: A_eff = A0 * <J_eff>
% Where <J_eff_i> = |<0|J_i|1>| for i=x,y and <0|J_z|0> for i=z
% So: A0 = A_eff * (2*gJ) / g_eff

A0_from_target = zeros(1,3);
A0_from_target(1) = A_eff_target_MHz(1) * (2*gJ) / g_eff_target(1);
A0_from_target(2) = A_eff_target_MHz(2) * (2*gJ) / g_eff_target(2);
A0_from_target(3) = A_eff_target_MHz(3) * (2*gJ) / g_eff_target(3);

fprintf('Required A0 from target A_eff:\n');
fprintf('  A0_x = %.2f MHz\n', A0_from_target(1));
fprintf('  A0_y = %.2f MHz\n', A0_from_target(2));
fprintf('  A0_z = %.2f MHz\n', A0_from_target(3));
fprintf('  Mean A0 = %.2f MHz\n', mean(A0_from_target));

A0_current_MHz = -125.9; % From line 137 of MF_Er_CaWO4_v1b.m
fprintf('\nCurrent A0 in MF_Er_CaWO4_v1b.m: %.2f MHz\n', A0_current_MHz);
fprintf('Correction factor needed: %.3f\n', mean(A0_from_target) / A0_current_MHz);

% Calculate what A_eff we get with current parameters
A_eff_current_MHz = A0_current_MHz * g_eff_1 / (2*gJ);
fprintf('\nEffective A with current parameters:\n');
fprintf('  A_eff = [%.1f, %.1f, %.1f] MHz\n', A_eff_current_MHz);
fprintf('  Target A_eff = [%.1f, %.1f, %.1f] MHz\n', A_eff_target_MHz);
fprintf('  Relative error: [%.2e, %.2e, %.2e]\n', ...
    (A_eff_current_MHz - A_eff_target_MHz)./A_eff_target_MHz);

%% Compare 16-level spectra at zero field
fprintf('\n--- 16-level Spectrum Comparison ---\n');

% Create nuclear spin operators
Iz_nuc = diag(I:-1:-I);
Ip_nuc = diag(sqrt((I - (I-1:-1:-I)) .* (I+1 + (I-1:-1:-I))), 1);
Im_nuc = Ip_nuc';
Ix_nuc = (Ip_nuc + Im_nuc)/2;
Iy_nuc = (Ip_nuc - Im_nuc)/2i;

% Method 1: CEF model (current implementation)
fprintf('\nCEF model (current A0 = %.2f MHz):\n', A0_current_MHz);
A0_current_meV = (A0_current_MHz/1000) * const.Gh2mV;

Hhf_cef = A0_current_meV * (kron(Jproj1.Jx, Ix_nuc) + ...
                              kron(Jproj1.Jy, Iy_nuc) + ...
                              kron(Jproj1.Jz, Iz_nuc));
[~, E_cef] = eig(Hhf_cef);
E_cef = sort(real(diag(E_cef)));
E_cef = E_cef - mean(E_cef); % Center

fprintf('  Energy range: [%.4f, %.4f] meV\n', min(E_cef), max(E_cef));
fprintf('  Energy spread: %.4f meV\n', max(E_cef) - min(E_cef));

% Method 2: Effective model
fprintf('\nEffective spin-1/2 model:\n');
[~,~,~,~,~,~,Jx_eff,Jy_eff,Jz_eff,Ix_eff,Iy_eff,Iz_eff] = spin_operators(0.5, I);

A_eff_meV = (A_eff_target_MHz/1000) * const.Gh2mV;
Hhf_eff = A_eff_meV(1)*Jx_eff*Ix_eff + ...
          A_eff_meV(2)*Jy_eff*Iy_eff + ...
          A_eff_meV(3)*Jz_eff*Iz_eff;
[~, E_eff] = eig(Hhf_eff);
E_eff = sort(real(diag(E_eff)));
E_eff = E_eff - mean(E_eff); % Center

fprintf('  Energy range: [%.4f, %.4f] meV\n', min(E_eff), max(E_eff));
fprintf('  Energy spread: %.4f meV\n', max(E_eff) - min(E_eff));

% Method 3: CEF model with corrected A0
fprintf('\nCEF model (corrected A0 = %.2f MHz):\n', mean(A0_from_target));
A0_corrected_meV = (mean(A0_from_target)/1000) * const.Gh2mV;

Hhf_cef_corr = A0_corrected_meV * (kron(Jproj1.Jx, Ix_nuc) + ...
                                    kron(Jproj1.Jy, Iy_nuc) + ...
                                    kron(Jproj1.Jz, Iz_nuc));
[~, E_cef_corr] = eig(Hhf_cef_corr);
E_cef_corr = sort(real(diag(E_cef_corr)));
E_cef_corr = E_cef_corr - mean(E_cef_corr); % Center

fprintf('  Energy range: [%.4f, %.4f] meV\n', min(E_cef_corr), max(E_cef_corr));
fprintf('  Energy spread: %.4f meV\n', max(E_cef_corr) - min(E_cef_corr));

% Compare spectra
fprintf('\nSpectrum comparison:\n');
fprintf('  RMS difference (CEF current vs Effective): %.4f meV\n', ...
    sqrt(mean((E_cef - E_eff).^2)));
fprintf('  RMS difference (CEF corrected vs Effective): %.4f meV\n', ...
    sqrt(mean((E_cef_corr - E_eff).^2)));

%% Plot comparison
figure('Position', [100 100 1200 800]);

subplot(2,2,1)
hold on
plot(1:16, E_eff*1000, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Effective model');
plot(1:16, E_cef*1000, 'rx-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'CEF (current A0)');
plot(1:16, E_cef_corr*1000, 'gs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'CEF (corrected A0)');
xlabel('Level index')
ylabel('Energy (μeV)')
title('16-level Hyperfine Spectrum at Zero Field')
legend('Location', 'best')
grid on
set(gca, 'FontSize', 12)

subplot(2,2,2)
hold on
plot(1:16, (E_cef - E_eff)*1000, 'rx-', 'LineWidth', 2, 'DisplayName', 'Current A0');
plot(1:16, (E_cef_corr - E_eff)*1000, 'gs-', 'LineWidth', 2, 'DisplayName', 'Corrected A0');
plot([0 17], [0 0], 'k--', 'LineWidth', 1)
xlabel('Level index')
ylabel('Energy difference (μeV)')
title('Difference from Effective Model')
legend('Location', 'best')
grid on
set(gca, 'FontSize', 12)

subplot(2,2,3)
bar(1:3, [g_eff_target; g_eff_1; g_eff_2]')
legend({'Target', 'projectDoublet', 'SW\_proj'}, 'Location', 'best')
xlabel('Component (x, y, z)')
ylabel('g-factor')
title('Effective g-factor Comparison')
set(gca, 'XTickLabel', {'x', 'y', 'z'}, 'FontSize', 12)
grid on

subplot(2,2,4)
bar(1:3, [A_eff_target_MHz; A_eff_current_MHz; ...
          mean(A0_from_target) * g_eff_1 / (2*gJ)]')
legend({'Target', 'Current (A0=-125.9 MHz)', 'Corrected'}, 'Location', 'best')
xlabel('Component (x, y, z)')
ylabel('A_{eff} (MHz)')
title('Effective Hyperfine Comparison')
set(gca, 'XTickLabel', {'x', 'y', 'z'}, 'FontSize', 12)
grid on

sgtitle('CEF vs Effective Model Diagnostic', 'FontSize', 14, 'FontWeight', 'bold')

fprintf('\n===================================================================\n');
fprintf('Diagnostic complete. See figure for visual comparison.\n');
fprintf('===================================================================\n');
