%% Analyze CEF_fitting.m behavior and identify issues
% This script tests the CEF_fitting optimization to understand:
% 1. Is the objective function correctly formulated?
% 2. Are we using consistent projection methods?
% 3. What does the optimization landscape look like?

clearvars

fprintf('=================================================================\n');
fprintf('CEF_fitting.m Analysis\n');
fprintf('=================================================================\n\n');

%% Setup constants and parameters (matching MF_Er_CaWO4_v1b.m)
const.hbar = 1.05457E-34;
const.muB = 9.274e-24;
const.muN = 5.05078e-27;
const.J2meV = 6.24151e+21;
const.Gh2mV = const.hbar * 2*pi * 10^9 * const.J2meV;
const.kB = 1.3806e-23;

L = 6; S = 3/2;
J = L + S;
I = 7/2;
gJ = gLande(L, S);

fprintf('Ion parameters:\n');
fprintf('  L = %g, S = %g, J = %g\n', L, S, J);
fprintf('  gJ = %.6f\n', gJ);
fprintf('  I = %g\n\n', I);

%% Target parameters from effective model
targetG = [8.3, 8.3, 1.26];
targetAeff_MHz = [-871.1, -871.1, -130.3];

fprintf('Target effective parameters:\n');
fprintf('  g_eff = [%.3f, %.3f, %.3f]\n', targetG);
fprintf('  A_eff = [%.1f, %.1f, %.1f] MHz\n\n', targetAeff_MHz);

%% Generate target 16-level spectrum from effective model
fprintf('Generating target 16-level spectrum from effective model...\n');

% Effective model operators
[~,~,~,~,~,~,Jx_eff,Jy_eff,Jz_eff,Ix_eff,Iy_eff,Iz_eff] = spin_operators(0.5, I);

% Build effective hyperfine Hamiltonian
A_eff_meV = (targetAeff_MHz/1000) * const.Gh2mV;
Hhf_eff = A_eff_meV(1)*Jx_eff*Ix_eff + ...
          A_eff_meV(2)*Jy_eff*Iy_eff + ...
          A_eff_meV(3)*Jz_eff*Iz_eff;

% Get eigenvalues
[~, E_eff] = eig(Hhf_eff);
eigenE_target = sort(real(diag(E_eff)));

fprintf('  Energy range: [%.6f, %.6f] meV\n', min(eigenE_target), max(eigenE_target));
fprintf('  Energy spread: %.6f meV\n\n', max(eigenE_target) - min(eigenE_target));

%% Current CEF parameters from MF_Er_CaWO4_v1b.m
B_current_cm = [-330.74 2212.59 196.501 5.57413 -8.07197 -217.553 81.2414];
B_current_meV = B_current_cm * 0.123983;

fprintf('Current CEF parameters (cm^-1):\n');
fprintf('  B = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n\n', B_current_cm);

%% Test 1: Evaluate current parameters with projectDoublet
fprintf('--- Test 1: Current CEF with projectDoublet (CEF_fitting method) ---\n');

[~,~,~,~,~,~,Jx,Jy,Jz,~,~,~] = spin_operators(J, 0);
Hcf_current = cf(J, B_current_meV, 0);
[proj1, Jproj1] = projectDoublet(Hcf_current, Jx, Jy, Jz);

% Calculate g-factors
g_eff_1 = zeros(1,3);
g_eff_1(1) = 2 * gJ * abs(Jproj1.Jx(1,2));
g_eff_1(2) = 2 * gJ * abs(Jproj1.Jy(1,2));
g_eff_1(3) = 2 * gJ * abs(Jproj1.Jz(1,1));

fprintf('Effective g-factors: [%.6f, %.6f, %.6f]\n', g_eff_1);
fprintf('Relative error:      [%.2e, %.2e, %.2e]\n', (g_eff_1 - targetG)./targetG);

% Calculate 16-level spectrum
Iz_nuc = diag(I:-1:-I);
Ip_nuc = diag(sqrt((I - (I-1:-1:-I)) .* (I+1 + (I-1:-1:-I))), 1);
Ix_nuc = (Ip_nuc + Ip_nuc')/2;
Iy_nuc = (Ip_nuc - Ip_nuc')/2i;

% Use corrected A0
A0_corrected_MHz = -250.6;
A0_meV = (A0_corrected_MHz/1000) * const.Gh2mV;

Hhf_cef1 = A0_meV * (kron(Jproj1.Jx, Ix_nuc) + ...
                      kron(Jproj1.Jy, Iy_nuc) + ...
                      kron(Jproj1.Jz, Iz_nuc));
[~, E_cef1] = eig(Hhf_cef1);
E_cef1 = sort(real(diag(E_cef1)));

fprintf('\nSpectrum comparison (projectDoublet):\n');
fprintf('  RMS error: %.6f meV\n', sqrt(mean((E_cef1 - eigenE_target).^2)));
fprintf('  Max error: %.6f meV\n', max(abs(E_cef1 - eigenE_target)));

%% Test 2: Evaluate current parameters with SW_proj
fprintf('\n--- Test 2: Current CEF with SW_proj (MF_Er_CaWO4_v1b method) ---\n');

% Setup ion structure
ion_sw.J = [J];
ion_sw.Jx = Jx;
ion_sw.Jy = Jy;
ion_sw.Jz = Jz;
ion_sw.Hcf = Hcf_current;
ion_sw.h4 = 0;
ion_sw.gLande = [gJ];
ion_sw.idx = 1;

% Zero field
params_sw.temp = 0.1;
params_sw.field = [0; 0; 0];

[~, ~, basis_sw, ~, ~, ~] = SW_proj(const, ion_sw, params_sw);

% Project operators
wav0_sw = basis_sw(:,:,1,1,1);
jx_sw = real(wav0_sw' * Jx * wav0_sw);
jy_sw = real(wav0_sw' * Jy * wav0_sw);
jz_sw = real(wav0_sw' * Jz * wav0_sw);

% Calculate g-factors
g_eff_2 = zeros(1,3);
g_eff_2(1) = 2 * gJ * abs(jx_sw(1,2));
g_eff_2(2) = 2 * gJ * abs(jy_sw(1,2));
g_eff_2(3) = 2 * gJ * abs(jz_sw(1,1));

fprintf('Effective g-factors: [%.6f, %.6f, %.6f]\n', g_eff_2);
fprintf('Relative error:      [%.2e, %.2e, %.2e]\n', (g_eff_2 - targetG)./targetG);

% Calculate 16-level spectrum
Hhf_cef2 = A0_meV * (kron(jx_sw, Ix_nuc) + ...
                      kron(jy_sw, Iy_nuc) + ...
                      kron(jz_sw, Iz_nuc));
[~, E_cef2] = eig(Hhf_cef2);
E_cef2 = sort(real(diag(E_cef2)));

fprintf('\nSpectrum comparison (SW_proj):\n');
fprintf('  RMS error: %.6f meV\n', sqrt(mean((E_cef2 - eigenE_target).^2)));
fprintf('  Max error: %.6f meV\n', max(abs(E_cef2 - eigenE_target)));

%% Test 3: Compare projection methods
fprintf('\n--- Test 3: Difference between projection methods ---\n');
fprintf('Δg_eff = [%.2e, %.2e, %.2e]\n', g_eff_2 - g_eff_1);
fprintf('Spectrum RMS difference: %.6f meV\n', sqrt(mean((E_cef2 - E_cef1).^2)));

%% Test 4: Evaluate objective function value
fprintf('\n--- Test 4: CEF_fitting objective function evaluation ---\n');

% Prepare inputs as in CEF_fitting
eigenE_normalized = eigenE_target - min(eigenE_target);
[~, ~, ~, Ix, Iy, Iz, ~, ~, ~, ~, ~, ~] = spin_operators(J, I);

% Weights from CEF_fitting
gTol = 1e-3;
wG = 1.0;
wE = 1.0;

% Calculate residuals for current parameters
% Using projectDoublet (as in CEF_fitting)
gres = (g_eff_1(:) - targetG(:)) ./ targetG(:);
gres = gres / gTol;

% Energy residual
Hhf_unit = kron(Jproj1.Jx, Ix) + kron(Jproj1.Jy, Iy) + kron(Jproj1.Jz, Iz);
evals = eig((Hhf_unit + Hhf_unit')/2);
evals = sort(real(evals));
evals = evals - mean(evals);

% Optimal A0 fit
A0_opt = (sum(evals .* eigenE_normalized(:))) / (sum(evals .* evals) + eps);
epred = A0_opt * evals;
scale = sqrt(mean(eigenE_normalized.^2)) + eps;
eres = (epred - eigenE_normalized(:)) / scale;

% Total residual
res_total = [wG * gres; wE * eres];
resnorm = sum(res_total.^2);

fprintf('Current CEF parameters:\n');
fprintf('  g-residual norm:      %.6f\n', norm(wG * gres));
fprintf('  Energy residual norm: %.6f\n', norm(wE * eres));
fprintf('  Total resnorm:        %.6f\n', resnorm);
fprintf('  Optimal A0 (in eigenE units): %.6f\n', A0_opt);

%% Test 5: Check excited state structure
fprintf('\n--- Test 5: CEF excited state structure ---\n');

[eigVec_cef, eigVal_cef] = eig(Hcf_current);
E_cef_levels = sort(real(diag(eigVal_cef)));
E_cef_levels = E_cef_levels - E_cef_levels(1); % Relative to ground state

fprintf('First 8 CEF levels (meV):\n');
for i = 1:min(8, length(E_cef_levels))
    fprintf('  Level %d: %.4f meV\n', i, E_cef_levels(i));
end

% Energy gaps
fprintf('\nKey energy gaps:\n');
fprintf('  Ground doublet splitting: %.4f meV (%.2f GHz)\n', ...
    E_cef_levels(2) - E_cef_levels(1), ...
    (E_cef_levels(2) - E_cef_levels(1)) / const.Gh2mV * 1000);
fprintf('  Gap to 1st excited doublet: %.4f meV (%.2f GHz)\n', ...
    E_cef_levels(3) - E_cef_levels(1), ...
    (E_cef_levels(3) - E_cef_levels(1)) / const.Gh2mV * 1000);

%% Visualize spectrum comparison
figure('Position', [100 100 1400 600]);

subplot(1,3,1)
hold on
plot(1:16, eigenE_target*1000, 'ko-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Target (effective)');
plot(1:16, E_cef1*1000, 'bx-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'CEF (projectDoublet)');
plot(1:16, E_cef2*1000, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'CEF (SW\_proj)');
xlabel('Level index')
ylabel('Energy (μeV)')
title('16-level Hyperfine Spectrum')
legend('Location', 'best')
grid on
set(gca, 'FontSize', 12)

subplot(1,3,2)
hold on
plot(1:16, (E_cef1 - eigenE_target)*1000, 'bx-', 'LineWidth', 2, 'DisplayName', 'projectDoublet');
plot(1:16, (E_cef2 - eigenE_target)*1000, 'rs-', 'LineWidth', 2, 'DisplayName', 'SW\_proj');
plot([0 17], [0 0], 'k--', 'LineWidth', 1)
xlabel('Level index')
ylabel('Error (μeV)')
title('Spectrum Error from Target')
legend('Location', 'best')
grid on
set(gca, 'FontSize', 12)

subplot(1,3,3)
bar([1:3; 4:6; 7:9]', [targetG; g_eff_1; g_eff_2]')
legend({'g_x', 'g_y', 'g_z'}, 'Location', 'best')
set(gca, 'XTickLabel', {'Target', 'projectDoublet', 'SW\_proj'})
ylabel('g-factor')
title('Effective g-factor Comparison')
grid on
set(gca, 'FontSize', 12)

sgtitle('CEF Parameters Evaluation', 'FontSize', 14, 'FontWeight', 'bold')

fprintf('\n=================================================================\n');
fprintf('Analysis complete.\n');
fprintf('=================================================================\n');
