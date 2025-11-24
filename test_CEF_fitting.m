%% Test script for updated CEF_fitting.m
% This script tests the CEF fitting with the target spectrum from the
% effective spin-1/2 model

clearvars
close all

fprintf('=================================================================\n');
fprintf('Testing CEF_fitting.m\n');
fprintf('=================================================================\n\n');

%% Setup constants
const.hbar = 1.05457E-34;
const.muB = 9.274e-24;
const.muN = 5.05078e-27;
const.kB = 1.3806e-23; % [J/K]
const.J2meV = 6.24151e+21;
const.Gh2mV = const.hbar * 2*pi * 10^9 * const.J2meV;

%% Target parameters from effective model
targetG = [8.3, 8.3, 1.26];
targetAeff_MHz = [-871.1, -871.1, -130.3];

fprintf('Target parameters from effective spin-1/2 model:\n');
fprintf('  g_eff = [%.3f, %.3f, %.3f]\n', targetG);
fprintf('  A_eff = [%.1f, %.1f, %.1f] MHz\n\n', targetAeff_MHz);

%% Generate target 16-level spectrum from effective model
fprintf('Generating target 16-level spectrum...\n');

I = 7/2;
[~,~,~,~,~,~,Jx_eff,Jy_eff,Jz_eff,Ix_eff,Iy_eff,Iz_eff] = spin_operators(0.5, I);

% Build effective hyperfine Hamiltonian
A_eff_meV = (targetAeff_MHz/1000) * const.Gh2mV;
Hhf_eff = A_eff_meV(1)*Jx_eff*Ix_eff + ...
          A_eff_meV(2)*Jy_eff*Iy_eff + ...
          A_eff_meV(3)*Jz_eff*Iz_eff;

% Get eigenvalues (these are the target to fit)
[~, E_eff] = eig(Hhf_eff);
eigenE_target = sort(real(diag(E_eff)));

fprintf('  Energy range: [%.6f, %.6f] meV\n', min(eigenE_target), max(eigenE_target));
fprintf('  Energy spread: %.6f meV\n\n', max(eigenE_target) - min(eigenE_target));

%% Test current CEF parameters (before optimization)
fprintf('=== Evaluating CURRENT CEF parameters (baseline) ===\n');

% Current CEF parameters from MF_Er_CaWO4_v1b.m
B_current_cm = [-330.74 2212.59 196.501 5.57413 -8.07197 -217.553 81.2414];
fprintf('Current CEF parameters:\n');
fprintf('  B = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', B_current_cm);

L = 6; S = 3/2;
J = L + S;
gJ = gLande(L, S);

% Build CEF Hamiltonian
B_current_meV = B_current_cm * 0.123983;
[~,~,~,~,~,~,Jx,Jy,Jz,~,~,~] = spin_operators(J, 0);
Hcf_current = cf(J, B_current_meV, 0);

% Use SW_proj (consistent with MF_Er_CaWO4_v1b.m)
ion_sw.J = [J];
ion_sw.Jx = Jx;
ion_sw.Jy = Jy;
ion_sw.Jz = Jz;
ion_sw.Hcf = Hcf_current;
ion_sw.h4 = 0;
ion_sw.gLande = [gJ];
ion_sw.idx = 1;

params_sw.temp = 0.1;
params_sw.field = [0; 0; 0];

[~, ~, basis_sw, ~, ~, ~] = SW_proj(const, ion_sw, params_sw);

% Project operators
wav0_sw = basis_sw(:,:,1,1,1);
jx_sw = real(wav0_sw' * Jx * wav0_sw);
jy_sw = real(wav0_sw' * Jy * wav0_sw);
jz_sw = real(wav0_sw' * Jz * wav0_sw);

% Calculate g-factors
g_eff_current = zeros(1,3);
g_eff_current(1) = 2 * gJ * abs(jx_sw(1,2));
g_eff_current(2) = 2 * gJ * abs(jy_sw(1,2));
g_eff_current(3) = 2 * gJ * abs(jz_sw(1,1));

fprintf('\nEffective g-factors:\n');
fprintf('  g_eff = [%.6f, %.6f, %.6f]\n', g_eff_current);
fprintf('  Target = [%.6f, %.6f, %.6f]\n', targetG);
fprintf('  Relative error: [%.2e, %.2e, %.2e]\n', ...
    (g_eff_current - targetG)./targetG);

% Calculate 16-level spectrum with corrected A0
Iz_nuc = diag(I:-1:-I);
Ip_nuc = diag(sqrt((I - (I-1:-1:-I)) .* (I+1 + (I-1:-1:-I))), 1);
Ix_nuc = (Ip_nuc + Ip_nuc')/2;
Iy_nuc = (Ip_nuc - Ip_nuc')/2i;

A0_corrected_MHz = -250.6;
A0_meV = (A0_corrected_MHz/1000) * const.Gh2mV;

Hhf_cef_current = A0_meV * (kron(jx_sw, Ix_nuc) + ...
                             kron(jy_sw, Iy_nuc) + ...
                             kron(jz_sw, Iz_nuc));
[~, E_cef_current] = eig(Hhf_cef_current);
E_cef_current = sort(real(diag(E_cef_current)));

fprintf('\n16-level spectrum comparison:\n');
rms_error_current = sqrt(mean((E_cef_current - eigenE_target).^2));
max_error_current = max(abs(E_cef_current - eigenE_target));
fprintf('  RMS error: %.6f meV (%.3f GHz)\n', rms_error_current, rms_error_current / const.Gh2mV * 1000);
fprintf('  Max error: %.6f meV (%.3f GHz)\n', max_error_current, max_error_current / const.Gh2mV * 1000);

%% Run CEF_fitting
fprintf('\n=== Running CEF_fitting optimization ===\n');
fprintf('This may take a few minutes...\n\n');

tic;
[B_fit_cm, g_eff_fit, stats] = CEF_fitting(eigenE_target);
elapsed = toc;

fprintf('\n=== Optimization Results ===\n');
fprintf('Elapsed time: %.1f seconds\n', elapsed);
fprintf('Exit flag: %d\n', stats.exitflag);
fprintf('Final resnorm: %.6g\n', stats.resnorm);
fprintf('Function evaluations: %d\n', stats.funcCount);

fprintf('\nFitted CEF parameters (cm^-1):\n');
fprintf('  B_fit = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', B_fit_cm);

fprintf('\nParameter changes from baseline:\n');
dB = B_fit_cm - B_current_cm;
fprintf('  ΔB = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', dB);
fprintf('  Relative change (rms): %.2f%%\n', 100 * sqrt(mean((dB./B_current_cm).^2)));

fprintf('\nFitted effective g-factors:\n');
fprintf('  g_fit = [%.6f, %.6f, %.6f]\n', g_eff_fit);
fprintf('  Target = [%.6f, %.6f, %.6f]\n', targetG);
fprintf('  Relative error: [%.2e, %.2e, %.2e]\n', ...
    (g_eff_fit - targetG)./targetG);

%% Evaluate fitted parameters
fprintf('\n=== Evaluating FITTED CEF parameters ===\n');

B_fit_meV = B_fit_cm * 0.123983;
Hcf_fit = cf(J, B_fit_meV, 0);

% Use SW_proj
ion_sw.Hcf = Hcf_fit;
[~, ~, basis_fit, ~, ~, ~] = SW_proj(const, ion_sw, params_sw);

% Project operators
wav0_fit = basis_fit(:,:,1,1,1);
jx_fit = real(wav0_fit' * Jx * wav0_fit);
jy_fit = real(wav0_fit' * Jy * wav0_fit);
jz_fit = real(wav0_fit' * Jz * wav0_fit);

% Calculate 16-level spectrum
Hhf_cef_fit = A0_meV * (kron(jx_fit, Ix_nuc) + ...
                         kron(jy_fit, Iy_nuc) + ...
                         kron(jz_fit, Iz_nuc));
[~, E_cef_fit] = eig(Hhf_cef_fit);
E_cef_fit = sort(real(diag(E_cef_fit)));

fprintf('16-level spectrum comparison:\n');
rms_error_fit = sqrt(mean((E_cef_fit - eigenE_target).^2));
max_error_fit = max(abs(E_cef_fit - eigenE_target));
fprintf('  RMS error: %.6f meV (%.3f GHz)\n', rms_error_fit, rms_error_fit / const.Gh2mV * 1000);
fprintf('  Max error: %.6f meV (%.3f GHz)\n', max_error_fit, max_error_fit / const.Gh2mV * 1000);

fprintf('\nImprovement:\n');
fprintf('  RMS error reduction: %.1f%%\n', 100 * (1 - rms_error_fit/rms_error_current));
fprintf('  Max error reduction: %.1f%%\n', 100 * (1 - max_error_fit/max_error_current));

%% Plot comparison
figure('Position', [100 100 1400 500]);

subplot(1,3,1)
hold on
plot(1:16, eigenE_target*1000, 'ko-', 'LineWidth', 2.5, 'MarkerSize', 10, 'DisplayName', 'Target (eff. model)');
plot(1:16, E_cef_current*1000, 'bs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Current CEF');
plot(1:16, E_cef_fit*1000, 'r^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Fitted CEF');
xlabel('Level index')
ylabel('Energy (μeV)')
title('16-level Hyperfine Spectrum')
legend('Location', 'best')
grid on
set(gca, 'FontSize', 12)

subplot(1,3,2)
hold on
plot(1:16, (E_cef_current - eigenE_target)*1000, 'bs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Current CEF');
plot(1:16, (E_cef_fit - eigenE_target)*1000, 'r^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Fitted CEF');
plot([0 17], [0 0], 'k--', 'LineWidth', 1)
xlabel('Level index')
ylabel('Error (μeV)')
title('Spectrum Error from Target')
legend('Location', 'best')
grid on
set(gca, 'FontSize', 12)

subplot(1,3,3)
categories = categorical({'x', 'y', 'z'});
bar(categories, [targetG; g_eff_current; g_eff_fit]')
legend({'Target', 'Current', 'Fitted'}, 'Location', 'best')
ylabel('g-factor')
title('Effective g-factors')
grid on
set(gca, 'FontSize', 12)

sgtitle('CEF Fitting Results', 'FontSize', 14, 'FontWeight', 'bold')

fprintf('\n=================================================================\n');
fprintf('Test complete. See figure for visual comparison.\n');
fprintf('=================================================================\n');
