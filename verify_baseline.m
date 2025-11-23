%% Quick verification of baseline CEF parameters
% This script quickly checks if the baseline parameters give reasonable g-factors
% Run this BEFORE running the full optimization to verify starting point

clearvars

fprintf('======================================\n');
fprintf('Baseline CEF Parameters Verification\n');
fprintf('======================================\n\n');

%% Constants
const.hbar = 1.05457E-34;
const.muB = 9.274e-24;
const.muN = 5.05078e-27;
const.J2meV = 6.24151e+21;
const.Gh2mV = const.hbar * 2*pi * 10^9 * const.J2meV;
const.kB = 1.3806e-23;

%% Ion parameters
L = 6; S = 3/2;
J = L + S;
I = 7/2;
gJ = gLande(L, S);

fprintf('Ion: Er3+ (L=%g, S=%g, J=%g, gJ=%.6f)\n\n', L, S, J, gJ);

%% Target
targetG = [8.3, 8.3, 1.26];
fprintf('Target g_eff = [%.3f, %.3f, %.3f]\n\n', targetG);

%% Baseline from MF_Er_CaWO4_v1b.m line 115
B_baseline_cm = [-979.402 1898.05 346.52 -272.893 -7.86205 -187.516 58.9021];
fprintf('Baseline CEF parameters (line 115):\n');
fprintf('  B = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n\n', B_baseline_cm);

%% Test with projectDoublet (fast)
fprintf('Testing with projectDoublet (zero-field)...\n');
B_meV = B_baseline_cm * 0.123983;
[~,~,~,~,~,~,Jx,Jy,Jz,~,~,~] = spin_operators(J, 0);
Hcf = cf(J, B_meV, 0);
[proj, Jproj] = projectDoublet(Hcf, Jx, Jy, Jz);

g_eff_proj = zeros(1,3);
g_eff_proj(1) = 2 * gJ * abs(Jproj.Jx(1,2));
g_eff_proj(2) = 2 * gJ * abs(Jproj.Jy(1,2));
g_eff_proj(3) = 2 * gJ * abs(Jproj.Jz(1,1));

fprintf('  g_eff = [%.6f, %.6f, %.6f]\n', g_eff_proj);
fprintf('  Error = [%.2e, %.2e, %.2e]\n', (g_eff_proj - targetG)./targetG);

%% Test with SW_proj (consistent with MF script)
fprintf('\nTesting with SW_proj (zero-field)...\n');

ion_sw.J = [J];
ion_sw.Jx = Jx;
ion_sw.Jy = Jy;
ion_sw.Jz = Jz;
ion_sw.Hcf = Hcf;
ion_sw.h4 = 0;
ion_sw.gLande = [gJ];
ion_sw.idx = 1;

params_sw.temp = 0.1;
params_sw.field = [0; 0; 0];

[~, ~, basis_sw, ~, ~, ~] = SW_proj(const, ion_sw, params_sw);

wav0 = basis_sw(:,:,1,1,1);
jx = real(wav0' * Jx * wav0);
jy = real(wav0' * Jy * wav0);
jz = real(wav0' * Jz * wav0);

g_eff_sw = zeros(1,3);
g_eff_sw(1) = 2 * gJ * abs(jx(1,2));
g_eff_sw(2) = 2 * gJ * abs(jy(1,2));
g_eff_sw(3) = 2 * gJ * abs(jz(1,1));

fprintf('  g_eff = [%.6f, %.6f, %.6f]\n', g_eff_sw);
fprintf('  Error = [%.2e, %.2e, %.2e]\n', (g_eff_sw - targetG)./targetG);

%% Summary
fprintf('\n======================================\n');
fprintf('Summary:\n');
fprintf('  projectDoublet: max|error| = %.2e\n', max(abs((g_eff_proj - targetG)./targetG)));
fprintf('  SW_proj:        max|error| = %.2e\n', max(abs((g_eff_sw - targetG)./targetG)));

if max(abs((g_eff_proj - targetG)./targetG)) < 0.1
    fprintf('\n✓ Baseline is GOOD - ready for optimization\n');
else
    fprintf('\n✗ Baseline is BAD - g-factors are far from target!\n');
    fprintf('  Check which parameters are actually active in MF_Er_CaWO4_v1b.m\n');
end
fprintf('======================================\n');
