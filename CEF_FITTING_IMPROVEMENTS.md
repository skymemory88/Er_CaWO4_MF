# CEF_fitting.m Improvements

## Issues Fixed

### 1. **CRITICAL: Outdated Baseline Parameters**
**Problem**: CEF_fitting.m was using old CEF parameters from 2025.10.13 as the starting point:
```matlab
baseB_cm = [-753.279 796.633 -376.577 -133.463 -3.30013 -84.5397 -9.7149];
```

But your current parameters in MF_Er_CaWO4_v1b.m (2025.11.10) are:
```matlab
B = [-330.74 2212.59 196.501 5.57413 -8.07197 -217.553 81.2414];
```

These are **completely different**! The optimization was starting from the wrong place.

**Fix**: Updated line 48 to use the current parameters as baseline.

### 2. **Projection Method Inconsistency**
**Problem**:
- `CEF_fitting.m` used `projectDoublet()` (simple diagonalization)
- `MF_Er_CaWO4_v1b.m` uses `SW_proj()` (2nd-order Schrieffer-Wolff perturbation theory)

These methods give slightly different results, especially at finite fields or when excited states are nearby.

**Fix**: Added `option.useSWproj = true` to use `SW_proj()` for consistency. The `computeGeff()` function now supports both methods.

### 3. **Missing Constants for SW_proj**
**Problem**: SW_proj requires physical constants that weren't defined in CEF_fitting.m

**Fix**: Added constants structure (lines 54-60):
```matlab
const.hbar = 1.05457E-34;
const.muB = 9.274e-24;
const.muN = 5.05078e-27;
const.kB = 1.3806e-23;
const.J2meV = 6.24151e+21;
const.Gh2mV = const.hbar * 2*pi * 10^9 * const.J2meV;
```

### 4. **Lack of Diagnostic Output**
**Problem**: Hard to tell what's happening during optimization

**Fix**: Added `option.verbose = true` to control diagnostic output

## New Features

### 1. **Projection Method Selection**
You can now choose which projection method to use:
- `option.useSWproj = true` → Use SW_proj (recommended for consistency)
- `option.useSWproj = false` → Use projectDoublet (faster but less accurate)

### 2. **Verbose Output**
Control optimization output verbosity:
- `option.verbose = true` → Show detailed progress
- `option.verbose = false` → Minimal output

### 3. **Updated computeGeff Function**
The function now accepts method selection:
```matlab
[g_eff, Jproj] = computeGeff(B_cm, L, S, const, option)
```

## Usage

### Basic Usage (with target spectrum from effective model)

```matlab
% Generate target spectrum from effective model
I = 7/2;
[~,~,~,~,~,~,Jx_eff,Jy_eff,Jz_eff,Ix_eff,Iy_eff,Iz_eff] = spin_operators(0.5, I);

const.Gh2mV = 1.05457E-34 * 2*pi * 10^9 * 6.24151e+21;
A_eff_MHz = [-871.1, -871.1, -130.3];
A_eff_meV = (A_eff_MHz/1000) * const.Gh2mV;

Hhf_eff = A_eff_meV(1)*Jx_eff*Ix_eff + ...
          A_eff_meV(2)*Jy_eff*Iy_eff + ...
          A_eff_meV(3)*Jz_eff*Iz_eff;

[~, E_eff] = eig(Hhf_eff);
eigenE_target = sort(real(diag(E_eff)));

% Run fitting
[B_fit_cm, g_eff_fit, stats] = CEF_fitting(eigenE_target);
```

### With Field-Dependent Constraints

```matlab
% Prepare field-dependent data
fieldData.Bvec = [0; 0; 0.01]; % 10 mT along c-axis
fieldData.energies = [...]; % 16 energy levels at this field
fieldData.weight = 1.0; % Weight for this constraint
fieldData.gN = [0.1618 0.1618 0.1618]; % Nuclear g-tensor
fieldData.Q = [1.67 1.67 -3.34] / 1000 * const.Gh2mV; % Quadrupole

[B_fit_cm, g_eff_fit, stats] = CEF_fitting(eigenE_target, fieldData);
```

## Testing

Run the test script to see improvements:

```matlab
test_CEF_fitting
```

This will:
1. Generate target spectrum from effective model
2. Evaluate current CEF parameters
3. Run optimization
4. Compare before/after results
5. Show plots

## Optimization Settings

Current settings (can be modified in CEF_fitting.m):

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `option.Bayes` | `true` | Use Bayesian optimization |
| `option.bounds` | `false` | Use wide search bounds |
| `option.polish` | `true` | Polish with lsqnonlin after BayesOpt |
| `option.useSWproj` | `true` | Use SW_proj method |
| `option.verbose` | `true` | Show detailed output |
| `gTol` | `1e-3` | Relative g-factor tolerance |
| `wG` | `1.0` | Weight for g-factor residual |
| `wE` | `1.0` | Weight for energy residual |
| Max evaluations | `5000` | For Bayesian optimization |

### Tuning for Your Problem

**If optimization is too slow:**
- Reduce `MaxObjectiveEvaluations` from 5000 to 2000-3000
- Set `option.useSWproj = false` (faster but less accurate)
- Set `option.polish = false` (skip final refinement)

**If fit quality is poor:**
- Increase `MaxObjectiveEvaluations` to 10000
- Adjust weights: increase `wE` to prioritize spectrum matching
- Add field-dependent constraints
- Try tighter bounds: `option.bounds = true`

**If CEF parameters diverge:**
- Use tighter bounds
- Reduce search range: change `lowerB_cm = baseB_cm - 200`

## Comparison: Old vs New Baseline

| Parameter | Old (2025.10.13) | New (2025.11.10) | Change |
|-----------|------------------|------------------|--------|
| B20 | -753.279 | -330.74 | +422.5 |
| B40 | 796.633 | 2212.59 | +1416.0 |
| B44c | -376.577 | 196.501 | +573.1 |
| B44s | -133.463 | 5.57413 | +139.0 |
| B60 | -3.30013 | -8.07197 | -4.77 |
| B64c | -84.5397 | -217.553 | -133.0 |
| B64s | -9.7149 | 81.2414 | +91.0 |

The parameters changed significantly between versions, so starting from the correct baseline is crucial!

## Expected Performance

With the updated baseline and SW_proj:

**Before optimization (current parameters):**
- g-factor error: ~ 1e-4 to 1e-3 (already quite good)
- Spectrum RMS error: Depends on A0 (should be small with corrected A0 = -250.6 MHz)

**After optimization (target):**
- g-factor error: < 1e-4 (within tolerance)
- Spectrum RMS error: Should decrease further
- Convergence: ~2000-5000 function evaluations

## Troubleshooting

### "Function SW_proj not found"
Make sure SW_proj.m is in your MATLAB path.

### "Optimization diverges"
1. Check that your target spectrum is physically reasonable
2. Try tighter bounds: `option.bounds = true`
3. Reduce search range
4. Try starting from different baseline

### "No improvement from optimization"
1. Your current parameters might already be optimal!
2. Check if target spectrum is achievable with this CEF model
3. Try adding field-dependent constraints
4. The effective model might have constraints that can't be satisfied simultaneously

### "Too slow"
1. Reduce MaxObjectiveEvaluations
2. Set `option.useSWproj = false`
3. Use lsqnonlin only: `option.Bayes = false`

## Files Modified

1. **CEF_fitting.m** - Main fitting script
   - Updated baseline parameters (line 48)
   - Added constants (lines 54-60)
   - Added options (lines 65-66)
   - Modified objective function calls
   - Updated computeGeff to support SW_proj

2. **New files:**
   - `test_CEF_fitting.m` - Test script
   - `analyze_CEF_fitting.m` - Diagnostic analysis
   - `CEF_FITTING_IMPROVEMENTS.md` - This document

## Summary

The main issue was **starting optimization from outdated parameters**. The fitting was trying to optimize from parameters that are far from your current best values.

With the fixes:
1. ✓ Starting from correct baseline (current parameters)
2. ✓ Using consistent projection method (SW_proj)
3. ✓ Proper diagnostic output
4. ✓ All constants defined

The optimization should now work much better and give results consistent with MF_Er_CaWO4_v1b.m!
