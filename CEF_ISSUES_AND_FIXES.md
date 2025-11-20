# CEF Fitting Issues and Solutions

## Problem Summary
The eigenstates and excitation spectrum from CEF parameters don't match the effective spin-1/2 model, despite getting roughly correct effective g-factors.

## Identified Issues

### Issue 1: **Incorrect Hyperfine Coupling Constant**

**Location:** `MF_Er_CaWO4_v1b.m:137-138`

**Current code:**
```matlab
A0_MHz = -125.9; % scalar hyperfine constant
A = (A0_MHz/1000) * const.Gh2mV * [1 1 1];
```

**Problem:**
The effective model uses strongly anisotropic hyperfine:
- `A_eff = [-871.1, -871.1, -130.3] MHz` (line 150)

The relationship between isotropic A0 and effective A_eff is:
```
A_eff_i = A0 * g_eff_i / (2 * gJ)
```

Therefore:
```
A0 = A_eff * (2 * gJ) / g_eff
```

For target values:
- `A0_x = -871.1 * 2*1.2 / 8.3 ≈ -251.8 MHz`
- `A0_y = -871.1 * 2*1.2 / 8.3 ≈ -251.8 MHz`
- `A0_z = -130.3 * 2*1.2 / 1.26 ≈ -248.3 MHz`
- **Mean A0 ≈ -250.6 MHz** (not -125.9 MHz!)

**Solution:**
```matlab
A0_MHz = -250.6; % Corrected scalar hyperfine constant
A = (A0_MHz/1000) * const.Gh2mV * [1 1 1];
```

### Issue 2: **Inconsistent Use of SW_proj vs projectDoublet**

**Locations:**
- `CEF_fitting.m:383` uses `projectDoublet(Hcf, Jx, Jy, Jz)` (zero-field only)
- `MF_Er_CaWO4_v1b.m:134` uses `SW_proj(const, ion, params)` (includes 2nd-order perturbation)

**Problem:**
- `projectDoublet()` simply diagonalizes the CEF Hamiltonian and takes the lowest two states
- `SW_proj()` performs Schrieffer-Wolff 2nd-order perturbation theory, which includes virtual transitions to excited states

These methods can give slightly different results, especially:
1. At finite magnetic fields
2. When excited states are relatively close in energy
3. When calculating field-dependent properties

**Impact:**
The projections might differ by small amounts, leading to slightly different effective operators and energy level structures.

**Solution:**
Use the **same projection method** in both files. Recommended approach:

**Option A (Simpler):** Modify `CEF_fitting.m` to also use `SW_proj` for consistency
**Option B (Current):** Keep `projectDoublet` but ensure it's also used in `MF_Er_CaWO4_v1b.m`

For consistency with the main script, I recommend **Option A**.

### Issue 3: **CEF Fitting Doesn't Optimize A0**

**Location:** `CEF_fitting.m:44`

**Current code:**
```matlab
reportHyperfine = false; % if true, print inferred A0 and A_eff
```

**Problem:**
The CEF fitting:
1. Optimizes CEF parameters to match `targetG = [8.3 8.3 1.26]`
2. Uses the zero-field 16-level spectrum as a constraint
3. But **does NOT** explicitly fit the hyperfine coupling A0

The variable `targetAeff_MHz = [-871.1, -871.1, -130.3]` (line 43) is defined but **not used in the optimization** - it's only used in the optional reporting function `inferA0AndAeff()`.

**Impact:**
The fitted CEF parameters might give the correct g-tensor but wrong hyperfine structure because A0 was never optimized.

**Solution:**
Either:
1. Calculate the correct A0 from the target A_eff and fitted g_eff
2. Or, modify the fitting to include A0 as a fitted parameter

## Recommended Fixes

### Fix 1: Correct A0 in MF_Er_CaWO4_v1b.m

**File:** `MF_Er_CaWO4_v1b.m`
**Lines:** 137-138

**Replace:**
```matlab
A0_MHz = -125.9; % scalar hyperfine constant for 4I15/2 of 167Er3+
A = (A0_MHz/1000) * const.Gh2mV * [1 1 1];
```

**With:**
```matlab
% Calculate A0 from target effective hyperfine and fitted g-factors
% A_eff = A0 * g_eff / (2*gJ)  =>  A0 = A_eff * (2*gJ) / g_eff
% Target: A_eff = [-871.1, -871.1, -130.3] MHz, g_eff = [8.3, 8.3, 1.26]
% For Er3+: gJ ≈ 1.2
% A0 ≈ mean([-871.1*2*1.2/8.3, -871.1*2*1.2/8.3, -130.3*2*1.2/1.26]) ≈ -250.6 MHz
A0_MHz = -250.6; % Corrected scalar hyperfine constant
A = (A0_MHz/1000) * const.Gh2mV * [1 1 1];
```

### Fix 2: Add A0 calculation helper function

Create a new helper function to automatically calculate the correct A0:

**File:** `calculate_A0_from_Aeff.m` (new file)

```matlab
function A0_MHz = calculate_A0_from_Aeff(A_eff_MHz, g_eff, gJ)
%CALCULATE_A0_FROM_AEFF Convert effective anisotropic A to isotropic A0
%   A0_MHz = calculate_A0_from_Aeff(A_eff_MHz, g_eff, gJ)
%
%   The relationship is: A_eff_i = A0 * g_eff_i / (2*gJ)
%   Therefore: A0 = A_eff_i * (2*gJ) / g_eff_i
%
%   This function calculates A0 for each component and returns the mean.

    A0_components = A_eff_MHz .* (2*gJ) ./ g_eff;
    A0_MHz = mean(A0_components);

    % Check consistency
    if std(A0_components) / abs(mean(A0_components)) > 0.05
        warning('A0 values from different components vary by more than 5%%: [%.2f, %.2f, %.2f] MHz', ...
            A0_components);
    end
end
```

### Fix 3: Update CEF_fitting.m to report/use correct A0

**File:** `CEF_fitting.m`
**Line:** 44

**Change:**
```matlab
reportHyperfine = false; % if true, print inferred A0 and A_eff
```

**To:**
```matlab
reportHyperfine = true; % Always report inferred A0 and A_eff for verification
```

Also, add code after line 205 to automatically save the recommended A0:

```matlab
% Save recommended A0 for use in MF_Er_CaWO4_v1b.m
stats.A0_recommended_MHz = A0_rec_MHz;
fprintf('\n*** IMPORTANT: Use A0 = %.2f MHz in MF_Er_CaWO4_v1b.m ***\n', A0_rec_MHz);
```

### Fix 4: Ensure consistent projection method (Optional but recommended)

**Option A:** Modify `MF_Er_CaWO4_v1b.m` to use `projectDoublet` instead of `SW_proj`

**Lines:** 132-135

**Replace:**
```matlab
for ii = 1:size(Bfield,2)
    params.field = Bfield(:,ii);
    [~, ham_E(:,:,ii,1), basis(:,:,ii,1), ~, ~, ~] = SW_proj(const, ion, params);
end
```

**With:**
```matlab
% For consistency with CEF_fitting, use projectDoublet at each field point
for ii = 1:size(Bfield,2)
    % Build total Hamiltonian including Zeeman
    Hz = -const.J2meV * const.muB * ion.gLande * ...
         (Bfield(1,ii)*ion.Jx + Bfield(2,ii)*ion.Jy + Bfield(3,ii)*ion.Jz);
    Htot = ion.Hcf + Hz;

    % Project onto lowest doublet
    [proj_temp, ~] = projectDoublet(Htot, ion.Jx, ion.Jy, ion.Jz);

    % Extract basis and energies
    basis(:,:,ii,1) = proj_temp.states;
    ham_E(:,:,ii,1) = diag(proj_temp.energies - min(proj_temp.energies));
end
```

**Option B:** Modify `CEF_fitting.m` to use `SW_proj` (more complex)

This would require restructuring the CEF_fitting optimization to use SW_proj, which is more involved.

## Testing the Fixes

After applying the fixes:

1. Run the diagnostic script: `diagnose_CEF_mismatch.m`
2. Check that:
   - The 16-level hyperfine spectrum matches between CEF and effective models
   - The g-factors are consistent (already working)
   - The excitation spectrum at various fields matches

3. Run `MF_Er_CaWO4_v1b.m` with:
   - `Options.CEF = true` (CEF model)
   - `Options.CEF = false` (effective model)

   And verify the transition energies match.

## Expected Impact

### Before fixes:
- ✓ g-factors: match (~8.3, ~8.3, ~1.26)
- ✗ Hyperfine spectrum: factor of ~2 error
- ✗ Transition energies: significant discrepancies

### After fixes:
- ✓ g-factors: match
- ✓ Hyperfine spectrum: should match within ~1-5%
- ✓ Transition energies: should match within numerical precision

## Summary

The main issue is **A0 is off by a factor of ~2** in `MF_Er_CaWO4_v1b.m` (line 137).

**Quick fix:** Change `A0_MHz = -125.9;` to `A0_MHz = -250.6;`

This should immediately improve agreement between the CEF and effective models.
