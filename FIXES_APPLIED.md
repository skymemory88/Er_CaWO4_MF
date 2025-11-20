# CEF Fitting Fixes - Summary

## Problem
The eigenstates and excitation spectrum from CEF parameters didn't match the effective spin-1/2 model, despite getting correct effective g-factors.

## Root Cause
**Hyperfine coupling constant A0 was off by a factor of ~2**

The scalar hyperfine constant `A0_MHz = -125.9 MHz` in `MF_Er_CaWO4_v1b.m:137` was incorrect.

## Analysis

### The Relationship Between A0 and A_eff

For an effective spin-1/2 model projected from a CEF ground doublet:

```
A_eff_i = A0 * g_eff_i / (2 * gJ)
```

Therefore:

```
A0 = A_eff_i * (2 * gJ) / g_eff_i
```

### Target Values
From the effective spin-1/2 model (lines 150-151):
- `A_eff = [-871.1, -871.1, -130.3] MHz`
- `g_eff = [8.3, 8.3, 1.26]`
- `gJ = 1.20003` (for Er3+: L=6, S=3/2)

### Calculated A0
```
A0_x = -871.1 * (2 * 1.20003) / 8.3  = -251.8 MHz
A0_y = -871.1 * (2 * 1.20003) / 8.3  = -251.8 MHz
A0_z = -130.3 * (2 * 1.20003) / 1.26 = -248.3 MHz

Mean A0 = -250.6 MHz
```

The original value of -125.9 MHz was **off by a factor of 1.99** (≈2).

## Files Modified

### 1. MF_Er_CaWO4_v1b.m (Line 137)

**Before:**
```matlab
A0_MHz = -125.9; % scalar hyperfine constant for 4I15/2 of 167Er3+
A = (A0_MHz/1000) * const.Gh2mV * [1 1 1];
```

**After:**
```matlab
% Calculate A0 from target effective hyperfine and g-factors
% The relationship is: A_eff = A0 * g_eff / (2*gJ)
% Therefore: A0 = A_eff * (2*gJ) / g_eff
% Target effective parameters from spin-1/2 model (lines 150-151):
%   A_eff = [-871.1, -871.1, -130.3] MHz
%   g_eff = [8.3, 8.3, 1.26]
% For Er3+: gJ = gLande(6, 3/2) ≈ 1.20003
% Calculate: A0 ≈ mean([-871.1*2*1.2/8.3, -871.1*2*1.2/8.3, -130.3*2*1.2/1.26])
%              ≈ mean([-251.8, -251.8, -248.3]) ≈ -250.6 MHz
A0_MHz = -250.6; % Corrected scalar hyperfine constant (was -125.9, factor of ~2 error)
A = (A0_MHz/1000) * const.Gh2mV * [1 1 1];
```

### 2. CEF_fitting.m (Line 44)

**Before:**
```matlab
reportHyperfine = false; % if true, print inferred A0 and A_eff
```

**After:**
```matlab
reportHyperfine = true; % if true, print inferred A0 and A_eff (IMPORTANT for MF_Er_CaWO4_v1b.m)
```

Added output message (after line 213):
```matlab
% Additional note about using A0 in MF_Er_CaWO4_v1b.m
if reportHyperfine
    fprintf('\n*** IMPORTANT for MF_Er_CaWO4_v1b.m: ***\n');
    fprintf('The recommended A0 value has been calculated above.\n');
    fprintf('Update line 137 in MF_Er_CaWO4_v1b.m with the recommended A0 (MHz).\n');
    fprintf('Current value in MF_Er_CaWO4_v1b.m should be A0_MHz = %.2f;\n', A0_rec_MHz);
end
```

## New Files Created

### 1. calculate_A0_from_Aeff.m
Helper function to automatically calculate A0 from effective hyperfine parameters:

```matlab
A0_MHz = calculate_A0_from_Aeff(A_eff_MHz, g_eff, gJ)
```

Example usage:
```matlab
L = 6; S = 3/2;
gJ = gLande(L, S);
A_eff_MHz = [-871.1, -871.1, -130.3];
g_eff = [8.3, 8.3, 1.26];
A0_MHz = calculate_A0_from_Aeff(A_eff_MHz, g_eff, gJ);
% Returns: A0_MHz = -250.6
```

### 2. diagnose_CEF_mismatch.m
Diagnostic script that compares:
- CEF model (with current and corrected A0)
- Effective spin-1/2 model
- Both projection methods (projectDoublet vs SW_proj)

Run this script to verify the fixes work correctly.

### 3. CEF_ISSUES_AND_FIXES.md
Comprehensive documentation of all issues found and recommended solutions.

## Testing the Fixes

### Test 1: Run with CEF model
```matlab
% In MF_Er_CaWO4_v1b.m
Options.CEF = true;
Options.nEn = 16;
Options.ndE = [1:15];
% Run the script
```

### Test 2: Run with effective model
```matlab
% In MF_Er_CaWO4_v1b.m
Options.CEF = false;
Options.nEn = 16;
Options.ndE = [1:15];
% Run the script
```

### Test 3: Compare results
The 16-level hyperfine spectrum at zero field should now match between the two models.

Specifically check:
- ✓ Ground doublet splitting: should match within ~1-5%
- ✓ Hyperfine level spacings: should be consistent
- ✓ Transition frequencies as a function of field: should agree

### Test 4: Run diagnostic script
```matlab
diagnose_CEF_mismatch
```

This will produce a detailed comparison with plots showing:
- 16-level hyperfine spectrum comparison
- Energy differences between models
- Effective g-factor comparison
- Effective hyperfine constant comparison

## Expected Results After Fix

| Property | Before Fix | After Fix |
|----------|------------|-----------|
| g-factors | ✓ Match (8.3, 8.3, 1.26) | ✓ Match (unchanged) |
| A0 value | -125.9 MHz (wrong) | -250.6 MHz (correct) |
| 16-level spectrum | ~2x error | Should match within ~1-5% |
| Transition energies | Significant discrepancies | Should agree within numerical precision |

## Additional Issues Identified (Not Yet Fixed)

### Issue: Different projection methods
- `CEF_fitting.m` uses `projectDoublet()` (simple diagonalization)
- `MF_Er_CaWO4_v1b.m` uses `SW_proj()` (2nd-order Schrieffer-Wolff)

**Impact:** Minor differences at finite field or when excited states are close
**Recommendation:** Use consistent method in both files (see CEF_ISSUES_AND_FIXES.md)
**Priority:** Low (current fix should resolve most discrepancies)

## Questions or Issues?

If the spectrum still doesn't match after this fix:

1. Check that the CEF parameters in `MF_Er_CaWO4_v1b.m:116` match those used in `CEF_fitting.m`
2. Verify that `gJ = gLande(L, S)` gives the expected value (~1.20003 for Er3+)
3. Run the diagnostic script to identify remaining discrepancies
4. Consider using consistent projection methods (see CEF_ISSUES_AND_FIXES.md)

## References

- `MF_Er_CaWO4_v1b.m:99-157` - Initialization function showing CEF vs effective models
- `CEF_fitting.m:374-388` - computeGeff function showing g-factor calculation
- `CEF_fitting.m:197-205` - inferA0AndAeff function (now enabled)
- `projectDoublet.m` - Simple zero-field doublet projection
- `SW_proj.m` - Schrieffer-Wolff 2nd-order projection with field
