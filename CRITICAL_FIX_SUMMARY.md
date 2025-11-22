# Critical Bug Fix Summary

## The Problem

The previous test run failed completely because CEF_fitting.m was using **the wrong baseline parameters**.

### What Happened

I mistakenly used **commented-out** parameters from line 116 of MF_Er_CaWO4_v1b.m:
```matlab
% Line 116 (COMMENTED, NOT ACTIVE):
B = [-330.74 2212.59 196.501 5.57413 -8.07197 -217.553 81.2414]
```

But the **actual active** parameters are on line 115:
```matlab
% Line 115 (ACTIVE):
B = [-979.402 1898.05 346.52 -272.893 -7.86205 -187.516 58.9021]
```

### Impact

The wrong baseline gave terrible g-factors:
```
Wrong baseline: g_eff = [6.936, 3.116, 1.264]
Target:         g_eff = [8.300, 8.300, 1.260]
Errors:                 [-16.4%, -62.5%, 0.3%]  ← VERY BAD!
```

The optimization:
- Started from a terrible place
- Got stuck in a local minimum
- Made things WORSE instead of better
- Took 7.6 hours due to slow SW_proj calls

## The Fixes

### 1. **Corrected Baseline Parameters** ✓
```matlab
% CEF_fitting.m line 49
baseB_cm = [-979.402 1898.05 346.52 -272.893 -7.86205 -187.516 58.9021];
```

### 2. **Switched to projectDoublet for Optimization** ✓
```matlab
option.useSWproj = false;  % Use fast projectDoublet during optimization
```
- SW_proj is ~100x slower
- projectDoublet is sufficient for optimization
- Can verify results with SW_proj afterward

### 3. **Improved Optimization Settings** ✓
```matlab
% Tighter tolerances (avoid getting stuck)
StepTolerance: 1e-3 → 1e-6
FunctionTolerance: 1e-3 → 1e-6

% Prioritize g-factor matching
wG: 1.0 → 100.0
```

### 4. **Added Verification Script** ✓
Created `verify_baseline.m` to quickly check if baseline parameters are good.

## What to Do Next

### Step 1: Verify the Baseline
Run this first to make sure the corrected baseline is good:
```matlab
verify_baseline
```

**Expected output:**
```
Target g_eff = [8.300, 8.300, 1.260]

Baseline CEF parameters (line 115):
  ...

Testing with projectDoublet:
  g_eff = [~8.3, ~8.3, ~1.26]
  Error = [<0.1, <0.1, <0.1]

✓ Baseline is GOOD - ready for optimization
```

If you see "✓ Baseline is GOOD", proceed to Step 2.

If you see "✗ Baseline is BAD", **STOP** and check which parameters are actually active in MF_Er_CaWO4_v1b.m line 115.

### Step 2: Run the Optimization
Once baseline is verified:
```matlab
test_CEF_fitting
```

**Expected behavior:**
- Should be MUCH faster (minutes instead of hours)
- Should start from reasonable g-factors
- Should improve the fit (not make it worse)
- Final g-factor errors should be < 1e-3

### Step 3: Verify with SW_proj
After optimization completes, you can verify the fitted parameters using SW_proj:
```matlab
% In MF_Er_CaWO4_v1b.m, set the fitted parameters from optimization
B = [fitted_values_from_CEF_fitting];  % line 115

% Run with CEF model
Options.CEF = true;
MF_Er_CaWO4_v1b

% Compare with effective model
Options.CEF = false;
MF_Er_CaWO4_v1b
```

The spectra should match much better than before.

## Why Did This Happen?

1. MF_Er_CaWO4_v1b.m has many commented CEF parameter sets from various sources
2. I grabbed the wrong commented line (116) instead of the active line (115)
3. The wrong parameters gave terrible g-factors but I didn't verify the baseline first
4. SW_proj being slow masked the problem (took so long the test seemed "normal")

## Prevention

**Always verify baseline parameters first!**

Use the new verification script:
```matlab
verify_baseline  % Quick check before running optimization
```

This will catch baseline parameter errors immediately instead of after hours of optimization.

## Summary of Changes

| Issue | Before | After |
|-------|--------|-------|
| Baseline source | Line 116 (commented) | Line 115 (active) |
| Baseline g-errors | [-16%, -62%, 0.3%] | Should be < 10% |
| Projection method | SW_proj (slow) | projectDoublet (fast) |
| Optimization time | 7.6 hours | ~minutes |
| Step tolerance | 1e-3 | 1e-6 |
| g-factor weight | 1.0 | 100.0 |
| Verification | None | verify_baseline.m |

## Files Modified

1. **CEF_fitting.m** - Corrected baseline + optimization improvements
2. **verify_baseline.m** - NEW: Quick baseline verification
3. **CRITICAL_FIX_SUMMARY.md** - This document

## Next Run Should...

✓ Start from good baseline parameters
✓ Be much faster (~minutes vs ~hours)
✓ Actually improve the fit
✓ Give g-factor errors < 1e-3
✓ Provide reasonable CEF parameters

If the next run still fails, check:
- Is line 115 really active in MF_Er_CaWO4_v1b.m?
- Does verify_baseline show "✓ Baseline is GOOD"?
- Are there any MATLAB errors or warnings?

---

**TL;DR**: Used wrong parameters from commented line instead of active line. Fixed. Run `verify_baseline` first, then `test_CEF_fitting`.
