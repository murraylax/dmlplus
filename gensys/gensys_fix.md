# gensysR Solver Bug: Diagnosis and Fix

## Problem

The `gensysR` package (installed from `murraylax/dmlplus/gensys/gensysR`) produces incorrect solution matrices `G` and `M`. The solver reports a unique solution (`nloose = 0`), but the solution **does not satisfy the model equations**.

## Evidence

### 1. Predetermined variables jump on impact

Capital stocks `K` and `K_G` are predetermined (their accumulation equations place them in `Gamma0` with only lagged terms in `Gamma1`). On impact (`x_{-1} = 0`), the accumulation equation `Gamma0["Kevo", ] %*% Msol[, "u_GI"]` should equal zero, since only `K` has a nonzero entry in `Gamma0["Kevo", ]`. Instead:

```
Gamma0["Kevo", ] %*% Msol[, "u_GI"] = 0.01759  (should be 0)
Gamma0["KGevo", ] %*% Msol[, "u_GI"] = 0.00273  (should be 0)
```

### 2. Nearly every equation is violated on impact

The full system residual `Gamma0 %*% Msol[, "u_GI"] - Upsilon[, "u_GI"]` is nonzero for almost all 70 equations, with residuals up to O(10^-2).

### 3. Residuals are NOT in the column space of Xi (Pi)

In a correct gensys solution, `Gamma0 %*% M - Upsilon` must lie in the column space of `Xi` (the expectation error matrix), because the only slack in the system comes from expectation errors. The orthogonal component should be machine-zero. Instead:

```
max(abs(orth_residual_M)) = 0.467  (should be ~1e-14)
```

### 4. G and M residuals are enormous

```
max(abs(Gamma0 %*% Gsol - Gamma1)) = 14.4
max(abs(Gamma0 %*% Msol - Upsilon)) = 4.5
```

## Root Cause

The bug is in the C++ source file `gensys.cpp` in the `gensysR` package. The solution formulas for `G` and `M` are incomplete — they omit the off-diagonal blocks from the QZ decomposition.

### Current (incorrect) code

```cpp
// G = Z_1 S_{11}^{-1} T_{11} Z_1'
mcdGsol = mcdZ1 * mcdS11.solve(mcdT11 * mcdZ1.adjoint());

// M = Z_1 S_{11}^{-1} (Psi1 - Pi1 * Pi2^{-1} * Psi2)
mcdMsol = mcdZ1 * mcdS11.solve(mcdPsi1 - mcdPi1 * mcdPi2.solve(mcdPsi2));
```

### What's missing

The QZ decomposition gives upper-triangular `S` and `T`, so the off-diagonal blocks `S_12` and `T_12` are generally nonzero. The formulas above ignore:

1. **In G**: The `T_12` block contribution: `+ T_12 * Z_2'`
2. **In M**: The `S_12` block coupling the stable and unstable blocks, and the `Z_2` contribution to `x_t` from the unstable block solution.

### Corrected formulas (unique solution case, nloose = 0)

```cpp
// Forward solution for unstable block: w_{2,t} = phi_z * z_t
Eigen::MatrixXcd phi_z = -mcdPi2.solve(mcdPsi2);

// G: transition matrix
// w_{1,t} = S_11^{-1} (T_11 w_{1,t-1} + T_12 w_{2,t-1})
// w_{2,t} = 0 (no lagged dependence in unique solution when z_t=0)
// x_t = Z_1 w_{1,t} + Z_2 w_{2,t}
mcdGsol = mcdZ1 * mcdS11.solve(mcdT11 * mcdZ1.adjoint() + mcdT12 * mcdZ2.adjoint());

// M: impact matrix
// w_{1,t} = S_11^{-1} (Psi1 - S_12 * phi_z - Pi1 * Pi2^{-1} * Psi2) z_t
// w_{2,t} = phi_z * z_t
// x_t = Z_1 w_{1,t} + Z_2 w_{2,t}
mcdMsol = mcdZ1 * mcdS11.solve(
    mcdPsi1 - mcdS12 * phi_z - mcdPi1 * mcdPi2.solve(mcdPsi2))
    + mcdZ2 * phi_z;
```

Note: `S_12` is already extracted in the existing code but never used — a telltale sign of the omission.

## How to Verify the Fix

After fixing gensysR, re-solve the model and check:

```r
# All of these should be near machine epsilon (~1e-14)
max(abs(nk$Gamma0 %*% nksol$Gsol - nk$Gamma1))  # should be ~0 or in col(Xi)
max(abs(nk$Gamma0 %*% nksol$Msol - nk$Upsilon))  # should be ~0 or in col(Xi)

# Predetermined variables should not jump on impact
nksol$Msol["K", "u_GI"]    # should be 0
nksol$Msol["K_G", "u_GI"]  # should be 0

# Orthogonal residual check
Xi <- nk$Xi
proj <- Xi %*% MASS::ginv(Xi)
orth_resid <- (diag(nrow(proj)) - proj) %*% (nk$Gamma0 %*% nksol$Msol - nk$Upsilon)
max(abs(orth_resid))  # should be ~1e-14
```

## Context

- The model is a New Keynesian DSGE with heterogeneous agents, government investment, and multiple fiscal instruments (70 variables, 70 equations, 21 shocks, 11 expectation errors).
- The gensysR package is a C++ (Eigen) implementation of the Sims (2002) gensys algorithm, called from R via Rcpp.
- Reference: Chris Sims' original MATLAB `gensys.m` at <https://sims.princeton.edu/yftp/gensys/>.