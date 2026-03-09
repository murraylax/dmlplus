# dmlplus

**dmlplus** is a Dynamic Macroeconomic Library written in C++. It provides tools for solving and simulating linear dynamic general equilibrium (DSGE) models, numerical optimization, and matrix/vector utilities. An R interface package (`gensysR`) is also included.

> This is an early work in progress.

**Author:** James Murray — james@murraylax.org  
**License:** GNU General Public License v3

---

## Dependencies

| Library | Purpose |
|---|---|
| [Eigen3](https://eigen.tuxfamily.org/) | Matrix and vector types used throughout |
| [GSL](https://www.gnu.org/software/gsl/) | Multivariate root finding and minimization |
| [LAPACK / LAPACKE](https://www.netlib.org/lapack/) | QZ decomposition via `LAPACKE_zgges` |
| [BLAS / CBLAS](https://www.netlib.org/blas/) | Required by LAPACK |

On Debian/Ubuntu these can be installed with:

```bash
sudo apt install libeigen3-dev libgsl-dev liblapacke-dev liblapack-dev libblas-dev
```

On Arch Linux, these can be installed with:

```bash
sudo pacman -S eigen gsl lapack blas
```

---

## Building

The project uses a `Makefile` and requires a C++23-capable compiler (e.g., `g++`).

### Build the static library

```bash
make dml
```

This compiles all modules and produces `libdml.a`.

### Build individual modules

This library includes the following modules:

- **`utils`** — Matrix/vector utilities for converting between Eigen, GSL, and LAPACK data structures, plus file I/O and timing functions
- **`qz`** — QZ (generalized Schur) decomposition wrapper around LAPACK's `zgges` function
- **`gensys`** — Linear DSGE model solver implementing the Sims (2002) method
- **`dml_multiroot`** — Multivariate root finder (wrapper around GSL)
- **`dml_multimin`** — Multivariate minimizer with optional box constraints (wrapper around GSL)

```bash
make utils.o
make qz.o
make gensys.o
make dml_multiroot.o
make dml_multimin.o
```

### Build and run test programs

Test programs demonstrate the library's capabilities:

- **`testgensys.out`** — Solves a simple 3-equation New Keynesian DSGE model with an IS curve, Phillips curve, and Taylor rule; displays system matrices, solution matrices (G, M, D), and impulse responses
- **`test_multiroot.out`** — Demonstrates multivariate root finding
- **`test_multimin.out`** — Demonstrates multivariate minimization

```bash
make test_multiroot.out   # numeric/test_multiroot.out
make test_multimin.out    # numeric/test_multimin.out
make testgensys.out       # gensys/testgensys.out
```

### Clean build artifacts

```bash
make clean
```

---

## Project Structure

```
dmlplus/
├── dml.h                    # Umbrella header — include this to use the full library
├── Makefile
├── utils/
│   ├── utils.h
│   └── utils.cpp
├── gensys/
│   ├── qz.h / qz.cpp
│   ├── gensys.h / gensys.cpp
│   ├── nksys.h / nksys.cpp
│   ├── gensysR.cpp          # Rcpp interface
│   ├── gensysR.R            # Companion R script
│   ├── run_eigen.cpp        # Standalone example
│   └── gensysR/             # R package (Rcpp)
├── numeric/
│   ├── dml_multiroot.h / dml_multiroot.cpp
│   ├── dml_multimin.h / dml_multimin.cpp
│   ├── test_multiroot.cpp
│   └── test_multimin.cpp
```

---

## Modules

### `dml.h` — Umbrella Header

Include a single header to pull in the entire library:

```cpp
#include <dml.h>
```

---

### `utils/` — Matrix and Vector Utilities

**Files:** `utils/utils.h`, `utils/utils.cpp`

A collection of utilities for converting between Eigen, GSL, and LAPACK data structures, writing output files, and timing code.

**Key functions:**

| Function | Description |
|---|---|
| `write_eigen_csv(mat, filepath)` | Write an `Eigen::MatrixXd` or `VectorXd` to a CSV file |
| `write_eigen_csv(mat, varnames, filepath)` | Write a matrix to CSV with column variable names |
| `copy_vector_to_gsl` / `copy_gsl_to_vector` | Convert between `Eigen::VectorXd` and `gsl_vector` |
| `copy_matrix_to_gsl` / `copy_gsl_to_matrix` | Convert between `Eigen::MatrixXd` and `gsl_matrix` |
| `copy_matrix_to_lapack_complex` / `copy_lapack_complex_to_matrix` | Convert between Eigen complex matrices and LAPACK arrays |
| `print_vector_to_rcode` / `print_matrix_to_rcode` | Generate R assignment code as a string (real or complex) |
| `start_timer()` / `stop_timer(start)` | Simple wall-clock timing utilities |

---

### `gensys/` — Linear DSGE Model Solver

Implements the **Sims (2002)** method for solving linear rational expectations models.

#### `qz.h` / `qz.cpp` — QZ Decomposition

A C++ wrapper around the LAPACK function `LAPACKE_zgges` (generalized Schur decomposition) operating on Eigen matrices.

```cpp
int qz(Eigen::MatrixXcd& Q, Eigen::MatrixXcd& Z,
       Eigen::MatrixXcd& S, Eigen::MatrixXcd& T,
       Eigen::VectorXd& Lambda,
       const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);
```

Returns the generalized Schur form of the matrix pencil $(A, B)$, with eigenvalues ordered so that unstable roots (modulus $\geq 1$) appear last.

Also provides `check_colspace(A, B)` to test whether the column space of `B` lies within the column space of `A` using SVD.

---

#### `gensys.h` / `gensys.cpp` — General System Solver

Solves a linear DSGE model in the form:

$$\Gamma_0 x_t = C + \Gamma_1 x_{t-1} + \Psi z_t + \Pi \eta_t$$

where $x_t$ is the vector of state variables, $z_t$ is a vector of exogenous shocks, and $\eta_t = x_t - E_{t-1} x_t$ is the vector of endogenous expectation errors.

The solution takes the form:

$$x_t = D + G x_{t-1} + M z_t$$

**Main functions:**

| Function | Description |
|---|---|
| `gensys(G, M, D, Γ₀, Γ₁, Ψ, Π, C)` | Solve the system; returns 0 (unique), > 0 (indeterminate), or -1 (no solution) |
| `gensys_qzdetails(...)` | Like `gensys` but also returns the QZ decomposition matrices and eigenvalues |
| `checksys(Γ₀, Γ₁, Ψ, Π, C)` | Check existence/uniqueness without computing the full solution |
| `gensys_irf(G, M, shock, idx, T)` | Compute impulse response for a single shock over `T` periods |
| `gensys_irf(G, M, T)` | Compute impulse responses for all shocks (magnitude 0.01) over `T` periods |
| `write_irf(IRF, varnames, filepath)` | Write impulse responses to a text file for import into R |
| `write_irf_to_csvfile(...)` | Write impulse responses to an open CSV file stream |
| `check_colspace(A, B)` / `check_colspace(A, B, tol)` | Test whether col(B) ⊆ col(A) |

**Return codes for `gensys` and `checksys`:**

| Return value | Meaning |
|---|---|
| `0` | Unique solution exists |
| `> 0` | Indeterminate — that many loose endogenous variables |
| `-1` | No solution exists |

---

#### `nksys.h` / `nksys.cpp` — New Keynesian DSGE Model

A fully specified New Keynesian DSGE model featuring **two agent types**: optimizing households and rule-of-thumb households. State variable and expectation error indices are defined as enums for readability.

Variables include consumption, labor, bonds, investment, capital, wages, the nominal interest rate, the real return on capital, inflation, Tobin's Q, fiscal variables, and various expectational terms.

---

#### `gensysR/` — R Package

An R package that wraps the `gensys` solver using **Rcpp**, making `gensys`, `gensys_qzdetails`, and `checksys` callable directly from R. See `gensys/gensysR/` for the package source.

---

#### `testgensys.cpp` — Test Program

A standalone test program that demonstrates the `gensys` solver with a simple 3-equation New Keynesian model.

**Model structure:**
- **Variables:** Output gap ($y_t$), inflation ($\pi_t$), interest rate ($r_t$), plus AR(1) demand and cost-push shocks
- **Equations:**
  1. IS curve: $y_t = E_t[y_{t+1}] - \sigma(r_t - E_t[\pi_{t+1}]) + \varepsilon_{y,t}$
  2. Phillips curve: $\pi_t = \beta E_t[\pi_{t+1}] + \kappa y_t + \varepsilon_{\pi,t}$
  3. Taylor rule: $r_t = \rho_r r_{t-1} + (1 - \rho_r)(\phi_\pi \pi_t + \phi_y y_t) + \varepsilon_{r,t}$
- **Shocks:** Demand shock ($\varepsilon_{y,t}$), cost-push shock ($\varepsilon_{\pi,t}$), monetary shock ($\varepsilon_{r,t}$)

**Output:**
- System matrices ($\Gamma_0, \Gamma_1, \Psi, \Pi$)
- Solution matrices ($G, M, D$), where $x_t = D + G x_{t-1} + M z_t$
- Solution status (unique, indeterminate, or non-existent)
- Impulse responses for 20 periods showing the dynamic effects of each shock

Build and run:
```bash
make testgensys.out
./gensys/testgensys.out
```

---

### `numeric/` — Numerical Methods

Eigen-based wrappers around GSL's multivariate solvers.

#### `dml_multiroot.h` / `dml_multiroot.cpp` — Multivariate Root Finder

Finds the zero of a vector-valued function $f : \mathbb{R}^n \to \mathbb{R}^n$ using GSL's `gsl_multiroot` routines.

```cpp
Eigen::VectorXd dml_multiroot(
    const Eigen::VectorXd& initial_guess,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void*)> func,
    const void* params,
    bool verbose = false);

// Overload with user-supplied Jacobian
Eigen::VectorXd dml_multiroot(
    const Eigen::VectorXd& initial_guess,
    std::function<Eigen::VectorXd(...)> func,
    std::function<Eigen::MatrixXd(...)> jacfunc,
    const void* params, bool verbose);
```

The `params` pointer allows passing arbitrary data to the objective function.

---

#### `dml_multimin.h` / `dml_multimin.cpp` — Multivariate Minimizer

Minimizes a scalar function $f : \mathbb{R}^n \to \mathbb{R}$ using GSL's `gsl_multimin` routines. Supports optional box constraints (element-wise lower and upper bounds).

```cpp
Eigen::VectorXd dml_multimin(
    const Eigen::VectorXd& initial_guess,
    std::function<double(const Eigen::VectorXd&, const void*)> func,
    const void* params,
    bool verbose = false);

// Overload with box constraints
Eigen::VectorXd dml_multimin(
    const Eigen::VectorXd& initial_guess,
    const Eigen::VectorXd& lower_bounds,
    const Eigen::VectorXd& upper_bounds,
    std::function<double(const Eigen::VectorXd&, const void*)> func,
    const void* params,
    bool verbose = false);
```

---

## Quick Start Example

```cpp
#include <dml.h>
#include <Eigen/Dense>
#include <iostream>

int main() {
    int nvar = 3, nshocks = 1, nendo = 1;

    Eigen::MatrixXd Gamma0 = ...; // nvar x nvar
    Eigen::MatrixXd Gamma1 = ...; // nvar x nvar
    Eigen::MatrixXd Psi    = ...; // nvar x nshocks
    Eigen::MatrixXd Pi     = ...; // nvar x nendo
    Eigen::VectorXd C      = Eigen::VectorXd::Zero(nvar);

    Eigen::MatrixXd G, M;
    Eigen::VectorXd D;

    int retval = gensys(G, M, D, Gamma0, Gamma1, Psi, Pi, C);

    if (retval == 0) {
        std::cout << "Unique solution found.\n";
        Eigen::MatrixXd irf = gensys_irf(G, M, 0.01, 0, 40);
    } else if (retval > 0) {
        std::cout << "Indeterminate: " << retval << " loose endogenous variables.\n";
    } else {
        std::cout << "No solution exists.\n";
    }
}
```

---

## Contact

James Murray — james@murraylax.org

