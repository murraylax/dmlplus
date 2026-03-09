/*
 * Test file for gensys() function
 * 
 * Tests the solution of a simple three-equation New Keynesian model:
 * - Output gap (y)
 * - Inflation (π)
 * - Interest rate (i)
 * 
 * Model equations:
 * 1. IS curve: y_t = E_t[y_{t+1}] - σ(r_t - E_t[π_{t+1}]) + ε_y,t
 * 2. Phillips curve: π_t = β*E_t[π_{t+1}] + κ*y_t + ε_π,t
 * 3. Taylor rule: r_t = ρ_r*r_{t-1} + (1 - ρ_r) * (φ_π*π_t + φ_y*y_t) + ε_r,t
 */

#include "gensys.h"
#include <iostream>
#include <iomanip>

// Enums for variable and shock indices
enum Variable { i_y, i_pi, i_r, i_a, i_u, i_Ey, i_Epi, nvar };  // Output gap, Inflation, Interest rate, demand shock, inflation shock, Expected future output gap, expected future inflation
enum Equation { i_IS, i_PC, i_Taylor, i_A, i_U, i_EY, i_EPI, neq };  // IS curve, Phillips curve, Taylor rule, demand shock, inflation shock, expectational error for y, expectational error for pi
enum Shock { i_shock_y, i_shock_pi, i_shock_r, nshocks };  // Demand, Supply, Monetary shocks
enum ForwardVar { i_eta_y, i_eta_pi, nendo };  // Forward-looking variables

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Testing gensys() with a simple 3-equation New Keynesian model\n";
    std::cout << "============================================================\n\n";

    // Model parameters
    double sigma = 1.0;     // Intertemporal elasticity of substitution
    double beta = 0.99;     // Discount factor
    double kappa = 0.1;     // Slope of Phillips curve
    double phi_pi = 1.5;    // Taylor rule response to inflation
    double phi_y = 0.1;     // Taylor rule response to output gap
    double rho_r = 0.7;     // Interest rate smoothing
    double rho_a = 0.7;     // AR(1) coefficient for demand shock
    double rho_u = 0.7;     // AR(1) coefficient for inflation shock

    std::cout << "Model parameters:\n";
    std::cout << "  σ (sigma) = " << sigma << " (intertemporal elasticity)\n";
    std::cout << "  β (beta)  = " << beta << " (discount factor)\n";
    std::cout << "  κ (kappa) = " << kappa << " (Phillips curve slope)\n";
    std::cout << "  φ_π       = " << phi_pi << " (Taylor rule inflation response)\n";
    std::cout << "  φ_y       = " << phi_y << " (Taylor rule output response)\n\n";
    std::cout << "  ρ_r       = " << rho_r << " (interest rate smoothing)\n\n";
    std::cout << "  ρ_a       = " << rho_a << " (demand shock AR(1) coefficient)\n\n";
    std::cout << "  ρ_u       = " << rho_u << " (inflation shock AR(1) coefficient)\n\n";

    // For Sims notation with timing t and t-1, we rewrite:
    // Actually, let's be more careful. The standard form is backward-looking in time.
    // Let me reformulate more carefully.

    // Standard Sims timing: Γ_0 x_t = Γ_1 x_{t-1} + Ψ z_t + Π η_t
    // But our model has E_t[x_{t+1}] terms, so we need expectational errors
    // η_t = x_t - E_{t-1}[x_t]
    
    // After some algebra, the canonical form becomes:
    // Γ_0 * x_t = Γ_1 * x_{t-1} + Ψ * z_t + Π * η_t
    
    Eigen::MatrixXd Gamma0(nvar, nvar);
    Eigen::MatrixXd Gamma1(nvar, nvar);
    Eigen::MatrixXd Psi(nvar, nshocks);
    Eigen::MatrixXd Pi(nvar, nendo);
    Eigen::VectorXd C(nvar);

    // Initialize matrices to zero
    Gamma0.setZero();
    Gamma1.setZero();
    Psi.setZero();
    Pi.setZero();
    C.setZero();

    // Equation 1: IS curve
    // y_t = E_t[y_{t+1}] - σ*(r_t - E_t[π_{t+1}]) + ε_y,t
    Gamma0(i_IS, i_y) = 1.0;                // y_t coefficient
    Gamma0(i_IS, i_r) = sigma;              // r_t coefficient
    Gamma0(i_IS, i_Ey) = -1.0;              // E_t[y_{t+1}] coefficient
    Gamma0(i_IS, i_Epi) = -1.0 * sigma;     // E_t[π_{t+1}] coefficient
    Gamma0(i_IS, i_a) = -1.0;               // demand shock

    // Equation 2: Phillips curve
    // π_t = β*E_t[π_{t+1}] + κ*y_t + ε_π,t
    Gamma0(i_PC, i_pi) = 1.0;               // π_t coefficient
    Gamma0(i_PC, i_y) = -1.0 * kappa;       // y_t coefficient
    Gamma0(i_PC, i_Epi) = -1.0 * beta;      // E_t[π_{t+1}] coefficient
    Gamma0(i_PC, i_u) = -1.0;               // inflation shock

    // Equation 3: Taylor rule 
    // r_t = ρ_r*r_{t-1} + (1 - ρ_r) * (φ_π*π_t + φ_y*y_t) + ε_r,t
    Gamma0(i_Taylor, i_r) = 1.0;            // r_t coefficient
    Gamma0(i_Taylor, i_y) = -1.0 * (1.0 - rho_r) * phi_y;   // y_t coefficient
    Gamma0(i_Taylor, i_pi) = -1.0 * (1.0 - rho_r) * phi_pi; // π_t coefficient
    Gamma1(i_Taylor, i_r) = rho_r;          // r_{t-1} coefficient
    Psi(i_Taylor, i_shock_r) = 1.0;         // ε_r,t shock

    // Equation 4: Expectational error for y_t
    // y_t = E_{t-1}[y_t] + η_{y,t}
    Gamma0(i_EY, i_y) = 1.0;           
    Gamma1(i_EY, i_Ey) = 1.0;           
    Pi(i_EY, i_eta_y) = 1.0;          

    // Equation 5: Expectational error for π_t
    // π_t = E_{t-1}[π_t] + η_{π,t}
    Gamma0(i_EPI, i_pi) = 1.0;
    Gamma1(i_EPI, i_Epi) = 1.0;
    Pi(i_EPI, i_eta_pi) = 1.0;

    // Equation 6: AR(1) demand shock
    // a_t = ρ_a * a_{t-1} + ε_a,t
    Gamma0(i_A, i_a) = 1.0;
    Gamma1(i_A, i_a) = rho_a;
    Psi(i_A, i_shock_y) = 1.0;

    // Equation 7: AR(1) cost-push shock
    // u_t = ρ_u * u_{t-1} + ε_u,t
    Gamma0(i_U, i_u) = 1.0;
    Gamma1(i_U, i_u) = rho_u;
    Psi(i_U, i_shock_pi) = 1.0;

    std::cout << "System matrices:\n\n";
    std::cout << "Γ_0 =\n" << Gamma0 << "\n\n";
    std::cout << "Γ_1 =\n" << Gamma1 << "\n\n";
    std::cout << "Ψ =\n" << Psi << "\n\n";
    std::cout << "Π =\n" << Pi << "\n\n";
    std::cout << "C =\n" << C << "\n\n";

    // Solution matrices
    Eigen::MatrixXd G(nvar, nvar);
    Eigen::MatrixXd M(nvar, nshocks);
    Eigen::VectorXd D(nvar);

    std::cout << "Calling gensys()...\n";
    int result = gensys(G, M, D, Gamma0, Gamma1, Psi, Pi, C);

    std::cout << "\nResult: " << result << "\n";
    if (result == 0) {
        std::cout << "  => Unique solution found!\n\n";

        std::cout << "Solution matrices:\n\n";
        std::cout << "G (state transition matrix) =\n" << G << "\n\n";
        std::cout << "M (shock impact matrix) =\n" << M << "\n\n";
        std::cout << "D (constant vector) =\n" << D << "\n\n";

        std::cout << "Solution: x_t = D + G * x_{t-1} + M * z_t\n";
        std::cout << "where x_t = [y_t, π_t, r_t, E_t y_{t+1}, E_t π_{t+1}]' and z_t = [ε_{y,t}, ε_{π,t}, ε_{r,t}]'\n\n";

        // Compute impulse responses
        size_t nirf = 20;
        std::cout << "Computing impulse responses for " << nirf << " periods...\n";
        Eigen::MatrixXd irf = gensys_irf(G, M, nirf);
        
        std::cout << "Impulse responses computed. Shape: " << irf.rows() << " x " << irf.cols() << "\n";
        std::cout << "(Rows = periods * shocks, Cols = variables + shock_id + time)\n\n";

        // Print first few periods for shock 0 (demand shock)
        std::cout << "Impulse response to demand shock (ε_y = 0.01):\n";
        std::cout << "Shock      Period      y_t        π_t        r_t\n";
        std::cout << "---------------------------------------------------\n";
        for(int t = 0; t < 10 && t < nirf; t++) {
            std::cout << std::setw(4) << irf(t, nvar) << "  ";
            std::cout << std::setw(4) << irf(t, nvar+1) << "  ";
            std::cout << std::setw(10) << irf(t, i_y) << " ";
            std::cout << std::setw(10) << irf(t, i_pi) << " ";
            std::cout << std::setw(10) << irf(t, i_r) << "\n";
        }

    } else if (result > 0) {
        std::cout << "  => Indeterminacy! " << result << " loose endogenous variable(s).\n";
    } else {
        std::cout << "  => No solution exists.\n";
    }

    return 0;
}
