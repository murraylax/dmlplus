#' checksysR
#'
#' Use Sims (2002) method for checking for a unique solution a linear dynamic general equilibrium model
#'
#' Model has the following form: \Gamma_0 x_t = C + \Gamma_1 x_t-1 + \Psi z_t + \Pi \eta_t
#'
#' where...
#'  x_t is a vector of variables
#'  z_t is a vector of exogenous shocks
#'  \eta_t is a vector of endogenous expectation errors: \eta_t = x_t - E_t-1 x_t
#'
#'  C is a vector of constants
#'  \Gamma_0 and \Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
#'  \Psi is a matrix, dimension nvar x nshocks
#'  \Pi is a matrix, dimension nvar x nendo
#'
#' The solution is given as,
#'
#' x_t = D_sol + G_sol x_t-1 + M_sol z_t
#'
#' The method includes a return value for existence and uniqueness of solution:
#'   Return value > 0: Indeterminacy, there are `return_value` number of loose endogenous variables, ignore outputs for D_{sol}, G_{sol}, and M_{sol}
#'   Return value = 0: Unique solution is found. There are 0 loose endogenous variables.
#'   Return value = -1: No solution exists.
#'
#' @param Gamma0 (input) The coefficient matrix \Gamma_0
#' @param Gamma1 (input) The coefficient matrix \Gamma_1
#' @param Psi (input) The coefficient matrix \Psi
#' @param Pi (input) The coefficient matrix \Pi
#' @param C (input) The constant vector C
#'
#' @return Integer equal to the number of loose parameters. =0 when there is a unique solution, =-1 for no solution, >0 for indeterminacy
#' @export
#' @name checksysR
NULL




#' gensysR
#'
#' Use Sims (2002) method for solving a linear dynamic general equilibrium model
#'
#' Model has the following form: \Gamma_0 x_t = C + \Gamma_1 x_t-1 + \Psi z_t + \Pi \eta_t
#'
#' where...
#'  x_t is a vector of variables
#'  z_t is a vector of exogenous shocks
#'  \eta_t is a vector of endogenous expectation errors: \eta_t = x_t - E_t-1 x_t
#'
#'  C is a vector of constants
#'  \Gamma_0 and \Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
#'  \Psi is a matrix, dimension nvar x nshocks
#'  \Pi is a matrix, dimension nvar x nendo
#'
#' The solution is given as,
#'
#' x_t = D_sol + G_sol x_t-1 + M_sol z_t
#'
#' The method includes a return value for existence and uniqueness of solution:
#'   Return value > 0: Indeterminacy, there are `return_value` number of loose endogenous variables, ignore outputs for D_{sol}, G_{sol}, and M_{sol}
#'   Return value = 0: Unique solution is found. There are 0 loose endogenous variables.
#'   Return value = -1: No solution exists.
#'
#' @param Gamma0 The coefficient matrix \Gamma_0
#' @param Gamma1 The coefficient matrix \Gamma_1
#' @param Psi The coefficient matrix \Psi
#' @param Pi The coefficient matrix \Pi
#' @param C The constant vector C
#'
#' @return A list with the solution and an integer specifying whether the solution exists and is unique.
#'   The elements of the list are as follows:
#'      Gsol The solution matrix G_{sol}
#'      Msol The solution matrix M_{sol}
#'      Dsol The solution vector D_{sol}
#'      nloose The number of loose parameters. nloose=0 when there is a unique solution, nloose=-1 for no solution, nloose>0 for indeterminacy
#' @useDynLib gensysR
#' @export
#' @name gensysR
NULL


#' gensysR_qzdetails
#'
#' Use Sims (2002) method for solving a linear dynamic general equilibrium model
#'
#' Model has the following form: \Gamma_0 x_t = C + \Gamma_1 x_t-1 + \Psi z_t + \Pi \eta_t
#'
#' where...
#'  x_t is a vector of variables
#'  z_t is a vector of exogenous shocks
#'  \eta_t is a vector of endogenous expectation errors: \eta_t = x_t - E_t-1 x_t
#'
#'  C is a vector of constants
#'  \Gamma_0 and \Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
#'  \Psi is a matrix, dimension nvar x nshocks
#'  \Pi is a matrix, dimension nvar x nendo
#'
#' The solution is given as,
#'
#' x_t = D_sol + G_sol x_t-1 + M_sol z_t
#'
#' The method includes a return value for existence and uniqueness of solution:
#'   Return value > 0: Indeterminacy, there are `return_value` number of loose endogenous variables, ignore outputs for D_{sol}, G_{sol}, and M_{sol}
#'   Return value = 0: Unique solution is found. There are 0 loose endogenous variables.
#'   Return value = -1: No solution exists.
#'
#' Solution involves the QZ decomposition of Gamma0, Gamma1, where
#'   Q' Gamma_0 Z = S
#'   Q' Gamma_1 Z = T
#'   Q'Q = Q Q' = Z'Z = Z Z' = I
#'   Eigenvalues, lambda_i = abs(T_ii / S_ii)
#'   And eigenvalues are ordered where all those at the bottom are >= 1.0 in magnitude
#' 
#' @param Gamma0 The coefficient matrix \Gamma_0
#' @param Gamma1 The coefficient matrix \Gamma_1
#' @param Psi The coefficient matrix \Psi
#' @param Pi The coefficient matrix \Pi
#' @param C The constant vector C
#'
#' @return A list with the solution and an integer specifying whether the solution exists and is unique.
#'   The elements of the list are as follows:
#'      Gsol The solution matrix G_{sol}
#'      Msol The solution matrix M_{sol}
#'      Dsol The solution vector D_{sol}
#'      S The QZ decomposition matrix S above
#'      T The QZ decomposition matrix T above
#'      Q The QZ decomposition matrix Q above
#'      Z The QZ decomposition matrix Z above
#'      Lambda The real magnitude of the eigenvalues
#'      nunstable The number of unstable eigenvalues
#'      nloose The number of loose parameters. nloose=0 when there is a unique solution, nloose=-1 for no solution, nloose>0 for indeterminacy
#' @useDynLib gensysR
#' @export
#' @name gensysR_qzdetails
NULL
