#ifndef __HEAT_EQN_H__
#define __HEAT_EQN_H__

// To help with portability of the program to multiple compilation environments.
#if (_WIN32 || __WIN32__)
    #define PRSIZET "Iu"
#elif (__STDC__ && __STDC_VERSION__ && (__STDC_VERSION__ >= 199901L))
    #define PRSIZET "zu"
#else
    #define PRSIZET "lu"
#endif

/**
 * Method for plotting the approximated data (using gnuplot, or whatever plotting software you prefer), 
 * from any of the finite difference methods below.
 *
 * output - name of file to write data output to.
 * script - gnuplot script to run to interpret data.
 * row - row dimension.
 * col - column dimension.
 * interior - whether only interior spatial points are included, e.g. spatial pts with values 0
 *            and col * dx are not interior pts.
 * dx - spatial step.
 * dt - timestep.
 * data - matrix of data to plot.
 */
void plot(char *output, char *script, size_t row, size_t col, int interior, double dx, double dt, double *data);

/**
 * Implicit and explicit finite difference schemes for solving the 1-d heat equation
 * with either Dirichlet or Neumann boundary conditions.
 * alpha - u_t = alpha * u_tt.
 * nx - number of spatial points.
 * nt - number of timesteps.
 * dt - timestep.
 */
void impheateqn_d(double alpha, size_t nx, size_t nt, double dt);
void expheateqn_d(double alpha, size_t nx, size_t nt, double dt);
void impheateqn_n(double alpha, size_t nx, size_t nt, double dt);
void expheateqn_n(double alpha, size_t nx, size_t nt, double dt);

#endif
