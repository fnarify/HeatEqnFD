# HeatEqnFD

This program provides a solution to the 1-d heat equation via a finite difference scheme. The two schemes present are 
either an explicit Dirichlet or implicit Neumann scheme. The program plots the output data using [gnuplot 4.2](http://www.gnuplot.info/download.html) 
or greater, gnuplot should either be on your path  or line 74 in heat_eqn.c should be uncommented and the respective path for your 
gnuplot executable should be entered.

The executable was compiled using mingw32 on Win7 64 bit.

To compile this program on *nix systems, from line 254 onwards in heat_eqn.c the scanf parameter needs to be changed to the 
relative format for the compiler. So on C99 compatible compilers it can be changed from %Iu -> %zu.

![Example Dirichlet (alpha=1, nx=100, nt=1000, dt=0.0005)](https://github.com/fnarify/HeatEqnFD/blob/master/exampleDir.png)
![Example Neumann (alpha=1, nx=100, nt=1000, dt=0.0005)](https://github.com/fnarify/HeatEqnFD/blob/master/exampleNeu.png)
