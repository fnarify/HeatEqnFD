# HeatEqnFD

This program provides a solution to the 1-d heat equation via a finite difference scheme. The two schemes present are 
either an explicit Dirichlet or implicit Neumann scheme. The program plots the output data using [gnuplot 4.2](http://www.gnuplot.info/download.html) or greater, gnuplot should either be on your path  or, the command in line 77 of heat_eqn.c should be changed to include the path to gnuplot/bin.

The executable was compiled using mingw32 on Win7 64 bit.

![Example Dirichlet (alpha=1, nx=100, nt=1000, dt=0.0005)](exampleDir.png)
![Example Neumann (alpha=1, nx=100, nt=1000, dt=0.0005)](exampleNeu.png)

# Installing and Running

## *nix
Assuming you have homebrew installed you can:
~~~
git clone git://github.com/fnarify/HeatEqnFd
cd HeatEqnFd
make heat_eqn
./heat_eqn
~~~

## Note
Due to the wxWidgets dependency not always installing on a gnuplot install, you may need to install gnuplot with this command to get 'wxt' terminal variants:
~~~
sudo apt-get install gnuplot-x11
~~~
Otherwise, remove the line **set terminal wxt size 600,600** from both plotn and plotd.

## Windows
Easiest way would be to download the zip of the project off github by clicking 'Clone or download' near the top left corner. Then run the executable, otherwise you can recompile it with whatever you might use for compiling C programs on Windows like Cygwin or mingw.
