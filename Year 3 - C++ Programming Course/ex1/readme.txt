Core Task:

I created a code that produces estimates of the 8-dimensional integral given in the manual, using a Monte Carlo technique. As expected, the error in the estimates falls off with increasing the number of Monte Carlo sample points. The exact behavior was investigated and from the gnuplot plot "MonteCarloPlot.png" it becomes evident that the error falls off as N^(-1/2). 

The code was able to produce a single estimate with N = 10^9 sample points in about 5 minutes. The estimate was 537.1866 +- 0.0012.

Using N = 10^7 the best estimate and error from nt = 25 estimates was 537.1853 +- 0.0119.

Supplementary Tasks:

The evaluations of the Cornu integrals were done using the class-based interface for the gsl_integration_qag routine in the cavlib header gslint_integ.hh. They were then used to plot the Cornu spiral and to investigate the Fresnel diffraction from a slit. The Cornu spiral plot is in the file "CornuSpiral.png". The relative amplitude plots for the cases of D=30, 50, 100cm are respectively in the files "AmplitudePlot30.png", "AmplitudePlot50.png", "AmplitudePlot100.png". The relative phase plots for the same cases are respectively in the files "PhasePlot30.png", "PhasePlot50.png", "PhasePlot100.png".
