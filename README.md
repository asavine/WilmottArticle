# Computation graphs, AAD and back-propagation
Companion code to the series of articles in Wilmott Magazine 

Antoine Savine, 2019

<p align="center"> 

![Wilmott November 2019](WilmottNov2019.jpg)

</p>

funWithGraphs.h contains the code snippets from Part I 

BlackScholes.h contains an implementation of the Black-Scholes formula. It relies on gaussians.h, which contains classic implementations of the Cumulative Normal Distribution and its inverse.

AAD.h contains the AAD framework developed in part II

dupireBarrier.h contains the pricing and risk code of part III. It relies on a number of utilities: matrix.h (a simple adapter class wrapping a vector with a matrix view) and interp.h (one and two dimensional linear and smooth-step interpolation). It also relies on random number generators, with base class written in random.h and two concrete implementation: L'Ecuyer's MRG32K3A (mrg32k3a.h) and Sobol (sobol.cpp and sobol.h).

The Excel files xl*.* implement the export of C++ functions to excel, as documented in the tutorial https://github.com/asavine/xlCppTutorial

In addition, the repo contains a prebuilt xll wArticle.xll, built with Visual Studio 2017 from the project file wArticle.vcxproj, a demonstration spreadsheet xlTest.xlsx and two redistributable installers that may be necessary to run the xll on a machine without an installation of Visual Studio.
