Optimally-weighted Herding is Bayesian Quadrature
=============

<img src="https://raw.githubusercontent.com/duvenaud/herding-paper/paper/figures/lovers2.png" width="300">

This repo contains the experiment source code as well as the latex source for a paper that mathematically relates two very different approaches to approximate integration.

Abstract:
---------

Herding and kernel herding are deterministic methods of choosing samples which summarise a probability distribution.  A related task is choosing samples for estimating integrals using Bayesian quadrature.  We show that the criterion minimised when selecting samples in kernel herding is equivalent to the posterior variance in Bayesian quadrature.  We then show that sequential Bayesian quadrature can be viewed as a weighted version of kernel herding which achieves performance superior to any other weighted herding method. We demonstrate empirically a rate of convergence faster than O(1/N).  Our results also imply an upper bound on the empirical error of the Bayesian quadrature estimate.

Contact the Authors:
------------------
[David Duvenaud](http://mlg.eng.cam.ac.uk/duvenaud/) (dduvenaud@seas.harvard.edu)

[Ferenc Huszar](http://mlg.eng.cam.ac.uk/ferenc/) (ferenc.huszar@gmail.com)

Running the code:
-------------------------
Running code/demo.m automatically reproduces most of the figures in the paper,
with some of the settings turned down to make the demo run fast.

If you want to exactly reproduce the results found in the paper, set

num_samples = 400;
num_queries = 10000;

The code used was optimized for legibility and simplicity, not speed.
Thus the herding implementation is O(N^3) instead of O(N^2), and the BQ implementation is O(N^4) instead of O(N^3).

Feel free to contact us if you have any questions.
