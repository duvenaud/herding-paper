This folder contains most of the code used to create the figures in the paper:

Optimally-Weighted Herding is Bayesian Quadrature
by Ferenc Huszar and David Duvenaud


Running demo.m automatically reproduces most of the figures in the paper,
with some of the settings turned down to make the demo run fast.

If you want to reproduce the results found in the paper, set

num_samples = 400;
num_queries = 10000;

The code used was optimized for legibility and simplicity, not speed.
Thus the herding implementation is O(N^3)
and the BQ implementation is O(N^4)

If you want a faster version of the code, or have any questions at all
contact me at dkd23@cam.ac.uk.


David Duvenaud
April 2012

