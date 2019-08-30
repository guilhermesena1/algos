### GMM
This program clusters high dimension data points by estimating the means, covariances and
emission probabilities of gaussian distributions that maximize the observed data likelihood.

Input consists of a first line containing number of samples (r) and number of dimensions (c)
The next r lines contain c numbers separated by spaces or tabs

### compile:

```
make all
```

### run:
```
./gmm <gmm_3d_4clusters.in
```

The number of clusters in the algorithm is defined by the Bayesian 
Information Criterion ([BIC](https://projecteuclid.org/euclid.aos/1176344136)).
We set a high upper bound of 50 clusters and try all possibilities from 1 to 50.
Because each attempt is independent, it is an embarrassingly parallel problem.

This implementation uses openmp for paralellization when available. You
can define the number of threads used in it according to the openmp
documentation by defining the `OMP_NUM_THREADS` variable, for instance, to
run the algorithm with 8 threads:

```
export OMP_NUM_THREADS=8
```

if openmp is not installed, you can disable paralellization on compiling 
by disabling the opt flag:

```
make OPT=0
```
