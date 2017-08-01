# Computation of the Demmler-Reinsch basis.

[![Build Status](https://travis-ci.org/ChrisRackauckas/ExamplePackage.jl.svg?branch=master)](https://travis-ci.org/ChrisRackauckas/ExamplePackage.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/9iuvdt0j0mw6au0k?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/examplepackage-jl)



The `Demmler-Reinsch basis` is provided for a given smoothness degree in a uniform grid. 
- `drbasis(nn,qq)`
- `nn`: Number of design points in the uniform grid. 
- `qq`: Smoothness degree of the basis.

### Details

The use of large numbers required by the basis is handled by the package Brobdingnag. The method assumes the grid is equidistant. Missing values are not supported.

### Value

A list object containing the following information.

- `eigenvalues`: estimadted eigenvalues.
- `eigenvectors`: estimated eigenvectors.
- `eigenvectorsQR`: orthonormal eigenvectors.
- `x`: equidistant grid used to build the basis.

### References

`Rosales F. (2016).` Empirical Bayesian Smoothing Splines for Signals with Correlated Errors: Methods and Applications
`Serra, P. and Krivobokova, T. (2015)` Adaptive Empirical Bayesian Smoothing Splines

### Authors
Francisco Rosales
John Barrera


## Intallation

Since this package is not registered, you must install it by cloning. To add this package, use:

```julia
Pkg.clone("https://github.com/johnkevinbarrera/drbasisprueba.jl")
```

## Trabajando con el pkg drbasisprueba

for use this package, yo should do:

```julia

using drbasisprueba


drbasisprueba.drbasis(500,6)
```
