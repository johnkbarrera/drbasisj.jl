# Computation of the Demmler-Reinsch basis.

[![Build Status](https://travis-ci.org/ChrisRackauckas/ExamplePackage.jl.svg?branch=master)](https://travis-ci.org/ChrisRackauckas/ExamplePackage.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/9iuvdt0j0mw6au0k?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/examplepackage-jl)



The `Demmler-Reinsch basis` is provided for a given smoothness degree in a uniform grid. 
- `drbasis(nn,qq)`
- `nn`: Number of design points in the uniform grid. 
- `qq`: Smoothness degree of the basis.

## Details

The use of large numbers required by the basis is handled by the package Brobdingnag. The method assumes the grid is equidistant. Missing values are not supported.

## Value

A list object containing the following information.

- `eigenvalues`: estimadted eigenvalues.
- `eigenvectors`: estimated eigenvectors.
- `eigenvectorsQR`: orthonormal eigenvectors.
- `x`: equidistant grid used to build the basis.

## References

`Rosales F. (2016).` Empirical Bayesian Smoothing Splines for Signals with Correlated Errors: Methods and Applications
`Serra, P. and Krivobokova, T. (2015)` Adaptive Empirical Bayesian Smoothing Splines

## Authors
       Francisco Rosales
       John Barrera

## Installation

As described in the manual, to [install unregistered packages][unregistered], use `Pkg.clone()` with the repository url:

```julia
Pkg.clone("https://github.com/johnkevinbarrera/drbasisj.jl")
Pkg.clone("https://github.com/johnkevinbarrera/drbasisj.jl")
```

Julia version 0.4 or higher is required (install instructions [here][version]).


## Usage

To use the function `drbasis` has to invoke the package with `using drbasisj `, now you can use the function as follows `drbasisj.drbasis(nn,qq) `

## Example

plot elements of the basis, I use [Gadfly][gadfly]

```julia
using drbasisj
using Gadfly

n = 100
Basis=[drbasisj.drbasis(n,1),
	drbasisj.drbasis(n,2),
	drbasisj.drbasis(n,3),
	drbasisj.drbasis(n,4),
	drbasisj.drbasis(n,5),
	drbasisj.drbasis(n,6)]


#Eigenvalues

set_default_plot_size(12cm, 8cm)
p1 = plot(x=1:n, y=Basis[1][3], Geom.line, Guide.Title("Eigenvalues (q=1)")) 
p2 = plot(x=1:n, y=Basis[2][3], Geom.line, Guide.Title("Eigenvalues (q=2)"))
p3 = plot(x=1:n, y=Basis[3][3], Geom.line, Guide.Title("Eigenvalues (q=3)")) 
p4 = plot(x=1:n, y=Basis[4][3], Geom.line, Guide.Title("Eigenvalues (q=4)"))
p5 = plot(x=1:n, y=Basis[5][3], Geom.line, Guide.Title("Eigenvalues (q=5)")) 
p6 = plot(x=1:n, y=Basis[6][3], Geom.line, Guide.Title("Eigenvalues (q=6)"))

set_default_plot_size(24cm, 24cm)
gridstack([p1 p2 ; p3 p4; p5 p6])
```
![Match function](https://user-images.githubusercontent.com/7105645/28856524-38521768-7709-11e7-9cd9-204097fa09a2.png)

```julia
#Eigenvector for q=3

set_default_plot_size(12cm, 8cm)
r1 = plot(x=1:n, y=Basis[1][2][:,1+3], Geom.line, Guide.Title("Eigenvector n.4")) 
r2 = plot(x=1:n, y=Basis[2][2][:,2+3], Geom.line, Guide.Title("Eigenvector n.5")) 
r3 = plot(x=1:n, y=Basis[3][2][:,3+3], Geom.line, Guide.Title("Eigenvector n.6")) 
r4 = plot(x=1:n, y=Basis[4][2][:,4+3], Geom.line, Guide.Title("Eigenvector n.7")) 
r5 = plot(x=1:n, y=Basis[5][2][:,5+3], Geom.line, Guide.Title("Eigenvector n.8")) 
r6 = plot(x=1:n, y=Basis[6][2][:,6+3], Geom.line, Guide.Title("Eigenvector n.9")) 

set_default_plot_size(24cm, 24cm)
gridstack([r1 r2 ; r3 r4; r5 r6])
```
![Match function](https://user-images.githubusercontent.com/7105645/28857079-5a76c228-770c-11e7-86d0-653543d18bb5.png)

```julia

#example of a smooth function in the Demmler-Reinsch basis

n = 500
Basis=[drbasisj.drbasis(n,1),
	drbasisj.drbasis(n,2),
	drbasisj.drbasis(n,3),
	drbasisj.drbasis(n,4),
	drbasisj.drbasis(n,5),
	drbasisj.drbasis(n,6)]
	
coef3 = vcat(zeros(3),(pi*(2:(n-2))).^(-3.1)).*(cos(2*(1:n)))
A3 = Basis[3][1]
mu = -A3*coef3

set_default_plot_size(16cm, 16cm)
plot(x=1:n, y = mu, Geom.line, Guide.ylabel("mu"))
```
![Match function](https://user-images.githubusercontent.com/7105645/28857138-ba8d65ae-770c-11e7-91e4-892159f0a3d8.png)

## Note
Remember Julia not uses lists, so I will use matrices and vectors.
Now, if you have an `a = drbasis(nn,qq)` object, you can access the elements as follows:

- eigenvectors : a[1]
- eigenvectorsQR : a[2]
- eigenvalues : a[3]
- x : a[4]

However, in the R language it is different.

## Testing

In a Julia session, run `Pkg.test("drbasisj")`.

## keyword
	Demmler-Reinsch basis
	non-parametric

[unregistered]:http://docs.julialang.org/en/release-0.5/manual/packages/#installing-unregistered-packages
[version]:http://julialang.org/downloads/platform.html
[gadfly]:http://gadflyjl.org/stable/
