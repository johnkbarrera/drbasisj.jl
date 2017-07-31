# drbasisj

[![Build Status](https://travis-ci.org/ChrisRackauckas/ExamplePackage.jl.svg?branch=master)](https://travis-ci.org/ChrisRackauckas/ExamplePackage.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/9iuvdt0j0mw6au0k?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/examplepackage-jl)
[![Coverage Status](https://coveralls.io/repos/github/ChrisRackauckas/ExamplePackage.jl/badge.svg?branch=master)](https://coveralls.io/github/ChrisRackauckas/ExamplePackage.jl?branch=master)


Creacion del repositorio

```julia
# Pkg.add("PkgDev")
using PkgDev
PkgDev.generate("drbasisprueba","MIT")

```

el codigo debe ir dentro de un modulo siempre


## Instalando el Pkg drbasisprueba

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
