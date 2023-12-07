<!-- Title -->
<h1 align="center">
  GFN1-xTB.jl
</h1>

<!-- description -->
<p align="center">
  <strong>GFN1-xTB single point energy calculator.</strong>
</p>


GFN1-xTB.jl is the tight-biniding quantum chemistry code. <br>
More about the method: https://doi.org/10.1021/acs.jctc.7b00118. <br>
This project is done as the part of Remote Computatinoal Project course of the EMTCCM Master 2023. <br>
Author: Zarko Ivkovic, Univeristy of Barcelona, zivkoviv7 AT alumnes.ub.edu <br>
December 2023

## Installation instructions

Recommended Julia: Stable release v1.9.3 or newer

To install, download the `GFN1-xTB.jl`
[source](https://github.com/zarkoivkovicc/GFN1-xTB.jl) with:

```
$ git clone https://github.com/zarkoivkovicc/GFN1-xTB.jl.git
```

Now change into the `GFN1-xTB.jl` directory with

```
$ cd GFN1-xTB.jl
```

To use GFN1-xTB.jl, you need to instantiate all dependencies with:

```
$ julia --project
julia> ]
(GFN1-xTB) pkg> instantiate
```

Alternatively, if you want to use this code from the main julia environment, you can download dependencies manually:

```
$ julia
julia> ]
(@v1.9) pkg> add ArgParse Memoization
```

If you decide to use this option, you don't have to type --project every time you want to run GFN1-xTB.jl

## Running tests

Tests are stored in the 'test' folder. To run the tests, use the following command from the main directory

```
$ julia --project test/runtests.jl
```

