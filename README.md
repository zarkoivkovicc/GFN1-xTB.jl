<!-- Title -->
<h1 align="center">
  GFN1-xTB.jl
</h1>

<!-- description -->
<p align="center">
  <strong>GFN1-xTB single point energy calculator.</strong>
</p>


GFN1-xTB.jl is the tight-biniding quantum chemistry code. 
More about the method: https://doi.org/10.1021/acs.jctc.7b00118.
This project is done as the part of Remote Computatinoal Project course of the EMTCCM Master 2023.
Author: Zarko Ivkovic, Univeristy of Barcelona, zivkoviv7 AT alumnes.ub.edu
December 2023

## Installation instructions

Recommended Julia: Stable release v1.9.4

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
(GFN1_xTB.jl) pkg> instantiate
```

## Running tests

Tests are stored in the 'test' folder. By running runtests.jl from the corresponding test environment, all test will be automatically run.

```
$ cd test
$ julia --project=test runtests.jl
```
