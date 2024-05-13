# ApproxMasEq

[![Build Status](https://github.com/d4v1d-cub/ApproxMasEq.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/d4v1d-cub/ApproxMasEq.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Instalation

For now, the package is not an official Julia package. Therefore, the user should manually add it by running:

```julia
import Pkg
Pkg.add("https://github.com/d4v1d-cub/ApproxMasEq.jl.git")
```

## Description

This is a package to perform the numerical integration of some Approximated Master Equations for the dynamics of binary variables in graphs. It contains two different methods: the Cavity Master Equation (CME) [[1]](#1) [[2]](#2) [[3]](#3) and the Conditional Dynamic Approximation (CDA) [[4]](#4). Provided a model and some network representing the interactions, the package estimates the local probabilities of observing a specific configuration of the variables at time $t$.

The two main functions implemented in the package so far are ```CME_KSAT``` and ```CDA_KSAT```, with the following arguments
```julia 
        (ratefunc::Function, rargs_cst, rarg_build::Function;
         graph::HGraph=build_empty_graph(), N::Int64=0, K::Int64=0,
         alpha::Union{Float64, Int64}=0.0, seed_g::Int64=rand(1:typemax(Int64)),
         links::Matrix{Int8}=Matrix{Int8}(undef, 0, 0), seed_l::Int64=rand(1:typemax(Int64)), 
         tspan::Vector{Float64}=[0.0, 1.0], p0::Float64=0.5, method=Tsit5, eth::Float64=1e-6,
         cbs_save::CallbackSet=CallbackSet(), dt_s::Float64=0.1, abstol::Float64=1e-6, reltol::Float64=1e-3)
```

## How to provide a model?

The information about the interactions goes into the transition rates. These are related to the first three arguments of the functions ```CME_KSAT``` and ```CDA_KSAT```:

```julia 
        ratefunc::Function            # function that computes the transition rate
        rargs_cst                     # constant arguments that 'ratefunc' will receive when evaluated
        rarg_build::Function          # function that computes the non-constant arguments to be obtained
                                      # at each integration step and then passed to 'ratefunc'
```

The structure of 'ratefunc' must fit the following

```julia 
        ratefunc(Ep::Int64, Em::Int64, rateargs...)
```

The package implements an example and gives the user access to it. The  

## How to provide the graph?

## How to numerically integrate?

## Accessing the results

## References

<a id="1">[1]</a> 
E. Aurell, G. D. Ferraro, E. Domínguez, and R. Mulet. A cavity master equation for the continuous time dynamics of discrete spins models. Physical Review E, 95:052119, 2017.

<a id="2">[2]</a> 
E. Aurell, E. Domínguez, D. Machado and R. Mulet. Exploring the diluted ferromagnetic p-spin model with a cavity master equation. Physical Review E, 97:050103(R), 2018

<a id="3">[3]</a> 
E. Aurell, E. Domínguez, D. Machado and R. Mulet. A theory non-equilibrium local search on random satisfaction problems. Physical Review Letters, 123:230602, 2019

<a id="4">[4]</a> 
D. Machado, R. Mulet and F. Ricci-Tersenghi. Improved mean-field dynamical equations are able to detect the two-step relaxation in glassy dynamics at low temperatures. Journal of Statistical Mechanics: Theory and Experiment, 2023:123301
