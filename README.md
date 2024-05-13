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
In its first version, the package works for models on hipergraphs where only one configuration in the factor nodes unsatisfies the interaction (K-SAT like). This contains as a particular case a model with pairwise interactions (2-SAT like) 

## How to provide a graph?

The package has its own implementation of hypergraph structures.

```julia 
mutable struct HGraph
    # This structure holds the hypergraph information
        N::Int64                                   # number of variable nodes
        M::Int64                                   # number of hyperedges
        K::Int64                                   # number of nodes in each hyperedge
        var_2_he::Array{Array{Int64, 1}, 1}        # list of hyperedges for each variable
        he_2_var::Matrix{Int64}                    # list of variables per hyperedge
        degrees::Vector{Int64}                     # list of variable nodes' degrees
        nodes_in::Array{Dict{Int64, Int64}, 1}     # dictionary with key=node_index and val=place_in_he
        #.......
end
```
that can be initialized with

```julia
build_HGraph(var_2_he::Array{Array{Int64, 1}, 1} , he_2_var::Matrix{Int64}, degrees::Vector{Int64})
```

There are a couple of implemented examples:

```julia
build_RR_HGraph(N::Int64, c::Int64, K::Int64, idum::Int64=rand(1:typemax(Int64)))          # Random Regular Graphs
build_ER_HGraph(N::Int64, c::Union{Float64, Int64}, K::Int64, idum::Int64=rand(1:typemax(Int64)))  # Erdös-Rényi
```
If the user does not provide a graph to the functions ```CME_KSAT``` or ```CDA_KSAT```, an empty graph is taken as default. The functions then expect to receive three numbers: $N$, $K$, and $\alpha$. These, together with the optional seed_g (seed for the random generator), are used to create an Erdös-Rényi hypergraph with size $N$, $K$ nodes per hyperedge and mean node connectivity $c=\alpha K$.

The couplings or links for the K-SAT boolean formula can be input via the parameter ```links::Matrix{Float64}```. The indexes are ```links[he, i]``` with $he=1, \ldots, M$ and $i=1, \ldots, K$, where $M$ is the number of hyperedges in the graph. If the user does not input any matrix of links, the default is to randomly choose one.

For a neat example see the end of section 

## How to provide a model?

The information about the interactions goes into the transition rates. These are related to the first three arguments of the functions ```CME_KSAT``` and ```CDA_KSAT```:

```julia 
        ratefunc::Function            # function that computes the transition rate
        rargs_cst                     # constant arguments that 'ratefunc' will receive when evaluated
        rarg_build::Function          # function that computes the non-constant arguments to be obtained
                                      # at each integration step and then passed to 'ratefunc'
```

The structure of ```ratefunc``` must fit the following

```julia 
        ratefunc(Ep::Int64, Em::Int64, rateargs...)
```
where ```Ep``` is the local energy when the variable is positive $(\sigma=1)$ and ```Em``` is the local energy when the variable is negative $(\sigma=-1)$. 


To build the rest of the arguments the user should use a function like 
```julia 
        rarg_build(graph::HGraph, st::State_CME, rargs_cst...)
        rarg_build(graph::HGraph, st::State_CDA, rargs_cst...)
```
The first argument is the graph, with all the associated information (see previous section).
The second argument is a ```struct``` with the following information:

```julia 
mutable struct State_CME
    p_cav::Array{Float64, 4}     # cavity conditional probabilities p_cav[hyperedge, site_cond, s_cond, chain]
    probi::Vector{Float64}       # single-site probabilities probi[node]
    pu_cond::Array{Float64, 3}   # cavity conditional probabilities of partially unsatisfied
                                 # hyperedges pu_cond[hyperedge, site_cond, s_cond]
    p_joint_u::Vector{Float64}   # probability of having the hyperedge unsatisfied p_joint_u[hyperedge]
end
```

```julia 
mutable struct State_CDA
    p_joint::Matrix{Float64}     # joint probabilities p_joint[hyperedge, chain]
    pu_cond::Array{Float64, 3}   # conditional probabilities of partially unsatisfied
                                 # hyperedges pu_cond[hyperedge, site_cond, s_cond]
    p_joint_u::Vector{Float64}   # probability of having the hyperedge unsatisfied p_joint_u[hyperedge]
end
```
All the other constant arguments of the transition rates must go into ``` rargs_cst... ```

The function ```rarg_build()``` must return a list of arguments ready to be inserted into the function ```ratefunc```

As an example, let us see the implementation of the Focused Metropolis Search algorithm for K-SAT optimization

```julia 
# This function computes the energy in a K-SAT formula, where there is only one unsatisfied
# configuration of the variables in a clause
function ener(p_joint_u::Vector{Float64})
    e = 0.0
    for he in eachindex(p_joint_u)
        e += p_joint_u[he]
    end
    return e
end

# This is the rate of FMS algorithm for KSAT
function rate_FMS_KSAT(Ep::Int64, Em::Int64, eta::Float64, avE::Float64, K::Number, N::Int64)
    dE = Em - Ep
    if dE > 0
        return Ep * N / K / avE * eta ^ dE
    else
        return Ep * N / K / avE
    end
end


# Each rate function needs some specific arguments. The general form of these functions 
function build_args_rate_FMS(graph::HGraph, st::State_CME, eta::Float64)
    avE = ener(st.p_joint_u)
    return eta, avE, graph.K, graph.N
end


# Each rate function needs some specific arguments. The general form of these functions 
function build_args_rate_FMS(graph::HGraph, st::State_CDA, eta::Float64)
    avE = ener(st.p_joint_u)
    return eta, avE, graph.K, graph.N
end
```

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
