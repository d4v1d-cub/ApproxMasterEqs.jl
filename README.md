# ApproxMasEq

[![Build Status](https://github.com/d4v1d-cub/ApproxMasEq.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/d4v1d-cub/ApproxMasEq.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Description

This is a package to perform the numerical integration of some Approximated Master Equations for the dynamics of binary variables in graphs. It contains two different methods: the Cavity Master Equation (CME) [[1]](#1) [[2]](#2) [[3]](#3) and the Conditional Dynamic Approximation (CDA) [[4]](#4). Provided a model and some network representing the interactions, the package estimates the local probabilities of observing a specific configuration of the variables at time $t$.

## How to provide a model?

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
