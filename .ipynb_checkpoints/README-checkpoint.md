# EO

[![Build Status](https://github.com/krizjona/EO.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/krizjona/EO.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Installation

The library is written in Julia. To run the interactive report you will need to install Julia 1.10 (`juliaup add 1.10` is the recommended way).

When running the library for the first time, run the commands `julia`, `]instantiate .`, `]activate .` in the main library directory, this will install all dependencies. 
To run the code any time after this calling `julia --project` will start julia and activate the installed dependencies.
To open the notebook environment, you can run it directly in VScode, or open it in browser by running `julia browse_notebooks.jl --project`.

## Homework 1

The interactive report is in the jupyter notebook file `TSP.ipynb`.
In case of problems with running the notebook an exported html version is also present.

If you wish to you can run the benchmark with the `run_benchmarks.ipynb` file, but the important results are already attached in the folder `benchmarks/TSP`.

## Homework 2

The interactive report is in the jupyter notebook file `report2.ipynb`.
In case of problems with running the notebook an exported html version is also present as `report2.html`.

If you wish to you can run the benchmark with the `hw2_benchmark.ipynb` file, but the important results are already attached in the folder `benchmarks/constr`.