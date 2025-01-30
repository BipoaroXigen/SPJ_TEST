# EO

This is a Julia library for evolutionary optimization.
The library was written mainly with focus on flexibility of a solver, whose behaviour (the evolution strategy) is defined by user supplied functions as arguments (eg, crossover function etc.).

## Installation

Install the package by running `using Pkg; Pkg.add("https://github.com/BipoaroXigen/SPJ_TEST")`
or `]add https://github.com/BipoaroXigen/SPJ_TEST` in REPL.


## Examples of use

We provide several jupyter notebooks which help the user see how the library can be utilized.
The notebooks can be found in the directory `notebooks`.

### Solving the TSP problem

The notebook `TSP.ipynb` demonstrates the use and performance of the library when applied to the TSP problem.
The solved TSP problems are from the web page: http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/ .
The notebook expects the instances in a folder `TSP_problems/`.

### Multiobjective optimization

The notebook `multiobjective.ipynb` demonstrates the use and performance of the library when applied to constrained analytical optimization problems.
And the definitions of some problems are already prepared in the code.

### Optimization of constrained functions

The notebook `constrained_optimization.ipynb` demonstrates the use and performance of the library when applied to constrained analytical optimization problems.
The solved TSP problems are from the web page: https://cw.fel.cvut.cz/wiki/_media/courses/a0m33eoa/cviceni/2006_problem_definitions_and_evaluation_criteria_for_the_cec_2006_special_session_on_constraint_real-parameter_optimization.pdf.
And the definitions of some problems are already prepared in the code.

### Algorithmic trading

The folder `trading` contains notebooks of an algorithmic trading project, built with this library.
