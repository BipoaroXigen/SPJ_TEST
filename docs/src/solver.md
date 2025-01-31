# Simple solver examples

Here, we show, how a solver can be set up by providing it functions which define its behavior.

All of the examples are available in the `zoo.ipynb` notebook in the `/notebooks` directory.

---

## Single objective optimization

A simple solver for finding the minimum of 10-dimensional parabole.
Any other user defined function, with its arguments, can be chosen.


```
using EO

dims = 10
pop_size = 100
objective_function  = enclose_arguments(f_sphere, zeros(dims))  # min at origin
initialization      = enclose_noargs(interval_real_initialization, dims, pop_size, objective_function, 0, 100)
selection           = enclose_arguments(EO.s_tournament, pop_size, 3)
crossover           = enclose_arguments(EO.cr_parent_sum, pop_size)
mutation            = enclose_arguments(gaussian_preturbation!, 0.25)
replacement         = EO.enclose_replacement(EO.r_keep_best_n, pop_size)
termination         = enclose_argument(iteration_termination, pop_size*100)
solution            = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination)
```

## Multiple objectives

With a small modification, the solver can be edited to minimize two (or more) objectives at once.

```
pop_size = 100
dimension = 2

F = typeof(enclose_arguments(EO.f_weighted_sum, ones(dimension)))
objective_function  = EO.MultiObjFunction{F}([enclose_arguments(EO.f_sphere, ones(dimension)), EO.f_rastrigin])
initialization      = enclose_noargs(interval_real_initialization, dimension, pop_size, objective_function, -100, 100)
selection           = enclose_arguments(EO.s_identity, pop_size*2)
crossover           = cr_single_point
mutation            = enclose_arguments(gaussian_preturbation!, 0.25)
replacement         = EO.enclose_replacement(EO.r_NSGA, pop_size)
termination         = enclose_argument(iteration_termination, 10000)

solution = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination)
```

## Genetic programming

The solver can be also, purely by providing different functions to the solver, used for symbolic regression.

```
pop_size = 300

#objective_function  = enclose_arguments(EO.f_function_diff_subset, x, y, 0.9)
objective_function  = enclose_arguments(EO.f_function_diff_squared, y, x)
initialization      = enclose_noargs(EO.expression_initialization, pop_size, objective_function, basis_functions, basis_variables)
selection           = enclose_arguments(EO.s_tournament, 100, 3)
#selection           = enclose_arguments(EO.s_greedy_overselection, pop_size, 0.16)
crossover           = enclose_arguments(EO.cr_subtree, pop_size, basis_functions, basis_variables)
#crossover           = enclose_arguments(EO.cr_GSGP, pop_size, basis_functions, basis_variables)
mutation            = enclose_arguments(EO.subtree_mutation!, basis_functions, basis_variables)
replacement         = enclose_replacement(EO.r_best_n_diverse, pop_size, 0.7)       # chosen fraction of the population will consist of the best n, the rest is random
termination         = enclose_argument(iteration_termination, pop_size*20)

solution = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination)

model = EO.Expr_parser(solution.top_coords);
```

## Nested Evolutionary Feature Selection


To be able to regress much more complex functions I came up with my own algorith based on https://dl.acm.org/doi/10.1145/2739480.2754693.


```
variables = [domain]
#= operations = Vector{Function}([+, -, *, protected_div, sin, cos, square, cube, logaritmus, power])
arities =    [2,2,2,2,1,1, 1, 1, 1, 2] =#
operations = Vector{Function}([+, -, *, EO.protected_div, EO.sin, EO.square, EO.cube, EO.logaritmus])
arities =    [2,2,2,2,1,1,1,1]
exprs = Vector{Any}([:x])

best_model, models_hist = EO.feature_synthesis(y, variables, operations, arities, exprs, 10, q=10, μ=3, max_depth=10);
```