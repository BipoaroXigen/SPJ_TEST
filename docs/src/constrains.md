# Solvers for constrained optimization

The EO library can be for example used to solve constrained optimizaation problems with the following algorithms.

The notebook `constrained_optimization.ipynb` demonstrates the use and performance of the library when applied to constrained analytical optimization problems.
The notebook contains more examples than this page.


### Stochastic ranking

```
pop_size = 100
objective(x) = (x.sequence[1]-10)^3 + (x.sequence[2]-20)^3
constraints = [EO.constraint(x -> -(x[1]-5)^2 - (x[2]-5)^2 + 100, .<=, 0.), EO.constraint(x -> (x[1]-6)^2 + (x[2]-5)^2 − 82.81, .<=, 0.),
					EO.constraint(x -> -x[1]+13, <=, 0.), EO.constraint(x -> x[1]-100, <=, 0.), 
					EO.constraint(x -> -x[2], <=, 0.), EO.constraint(x -> x[2]-100, <=, 0.)]

objective_function  = objective
initialization      = enclose_noargs(interval_real_initialization, dimension, pop_size, objective_function, 0, 100)
selection           = enclose_arguments(EO.s_stochastic_tournament, 30, 3, constraints, 0.4)
crossover           = enclose_arguments(EO.cr_parent_sum, pop_size)
mutation            = enclose_arguments(gaussian_preturbation!, 0.25)
replacement         = EO.enclose_replacement(EO.r_keep_best_n_stoch, pop_size, constraints, 0.4)
termination         = enclose_argument(iteration_termination, pop_size*100)

solution = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination; constraints=constraints)
```

### MOEA

Multi objective ES, where one objective dimension is the objective itself and the second objective is sum of all constraint violations.

This approach has a problem caused by the nature of the pareto front, where a solution with great objective value, but terrible constraint violation can be considered as dominating.

```
pop_size = 100
objective(x) = (x.sequence[1]-10)^3 + (x.sequence[2]-20)^3
constraints = [EO.constraint(x -> -(x[1]-5)^2 - (x[2]-5)^2 + 100, .<=, 0.), EO.constraint(x -> (x[1]-6)^2 + (x[2]-5)^2 − 82.81, .<=, 0.),
					EO.constraint(x -> -x[1]+13, <=, 0.), EO.constraint(x -> x[1]-100, <=, 0.), 
					EO.constraint(x -> -x[2], <=, 0.), EO.constraint(x -> x[2]-100, <=, 0.)]
f_pens(x::RealChromosome, constraints::Vector{EO.Constraint}) = sum(EO.Gs(x, constraints))

F = typeof(enclose_arguments(EO.f_weighted_sum, ones(dimension)))
objective_function  = EO.MultiObjFunction{F}([objective, enclose_arguments(f_pens, constraints)])
initialization      = enclose_noargs(interval_real_initialization, dimension, pop_size, objective_function, 0, 100)
selection           = enclose_arguments(EO.s_tournament, 30, 5)
crossover           = enclose_arguments(EO.cr_parent_sum, pop_size)
mutation            = enclose_arguments(gaussian_preturbation!, 0.25)
replacement         = EO.enclose_replacement(EO.r_NSGA, pop_size)
termination         = enclose_argument(iteration_termination, pop_size*100)

solution = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination)
```

### Adaptive penalization

The objective function is penalized by the sum of constraint violations, multiplied by rg(t) which controls the exploration exploitation tradeoff.

```

pop_size = 100
objective(x) = (x.sequence[1]-10)^3 + (x.sequence[2]-20)^3
constraints = [EO.constraint(x -> -(x[1]-5)^2 - (x[2]-5)^2 + 100, .<=, 0.), EO.constraint(x -> (x[1]-6)^2 + (x[2]-5)^2 − 82.81, .<=, 0.),
					EO.constraint(x -> -x[1]+13, <=, 0.), EO.constraint(x -> x[1]-100, <=, 0.), 
					EO.constraint(x -> -x[2], <=, 0.), EO.constraint(x -> x[2]-100, <=, 0.)]

objective_function  = objective
initialization      = enclose_noargs(interval_real_initialization, dimension, pop_size, objective_function, 0, 100)
selection           = enclose_arguments(EO.s_tournament, 30, 3)
crossover           = enclose_arguments(EO.cr_parent_sum, pop_size)
mutation            = enclose_arguments(gaussian_preturbation!, 0.25)
replacement         = EO.enclose_replacement(EO.r_keep_best_n, pop_size)
termination         = enclose_argument(iteration_termination, pop_size*100)

solution = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination; constraints=constraints, penalty=penatly);
```

### NSGA-II with modified tournament operator

The best of here implemented multi objective approaches. The violation of constraints is also added into the binary tournament, preffering less violating solutions. This algorithm mas much denser pareto front around the looked for optimum.

```

pop_size = 100
objective(x) = (x.sequence[1]-10)^3 + (x.sequence[2]-20)^3
constraints = [EO.constraint(x -> -(x[1]-5)^2 - (x[2]-5)^2 + 100, .<=, 0.), EO.constraint(x -> (x[1]-6)^2 + (x[2]-5)^2 − 82.81, .<=, 0.),
					EO.constraint(x -> -x[1]+13, <=, 0.), EO.constraint(x -> x[1]-100, <=, 0.), 
					EO.constraint(x -> -x[2], <=, 0.), EO.constraint(x -> x[2]-100, <=, 0.)]

f_pens(x::RealChromosome, constraints::Vector{EO.Constraint}) = sum(EO.Gs(x, constraints))

F = typeof(enclose_arguments(EO.f_weighted_sum, ones(dimension)))
objective_function  = EO.MultiObjFunction{F}([objective, enclose_arguments(f_pens, constraints)])
initialization      = enclose_noargs(interval_real_initialization, dimension, pop_size, objective_function, 0, 100)
selection           = enclose_arguments(EO.s_tournament, 30, 3)
crossover           = enclose_arguments(EO.cr_parent_sum, pop_size)
mutation            = enclose_arguments(gaussian_preturbation!, 0.25)
replacement         = EO.enclose_replacement(EO.r_cNSGA, pop_size, constraints)
termination         = enclose_argument(iteration_termination, pop_size*100)

solution = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination)

```