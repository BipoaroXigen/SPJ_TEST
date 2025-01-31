# Traveling salesman problem

Here, we show, how a solver can be set up to solve the Traveling salesman (TSP) problem.

All of the examples are available in the `TSP.ipynb` notebook in the `/notebooks` directory.
The solved TSP problems are from the web page: http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/ .
The notebook expects the instances in a folder `TSP_problems/`.

---

### Local search solver

The solver minimizes the length of a sequence of vertices of a graph G, which is defined by a distance matrix.


```
G = # distance matrix for the problem graph
pop_size = 1
dimension = size(G, 1)

objective_function  = enclose_arguments(EO.f_dist_sum, G)
initialization      = enclose_noargs(TSP_initialization, dimension, pop_size, objective_function)   # initialize random vertex sequences
selection           = EO.s_identity
crossover           = identity
mutation            = [enclose_arguments(EO.order_switch!, G), enclose_arguments(EO.pair_switch!, G), enclose_arguments(EO.weaklink_preturbation!, G)]
replacement         = EO.enclose_replacement(r_keep_best_n, pop_size)
termination         = enclose_argument(iteration_termination, 50000)    # maximal number of objective function calls
```

### Evolutionary strategy

Setting the population size above 1 and providing the solver with `selection` and `crossover` functions, the local search can be changed into an evolutionary algorithm.

```
G = # distance matrix for the problem graph
pop_size = 100
dimension = size(G, 1)

objective_function  = enclose_arguments(EO.f_dist_sum, G)
initialization      = enclose_noargs(TSP_initialization, dimension, pop_size, objective_function)
selection           = enclose_arguments(s_tournament, pop_size*3, round(Int, pop_size/3))
crossover           = [EO.cr_ordered, EO.cr_subtour, EO.cr_edge_recombination]
mutation            = enclose_arguments(EO.order_switch!, G)
replacement         = EO.enclose_replacement(r_keep_best_n, pop_size)
termination         = enclose_argument(iteration_termination, 5000#= 0 =#)
```

---

Showcase of memetic algorithm, heuristic initialization of initial population and results of benchmarks can be found in the `TSP.ipynb` notebook.