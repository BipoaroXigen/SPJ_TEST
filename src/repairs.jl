"""repair chromosomes in population to be a feasible soluiton"""
function TSP_repair(p::Population, g::TSPgraph)::Population
    for x in p.population
        for i in 1:2:x.length
            # if no edge exists between x[i] and x[i+1] do something about it 
        end
    end
end