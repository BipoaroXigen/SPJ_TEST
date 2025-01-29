function basic_preturbation!(x::BinaryChromosome, p::Real)
	indexes = rand(x.length) .< p
	x.sequence[indexes] = .!x.sequence[indexes]
end

function gaussian_preturbation!(x::RealChromosome, p::Real)
	change = randn(x.length).*p
	x.sequence = x.sequence .+ change
end

function gaussian_preturbation!(x::RealChromosome, p::Vector{<:Real})
	change = randn(x.length).*p
	x.sequence = x.sequence .+ change
end

### TSP

function weaklink_preturbation!(x::RealChromosome, G::Matrix{Float64})
	# find two longest edges and switch the ordering between
	a::Int = 1
	b::Int = 1
	max::Float64 = -Inf

	for i in 1:x.length-1
		d = G[x.sequence[i], x.sequence[i+1]]
		if d > max
			max = d
			b = a
			a = i+1
		end
	end
    a = a-1
    x.sequence[b:a] = shuffle(x.sequence[b:a])
end

function order_switch!(x::RealChromosome, G::Matrix{Float64})
	# find two longest edges and switch the ordering between
	idx = sort(rand(1:x.length, 2))
	a = idx[1]
	b = idx[2]
    x.sequence[a:b] = reverse(x.sequence[a:b])
end

function pair_switch!(x::RealChromosome, G::Matrix{Float64})
	# find two longest edges and switch the ordering between
	idx = sort(rand(1:x.length, 2))
	a = idx[1]
	b = idx[2]

	tmp = x.sequence[a]
	x.sequence[a] = x.sequence[b]
	x.sequence[b] = tmp
end

function subtree_mutation!(x::ExprChromosome, basis_functions::Vector{Function}, basis_variables::Vector{Any})
	#rand() >= 0.05 && return

	d = 2
	root = random_expression_grow(basis_functions, basis_variables, d, 0)	# AST with depths at most d

	change_node = rand() < 0.9

	if change_node
		subexps = get_subexpressions(x.sequence)
		subexp_cnt = length(subexps)
		if subexp_cnt > 1 && root isa Expr
			rand_subexp_i = rand(1:subexp_cnt)				# don't replace the root i think it could fuck up the tree
			x.sequence = replace_subexp!(x.sequence, rand_subexp_i, root)
		else
			leaves = flatten(get_leaves(x.sequence))
			rand_leaf_i = rand(1:length(leaves))
			x.sequence = replace_leaf!(x.sequence, rand_leaf_i, root)
		end
	else
		leaves = flatten(get_leaves(x.sequence))
		rand_leaf_i = rand(1:length(leaves))
		x.sequence = replace_leaf!(x.sequence, rand_leaf_i, root)
	end

	x.length = get_ast_len(x.sequence)
end

function no_mutation!(x::T, p) where T<:Chromosome
	return
end