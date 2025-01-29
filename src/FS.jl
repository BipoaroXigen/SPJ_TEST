using Lasso, Statistics

struct Model
    coefficients::Vector{Float64}
    expressions::Vector{Any}
    model
end

#= 
function Base.show(io::IO, m::Model)
    println(io, "coefficients: $(m.coefficients)")
    println(io, "expressions: $(m.expressions)")
end =#

function (m::Model)(x::Vector{Float64})
    res = zeros(length(x))
    for (i, e) in enumerate(m.expressions)
        res += EO.Expr_parser(e, Dict(:x => 1, :a => 2, :b => 3, :c =>4, :d => 5)).(x).*m.coefficients[i]
    end
    return res
end

function (m::Model)(x::Float64)
    res = 0.0
    for (i, e) in enumerate(m.expressions)
        if e isa Expr
            res += EO.Expr_parser(e, Dict(:x => 1, :a => 2, :b => 3, :c =>4, :d => 5))(x)*m.coefficients[i]
        end
    end
    return res
end

function generate_feature(current_features, operations, arities, current_exceptions, max_depth, importance)
    i = rand(1:length(arities))
    operation = operations[i]
    arity = arities[i]
    features = []
    exprs = []

    for input in 1:arity
        #feature_i = findmax((importance)[rand(1:length(importance), 2)])[2]
        feature_i = rand(1:length(current_features))
        push!(features, current_features[feature_i])
        push!(exprs, current_exceptions[feature_i])
    end

    res = operation.(features...)
    exp = :($(operation)($(exprs...)))

    #is_same = [res == feature for feature in current_features]

    cors = [ cor(feature, res) for feature in features]

    # is constants
    while all(x->x==res[1], res) || get_ast_len(exp) > max_depth || any(cors .>= 0.95)
        i = rand(1:length(arities))
        operation = operations[i]
        arity = arities[i]
        features = []
        exprs = []
    
        for input in 1:arity
            #feature_i = findmax((importance)[rand(1:length(importance), 2)])[2]
            feature_i = rand(1:length(current_features))
            push!(features, current_features[feature_i])
            push!(exprs, current_exceptions[feature_i])
        end
        res = operation.(features...)
        exp = :($(operation)($(exprs...)))

        cors = [ cor(feature, res) for feature in features]
    end

    # apply operation on features
    return res, exp
end


function feature_synthesis(y::Vector{Float64}, X::Vector{Vector{Float64}}, operations::Vector{Function}, arities::Vector{Int}, exprs::Vector{Any}, iters::Int; q=3, μ=3, max_depth=5)
    @assert length(arities) == length(operations) "the number of operations is different than the number of their arities"

    variables = X

    p = length(X)   # dimensionality of the problem
    M = p + q + μ

    # initialize features
    features = variables
    for i in 1:μ
        new_feature, e = generate_feature(features, operations, arities, exprs, max_depth, zeros(length(features)))
        push!(features, new_feature)
        push!(exprs, e)
    end

    e_min = Inf
    best_model = nothing
    #best_exprs = nothing
    models = Vector{Model}()
    i = 1
    while i < iters


        crossval_i = rand(1:length(y))

        F = reduce(hcat,features)
        #F[F.==Inf] .= 1e10

        # fit the model and evaluate it
        model = fit(LassoPath, F, y; cd_tol=1e10, criterion=:obj, algorithm=CovarianceCoordinateDescent) # features inputted as matrix
        push!(models, Model(Array(model.coefs[:, findmax(model.pct_dev)[2]]), exprs, model))

        λ_best_i = findmax(model.pct_dev)[2]
        ỹ = predict(model)[:, λ_best_i]
        e = (ỹ .- y).^2 |> sum
        if e < e_min
            #best_model = model
            best_model = Model(Array(model.coefs[:, findmax(model.pct_dev)[2]]), exprs, model)
            #best_exprs = copy(exprs)
            e_min = e
        end

        # pro kazdou lambdu, kdyz coef featury neni 0 procti pct_dev lambdy do importance featury
        α = model.coefs
        importance = zeros(length(features))
        for j in 1:length(model.pct_dev)
            useful = α[:, j] .!= 0
            importance[useful] .+= model.pct_dev[j]^2 
        end
        importance[1:p] .= maximum(importance)

        # fixed size model:
        order = sortperm(importance, rev=true)
        keep_i = order .<= p + q

        features = features[keep_i]
        exprs = exprs[keep_i]

        b = model.b0[λ_best_i]

        features = push!(features, ỹ)
        m_e = "model at $i"

        exprs = push!(exprs, m_e)

        # make the next generation of potential features
        for i in 1:μ
            new_feature, e = generate_feature(features, operations, arities, exprs, max_depth, importance)
            push!(features, new_feature)
            push!(exprs, e)
        end

        i += 1
    end

    return best_model, models
end


### HELPER FUNCTIONS ###

function extract_features(models::Vector{Model}, x::Vector{Float64})::Dict
    features = []
    for m in models
        for e in m.expressions
            if e in features
                continue
            end
            push!(features, e)
        end
    end

    fs = Dict()

    for feature in features
        fs[feature] = feature_eval(feature, models, fs, float.(x))
    end

    return fs
end

function feature_eval(f::Expr, models::Vector{EO.Model}, fs, x::Vector{Float64})
    operation = f.args[1]
    arguments = f.args[2:end]
    for (i, arg) in enumerate(f.args[2:end])
        if arg in keys(fs)
            arguments[i] = fs[arg]
        else
            value = f.args[1].(f.args[2:end])
            fs[arg] = value
            arguments[i] = value
        end
    end
    return operation.(arguments...)
end
function feature_eval(f::Symbol, models::Vector{EO.Model}, fs, x::Vector{Float64})
    return x
end
function feature_eval(f::String, models::Vector{EO.Model}, fs, x::Vector{Float64})
    i = parse(Int64, split(f)[3])
    model = models[i]
    res = zero(x)

    features = Vector{Vector{Float64}}()

    for (i, e) in enumerate(model.expressions)
        push!(features, fs[e])
    end
    λ_best_i = findmax(model.model.pct_dev)[2]
    res = predict(model.model, reduce(hcat,features))

    res[res .== Inf]  .= 1e10
    res[res .== -Inf] .= -1e10
    res[res .== NaN]  .= 0.0

    return res[:, λ_best_i]
end

function prediction(m::EO.Model, x::Vector{Float64}, fs::Dict; λ=1e3, α=1e3)::Vector{Float64}

    # get model variables values from computed features
    F = Vector{Vector{Float64}}()
    for e in m.expressions
        push!(F, fs[e])
    end

    # predict
    λ_best_i = findmax(m.model.pct_dev)[2]
    res = predict(m.model, reduce(hcat,F))

    # numerical stability
    res[res .>= λ] .= λ
    res[res .<= α] .= α
    #res[res .== NaN] .= 0.0

    # return the prediction with highest R^2
    return res[:, λ_best_i]
end

function simulate_trading(res::Vector{Vector{Vector{Float64}}}, price::Vector{Float64}, training_width::Int64, prediction_width::Int64, timestep::Int64)
    diffs = Vector{Float64}()

    balance = 1 # unit is money
    stock_balance = 0
    volatility = 0.3          # speed of fund moving, 1 means all cash moves from one asset to the other instantly
    required_certainty = 0
    λ = 0.1

    snapshots_s = Vector{Float64}()
    snapshots_m = Vector{Float64}()
    stock_value = Vector{Float64}()
    stock_value_decim = Vector{Float64}()
    actions = Vector{Float64}()
    buy  = zeros(length(res))
    sell = zeros(length(res))
    hold = zeros(length(res))

    smoothed_price = EO.signal_mean(float.(price))

    prev_diff = nothing

    i = 1
    for prediction in res
        push!(snapshots_m, balance)
        push!(snapshots_s, stock_balance)

        #update balance
        previous_price = #= smoothed_ =#price[(2+timestep*(i-1))]
        current_price = #= smoothed_ =#price[(2+prediction_width+training_width+timestep*(i-1))]

        change = current_price/previous_price
        push!(stock_value, stock_balance*current_price)
        push!(stock_value_decim, current_price)

        pr = EO.signal_mean(mean(prediction))
        last_known_val = pr[training_width]
        predicted_val  = pr[training_width+prediction_width]

        diff = 10*(predicted_val - last_known_val)
        fell = false
        #diff = current_price - previous_price
        if isnothing(prev_diff)
            prev_diff = diff
        else
            diff = (1-λ)*prev_diff + λ*diff
            fell = diff < prev_diff
            prev_diff = diff
        end
        push!(diffs, diff)
        if abs(diff) > required_certainty
            if diff < 0         # the value will go down
                #sell
                trade_amount = stock_balance*volatility
                stock_balance = stock_balance - trade_amount

                balance = balance + trade_amount*current_price

                push!(actions, -1)
                sell[i] = 1
            elseif diff > 0
                #buy
                trade_amount = balance*volatility
                balance = balance - trade_amount
                
                stock_balance = stock_balance + trade_amount/current_price
                
                push!(actions, 1)
                buy[i] = 1
            end
    #=         if fell
                #sell a little
                trade_amount = stock_balance*0.1
                stock_balance = stock_balance - trade_amount
                
                balance = balance + trade_amount*current_price
                
                push!(actions, -1)
                sell[i] = 1
            else
                #buy a little
                trade_amount = balance*0.1
                balance = balance - trade_amount
                
                stock_balance = stock_balance + trade_amount/current_price
                
                push!(actions, 1)
                buy[i] = 1
            end =#
        else
            push!(actions, 0)
            hold[i] = 1
        end
        i+=1
    end

    my_returns = (stock_value+snapshots_m)[end]/(stock_value+snapshots_m)[1]
    stock_returns = stock_value_decim[end]/stock_value_decim[1]

    return my_returns, stock_returns

end