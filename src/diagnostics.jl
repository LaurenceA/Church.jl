resample(n::Int, thin::Int, args::(RV...)) = begin
    results = map(arg -> typeof(value(arg))[value(arg)], args)
    for i = 1:iceil(n/thin)
        resample(thin)
        for i = 1:length(args)
            push!(results[i], value(args[i]))
        end
    end
    results
end
resample(n::Int, thin::Int, arg::RV) = resample(n, thin, (arg,))[1]
