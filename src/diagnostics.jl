resample(n::Int, skip::Int, args::RV...) = begin
    results = Vector[[value(arg)] for arg in args]
    for i = 1:n
        resample(skip)
        for i = 1:length(args)
            push!(results[i], value(args[i]))
        end
    end
    results
end
