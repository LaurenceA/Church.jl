using Church
using Distributions

mem(f::Function) = begin
    d = Dict()
    nf(args...) = begin
        if !haskey(d, args)
            d[args] = f(args...)
        end
        d[args]
    end
    nf
end

dp(concentration::Real, base_measure::Function) = begin
    sticks = mem(i::Int -> sample(Beta(1., concentration)))  #Dict{Int, Sample}
    atoms  = mem(i::Int -> base_measure())                   #Dict{Int, Sample}
    loop(i::Int) = 
        @branch(sample(Bernoulli[sticks(i)]), atoms(i), loop(i+1))
    () -> loop(0)
end

dp_mixture(concentration::Real, base_measure::Function, parameter::Function) = begin
    dp_ = dp(concentration, base_measure)
    () -> parameter(dp_())
end

d = dp(1., () -> sample(Normal()))
ds = [d() for i = 1:20]

f = dp_mixture(1., () -> sample(Normal()), x -> sample(Normal[x, 0.1]))
fs = [f() for i = 1:20]

