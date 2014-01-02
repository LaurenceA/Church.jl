using Church

dp(concentration::Real, base_measure::Function) = begin
    sticks = Mem(i::Int -> sample(beta(1., concentration)), Dict())
    atoms  = Mem(i::Int -> base_measure(), Dict())
    loop(i::Int) = 
        @If(sample(bernoulli[sticks[i]]), atoms[i], loop(i+1))
    () -> loop(1)
end

dp_mixture(concentration::Real, base_measure::Function, parameter::Function) = begin
    dp_ = dp(concentration, base_measure)
    () -> parameter(dp_())
end

d = dp(1., () -> sample(normal()))
ds = [d() for i = 1:20]

#f = dp_mixture(1., () -> sample(normal()), x -> sample(normal[x, 0.1]))
#fs = [f() for i = 1:20]
