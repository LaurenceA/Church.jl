using Base.Test
using Church
using Distributions

#Geometric distribution

geom() = @branch(sample(Bernoulli()), 1+geom(), 1)
test_geom() = begin
    a = geom()
    condition(Normal[a, 3], 15)
    iter=10^6
    as = zeros(iter)
    for i = 1:10^6
        resample()
        as[i] = value(a)
    end
    as
end
as = test_geom()

ms = mean(as)
vs = sqrt(var(as))

#Analytic
n = 30
xs = 1:n
prior = 0.5.^xs 
likelihood = pdf(Normal(15, 3), xs)
posterior = prior .* likelihood
norm_posterior = posterior/sum(posterior)
ma = sum(norm_posterior.*xs)
va = sqrt(sum(norm_posterior.*xs.^2) - ma^2)

@test ma*0.95 < ms < ma*1.05
@test va*0.95 < vs < va*1.05

gc_church()

#Normal distribution
iter = 10^6
as = zeros(iter)
test_norm() = begin
    a = sample(Normal[0, 2])
    b = condition(Normal[a, 1], 3)
    iter=10^6
    as = zeros(iter)
    for i = 1:iter
        resample()
        as[i] = value(a)
    end
    as
end
as = test_norm()
ms = mean(as)
vs = var(as)

#Analytic Moments
va = 1/(1/2^2 + 1/1)
ma = (3/1)/(1/2^2+1/1)

@test 0.95*ma < ms < 1.05*ma
@test 0.95*va < vs < 1.05*va
