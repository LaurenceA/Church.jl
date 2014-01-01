using Base.Test
using Church

#Geometric distribution

geom() = If(() -> 1+geom(), () -> 1)[sample(bernoulli())]
test_geom() = begin
    a = geom()
    condition(normal[a, 3], 15)
    iter=10^5
    as = zeros(iter)
    for i = 1:iter
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
likelihood = pdf(normal(15, 3), xs)
posterior = prior .* likelihood
norm_posterior = posterior/sum(posterior)
ma = sum(norm_posterior.*xs)
va = sqrt(sum(norm_posterior.*xs.^2) - ma^2)

println((ms, ma))
println((vs, va))

@test ma*0.95 < ms < ma*1.05
@test va*0.95 < vs < va*1.05

#gc_church()

##Normal distribution
#iter = 10^6
#as = zeros(iter)
#test_norm() = begin
#    a = sample(normal(0, 2))
#    b = condition(normal[a, 1], 3)
#    iter=10^6
#    as = zeros(iter)
#    for i = 1:iter
#        resample()
#        as[i] = value(a)
#    end
#    as
#end
#as = test_norm()
#ms = mean(as)
#vs = var(as)
#
##Analytic Moments
#va = 1/(1/2^2 + 1/1)
#ma = (3/1)/(1/2^2+1/1)
#
#@test 0.95*ma < ms < 1.05*ma
#@test 0.95*va < vs < 1.05*va
#
#as = [sample(normal()) for i = 1:2]
#ind1 = sample(bernoulli()) + 1
#ind2 = sample(bernoulli()) + 1
#b = as[ind1]
#c = as[ind2]
#
#d = b + c
