# Church

[![Build Status](https://travis-ci.org/LaurenceA/Church.jl.png)](https://travis-ci.org/LaurenceA/Church.jl)

Chuch.jl aims to make it easy for anyone to perform inference in complex, and simple, probabilistic models.

Constructing a model
----------------
This is easiest to describe by example.
```julia
using Church

#You can define a random variable, simply by calling a function named according to the distribution.
a = normal()

#You can apply usual arithemetic operators to samples.
#Here we square of a.
b = a*a

#To apply other functions to random variables, you must first "lift" the function,
#so the function knows how to deal with random variables.
#The second argument to lift is the number of input arguments.
@lift(cosh, 1)
c = cosh(a)

#You can define random variables which depend on previous random variables,
#Here, we define a Gaussian with mean a, and stdev cosh(a)
d = normal(a, c)

for i = 1:5
  for j = 1:5
    #Resample performs a single MCMC step.
    resample()
  end
  #value(a) returns the value of a for the current sample.
  @printf("a:% .3f, b:% .3f, c:% .3f, d:% .3f", value(a), value(b), value(c), value(d)); println()
end

#Prints:
#a: 0.978, b: 0.957, c: 1.518, d:-0.240
#a: 0.262, b: 0.069, c: 1.035, d:-1.598
#a: 0.262, b: 0.069, c: 1.035, d:-1.065
#a: 0.262, b: 0.069, c: 1.035, d:-0.458
#a:-0.182, b: 0.033, c: 1.017, d:-0.112
```

Conditioning
------------
So far, we haven't done anything interesting - you could sample `a` and `b` in the previous sections by simply using `a = rand(Normal(0, 1))` and `b = a*a`.
In Church.jl, you can condition these draws on known data.
For instance, to sample P(a, b| c=10), where c ~ Normal(b, 0.1),
```julia
using Church

@lift(abs, 1)

a = normal(0, 1)
b = abs(a)

#Lets say we know that c=3 was drawn from a Gaussian with mean b, and stdev 0.1.
c = normal(b, 0.1; condition=3)

#Now that we're doing inference, we need to perform many sampling steps, 
#for the model to converge to the correct distribution.  This is known as burn-in.
for i = 1:1000
  resample()
end
println(value(a))

#Prints:
#-2.7382930822004345
```

Mixture Model
-------------
What if we don't know the
```julia
using Church

#Generate some data
data = [randn(10)+6, randn(10)-6]

#The model parameters
K = 2
ms = [normal(0, 10) for i = 1:K]
vs = [gamma(2, 2) for i = 1:K]
ps = dirichlet([1.,1.])

#Which mixture component does each data item belong to?
ks = [categorical(ps) for i = 1:length(data)]

for i = 1:length(data)
  #Condition on the data.
  normal(ms[ks[i]], vs[ks[i]]; condition=data[i])
end

for i = 1:10^4
  resample()
end
@printf("m1:% .3f, m2:% .3f, v1:% .3f, v2:% .3f", value(ms[1]), value(ms[2]), value(vs[1]), value(vs[2]))
println()
map(x -> print(value(x)), ks)
println()
println((value(ps)[1], value(ps)[2]))

#Prints:
#m1:-5.575, m2: 6.024, v1: 0.947, v2: 0.940
#22222222221111111111
#(0.43393275606877524,0.5660672439312248)
```

Varying the number of components.
If statements
-------------
If statements are useful, as they allow model selection
```julia
using Distributions
using Church

#Generate data.
data = randn(10)

#The distribution could be a Normal, or a Gamma.
dist = @If(bernoulli(), normal(0, 1), normal(0, 2))

for i = 1:10
  condition(dist, data[i])
end

for i = 1:5
  resample()
  println(value(dist))
end

#Prints:
Normal( μ=0.0 σ=1.0 )
Normal( μ=0.0 σ=1.0 )
Normal( μ=0.0 σ=1.0 )
Normal( μ=0.0 σ=1.0 )
Normal( μ=0.0 σ=1.0 )
```
