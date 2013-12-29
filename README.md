# Church

[![Build Status](https://travis-ci.org/LaurenceA/Church.jl.png)](https://travis-ci.org/LaurenceA/Church.jl)

Chuch.jl aims to make it easy for anyone to perform inference in complex, and simple, probabilistic models.

Random variables
----------------
You can construct a random variable by calling `sample` on a distribution, from the `Distributions` package, for example,
```julia
using Distributions
using Church

#Define the random variable a, with distribution N(0, 1).
a = sample(Normal(0, 1))

for i = 1:5
  #Resample performs a single MCMC step.
  resample()
  #value(a) gives the value of a for the current sample.
  println(value(a))
end

#Prints:
#-0.1029734349357338
#-1.015922485444797
#-2.461967392570218
# 0.7876594234354186
#-1.6596582450563369
```

Deterministic functions
-----------------------
You can construct deterministic functions of random variables by replacing the brackets in the function call by square brackets.  
For instance, you use `abs[a]` instead of `abs(a)`.
The reason for this is straightforward, `a` has type `Sample`, so `abs` does note know how to deal with it directly.
Calling `abs[a]` allows Church.jl to handle the function call correctly.
Extending the model we wrote previously, we could use,
```julia
using Distributions
using Church

a = sample(Normal(0, 1))
b = abs[a]

for i = 1:5
  resample()
  #Note that value also works for deterministic functions of random variables.
  println(value(b))
end

#Prints:
# 0.8135278459720077
# 0.5929276217773704
# 0.2504229036079722
# 0.20016542444988056
# 1.998632465314105
```

Conditioning
------------
So far, we haven't done anything interesting - you can sample `a` and `b` in the previous sections by simply using `a = rand(Normal(0, 1))` and `b = abs(a)`.
Church.jl is interesting because you can condition these draws on known data.
For instance, to sample P(a, b| c=10), where c ~ Normal(b, 0.1),
```julia
using Distributions
using Church

a = sample(Normal(0, 1))
b = abs[a]
#Note that Normal does not know how to handle b::Sample, so Normal must be invoked with square brackets.
condition(Normal[b, 0.1], 3)

#Now that we're doing inference, we need to perform many sampling steps, 
#for the model to converge to the correct distribution.  This is known as burn-in.
for i = 1:1000
  resample()
end
println(value(a))

#Prints:
#-2.7382930822004345
```

Constructing more complex models
--------------------------------
One very important aspect of the Church.jl framework is that you can use arbitrary programming constructs to generate your model, for instance,
```julia
using Distributions
using Church

#Generate some data
data = randn(10)

#Define the model parameters
m = sample(Normal(0, 10))
v = sample(Gamma(2, 2))

for i = 1:10
  #Condition on the data.
  condition(Normal[m, v], data[i])
end

for i = 1:1000
  resample()
end
println((value(m), value(v)))

#Prints:
#(1.1352555080470663,0.8533067108422825)
```

If statements
-------------
If statements are useful, as they allow model selection
```julia
using Distributions
using Church

#Generate data.
data = randn(10)

#The distribution could be a Normal, or a Gamma.
dist = @branch(sample(Bernoulli()), Normal(0, 1), Gamma(1, 1))

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
