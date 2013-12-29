# Church

[![Build Status](https://travis-ci.org/LaurenceA/Church.jl.png)](https://travis-ci.org/LaurenceA/Church.jl)

Church.jl allows you perform inference on complex, and simple, probabilistic models.

Random variables
----------------
To create a random variable, we call`sample(dist)`, for instance,
```
using Distributions
using Church

#Construct random variable a, with a ~ N(0, 1).
a = sample(Normal(0, 1))

for i = 1:5
  #resample does one MCMC step
  resample()
  #value(a) returns the value of a in the current sample.
  println(value(a))
end

#Prints
```

Deterministic variables
-----------------------
You might also want to construct variables that are deterministic functions of random variables.
To do this, you call the deterministic function with square brackets,`abs[a]`, rather than brackets`abs(a).
The reason for this is straightforward, if `a = sample(Normal(0, 1))`, then `a::Sample`, and `abs` does not have a method that can deal with `Sample`.
The square brackets allow Church.jl to intercept the call and handle it correctly.
```
using Distributions
using Church

a = sample(Normal(0, 1))
#Finding the |a| requires us to use square brackets for the function call.
b = abs[a]
#Some operators can be used directly on samples
c = a*a

for i = 1:5
  #resample does one MCMC step
  resample()
  #value also works for deterministic functions.
  println(value(c))
end

#Prints
```

Conditioning
------------
Everything we have seen so far could be implemented using `rand`, for instance, we could sample `b` using `a = rand(Normal(0, 1))` and `b = abs(a).
Church.jl offers you the ability to condition random variables on observations.
For instance, to sample P(a, b| c=3) where c~Normal(b, 0.1),
```
using Distributions
using Church

a = sample(Normal(0, 1))
b = abs[a]
#Note that Normal doesn't know how to handle values of type Sample, 
i#so it must also be invoked using square brackets.
condition(Normal[b, 0.1])

#Now that we're doing inference we need a "burn-in" period,
#in which we allow the sampling process to reach equilibrium.
for i = 1:1000
  resample()
end

println(b)

#Prints
```

Constructing complex models
---------------------------
```
using Distributions
using Church

#Generate data
n = 30
data = 0.1*randn(n) + randn(n).^2

#Define the model parameters
sigma = sample(Gamma(1, 1))
m     = sample(Normal(0, 1))
v     = sample(Gamma(1, 1))

for i = 1:n
  l = sample(Normal[m, v])
  condition(Normal[l, sigma], data[i])
end

for i = 1:10^4
  resample()
end

println((value(sigma), value(m), value(v)))

#Prints
```


