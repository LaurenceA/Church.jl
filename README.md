# Church

[![Build Status](https://travis-ci.org/LaurenceA/Church.jl.png)](https://travis-ci.org/LaurenceA/Church.jl)

Chuch.jl aims to make it easy for anyone to perform MCMC inference in complex, and simple, probabilistic models.
We aim to be:
 - Fast
   - Church.jl considers only the bits of the model that need to be updated during MCMC.
 - Practical
   - Easy to install and use.
 - General
   - Define new, arbitrarily complex distributions, within Church.jl.
   - Allows you to define stochastic, recursive, and non-parametric distributions.
   - Uses all distributions available in Distributions.jl.
   - Includes a probabilistic garbage collector - so you define large models without running out of memory.

Getting started
---------------
To install, use
```julia
Pkg.add("Church")
```
To load, use
```julia
using Church
```

Constructing a model
---------------------
You can define a random variable simply by calling a function,
```julia
a = normal()
```
Church.jl contains all the distributions in Distributions.jl - with function names that are just lowercase versions of the distribution names.

You can apply standard operators (e.g. `+`) directly to samples,
```julia
b = normal() * normal()
```

To apply other functions to random variables, you must first "lift" the function.
This overloads the function, so that it can deal with random variables.
For instance,
```julia
@lift(cosh, 1)
c = cosh(normal())
```
Note that the second argument to lift is the number of arguments, which, in this case, is 1.

Combining these, we could write,
```julia
a = normal()
b = a * a
@lift(cosh, 1)
c = cosh(a)
d = normal(a, c)
```
To sample these random variables, we use `resample()`, which performs a single MCMC step, and use `value(a)` to report the value of `a` for the current sample.
Note that `value(a)` is ONLY provided for recording or printing the value of samples.
Any other is liable to give meaningless quantities in the best case, or cause the algorithm to no longer sample the correct distribution in the worst case.
```julia
for i = 1:5
  #Do 5 MCMC steps.
  resample(5)
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
In Church.jl, you can condition these draws on known data, using the keyword argument `condition`,
```julia
normal(1, 1; condition=3)
```
In a more complete example, to sample P(a| c=10), where c ~ Normal(b, 0.1), and a ~ Normal(0, 1),
```julia
using Church

@lift(abs, 1)

a = normal(0, 1)
b = abs(a)
c = normal(b, 0.1; condition=3)

#Now that we're doing inference, we need to perform many sampling steps, 
#for the model to converge to the correct distribution.  This is known as burn-in.
resample(10^3)
println(value(a))

#Prints:
#-2.7382930822004345
```

Examples: mixture model - fixed number of components.
-------------------------------
```julia
using Church

#Generate some data
data = [randn(10)+6, randn(10)-6]

#The model parameters
K = 2
ms = [normal(0, 10) for i = 1:K]
vs = [gamma(2, 2) for i = 1:K]
ps = dirichlet(ones(K))

#Which mixture component does each data item belong to?
ks = [categorical(ps) for i = 1:length(data)]

for i = 1:length(data)
  #Condition on the data.
  normal(ms[ks[i]], vs[ks[i]]; condition=data[i])
end

resample(10^4)
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
Note that the inferred parameters are sensible, given the data.

Garbage collecting unused mixture components
------------------------------------------
We might want to define a mixture model with a variable number of components, for instance,
```julia
K = poisson(3)
ms = [normal(0, 10) for i = 1:K]
```
However, you cannot do this, because the list comprehension needs K to be an integer, not a sample.
Instead, you can use a large number of components, then exploit the lazy datastructures and garbage collector in Church.jl to avoid instantiating unused mixture components.
For instance,
```julia
using Church
using Distributions

#Generate some data
data = [randn(10)+6, randn(10)-6]

#The model parameters
ms = Mem((i::Int) -> normal(0, 10))
vs = Mem((i::Int) -> gamma(2, 2))
ps = dirichlet(10, 1.; sampler=(d,v)->Dirichlet(3*v+0.001))

#Which mixture component does each data item belong to?
ks = [categorical(ps) for i = 1:length(data)]


#Condition on the data.
for i = 1:length(data)
    normal(ms[ks[i]], vs[ks[i]]; condition=data[i])
end

for i = 1:10^3
  resample(10^3)
  gc_church()
end
map(x -> print(value(x)), ks)

#Prints:
#89998788884444444444
```
So the model is only using 4 components.
Looking at the value of `ms`, we see that the parameters for the other components have not been instansiated.
The other components, that have been created at some point during the sampling, have been cleaned up.

Defining new distributions
--------------------------
In Church.jl, a distribution is just a function, so to define a new distributioon, we just need to define a function.  A very simple example, is a mixture of normal distributions with different standard deviations,
```julia
gsm(m::Real) = normal(m, gamma(2, 2))
```
We can also allow conditioning on new distributions, if the final call also allows conditioning.  In this case,
```julia
gsm(m::Real; condition=nocond) = normal(m, gamma(2, 2); condition=condition)
```
Now `gsm` can be conditioned just like any other distribution.  Note that you should use nocond as the default value of condition - this is the special value indicating that the distribution is not conditioned.

In another example, we could use recursion to write down a distribution,
```julia
geom(p::Real) = @If(bernoulli(p), 1+geom(), 1).
```
The macro `@If` returns `1+geom()` if `bernoulli(p)` is `1`, and returns `1` if `bernoulli(p)` is `0`.
It only evaluates its arguments as necassery, so we do not get infinite recursive calls to `geom`.

Finally, you can write down random distributions.
For instance, the dirichlet distribution can be thought of as returning a vector, whose elements are positive and sum to 1, or it can be thought of as returning a categorical distribution.
The `dirichlet` distribution returns a vector, in Church.jl.
However, we could define `fdirichlet`, which does return a distribution,
```julia
fdirichlet(args...) = begin
    ps = dirichlet(args...)
    () -> categorical(ps)
end

#dir is a categorical distribution
dir = fdirichlet(9, 0.1)

for i = 1:10
    print(value(dir()))
end

#Prints
#4444448444
```

This mechanism allows you to write down a dirichlet process, which, again, is a distribution over random distributions,
```julia
dp(concentration::Real, base_measure::Function) = begin
    sticks = Mem(i::Int -> beta(1., concentration))
    atoms  = Mem(i::Int -> base_measure())
    loop(i::Int) = 
        @If(bernoulli(sticks[i]), atoms[i], loop(i+1))
    d = () -> loop(1)
end

#d is a DP
d = dp(1., normal)

for i = 1:10
    println(value(d()))
end

#Prints
#-0.8241868921220118
#1.2419259288985225
#0.7481036918602126
#1.4163485859423797
#1.2419259288985225
#1.2419259288985225
#1.2419259288985225
#1.2419259288985225
#1.2419259288985225
#1.2419259288985225
```
