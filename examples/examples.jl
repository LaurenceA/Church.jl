reload("church.jl")
background_sampler()

#Define Gaussian variables.
a = Normal[0, 1]
b = Normal[3, 1]

#Define a Gaussian variable, who's mean is the sum of a and b
d = Normal[a+b, 1]
#Condition d to be 10
condition(d, 10)

#Way to define a geometric distribution.
#Branch takes a {0,1} variable, and two lambdas.
#It only evaluates the lambda when necassery, so doesn't run off to infinity.
geom() = Branch(Bernoulli[], () -> (+)[1, geom()], ()->0)
g1 = geom()
g2 = geom()
g3 = geom()
#Condition the geometric distributions to have a sum of around 30.
#The geoms themselves are not distributions, so cannot be conditioned.
condition(Normal[g1+g2+g3, 1.], 30.)

#Demonstrate CRP etc.

