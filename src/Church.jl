module Church
using Distributions
include("If.jl")

export sample, condition, value, resample, gc_church, background_sampler, n_samples, pdf

#Types
type Sample{V}
    det::Any
    dist::Distribution
    value::V
    logp::Float64
    deps::Vector
    conditioned::Bool
end
type Det{T <: Tuple}
    f::Function
    args::T
    deps::Vector
    conditioned::Bool
end
type GetIndex{S, T <: Tuple}
    struct::S
    args::T
    deps::Vector
    conditioned::Bool
end
type Uneval
end
type NoCondition
end
const uneval = Uneval()
SDG = Union(Sample, Det, GetIndex)
SG = Union(Sample, GetIndex)
WSDG = Union(WeakRef, Sample, Det, GetIndex)

#Store declared samples/branches
const samples = Array(Union(Sample, WeakRef), 0)
#const nodes = Array(WSDG, 0)

#Show functions
import Base.show
show(io::IO, s::Sample) = print(io, "Sample", tuple(s.dist, s.value))
show(io::IO, b::GetIndex) = print(io, "GetIndex")

#Constructors
import Distributions.sample
sample(det) = begin
    dist = value(det)
    val = rand(dist)
    lp = logp(dist, val)
    s = Sample(det, dist, val, lp, Array(WSDG, 0), false)
    push!(samples, s)
#    push!(nodes, s)
    add_dep(s, det)
    s
end
condition(det, val) = begin
    dist = value(det)
    lp = logp(dist, val)
    s = Sample(det, dist, val, lp, Array(WSDG, 0), false)
    #push!(samples, s)
    #push!(nodes, s)
    add_dep(s, det)
    s
end
import Base.getindex
#Construct Det
getindex(f::Function, args...) = begin
    d = Det(f, args, Array(WSDG, 0), false)
#    push!(nodes, d)
    for arg in args
        add_dep(d, arg)
    end
    d
end
#Construct GetIndex
getindex(s, a1, a2, a3, a4::SDG, args...) = 
    GI(s, a1, a2, a3, a4, args...)
getindex(s, a1, a2, a3::SDG, args...) = 
    GI(s, a1, a2, a3, args...)
getindex(s, a1, a2::SDG, args...) = 
    GI(s, a1, a2, args...)
getindex(s, a1::SDG, args...) = 
    GI(s, a1, args...)
GI(struct, args...) = begin
    g = GetIndex(struct, args, Array(WSDG, 0), false)
    add_dep(g, struct[map(value, args)...])
    for arg in args
        add_dep(g, arg)
    end
    g
end

#Make the sample that you're currently creating a dependent of previous samples.
#Required by sample and condition.
add_dep(s::SDG, arg::SDG) =
    if !in(s, arg.deps)
        push!(arg.deps, s)
    end
add_dep(s::SDG, arg) = nothing

#Get the value of an expression.
value(det::Det) = begin
    det.f(map(value, det.args)...)
end
value(s::Sample) = s.value
value(s) = s
value(g::GetIndex) = begin
    res = g.struct[map(value, g.args)...]
    #Add g as a dependent of struct [arg]
    add_dep(g, res)
    value(res)
end

#Don't care whether cdf or pdf is needed.
logp(d::ContinuousDistribution, x) = logpdf(d, x)
logp(d::DiscreteDistribution, x) = logpmf(d, x)
#old_logp(s::Sample) = s.logp

#If can follow args back to source,
#deps_inner(s::Sample) = Sample[s]
#deps_inner(g::Union(GetIndex, Det)) = vcat(map(deps_inner, g.deps)...)
#deps_outer(s::Sample) = collect(Set(vcat(map(deps_inner, s.deps)...)...))
deps_inner(s::Sample, deps::Vector{Sample}) = begin
    if !in(s, deps)
        push!(deps, s)
    end
    nothing
end
deps_inner(s::Union(Det, GetIndex), deps::Vector{Sample}) =
    deps_recurse(s, deps)
deps_recurse(s, deps::Vector{Sample}) = begin
    for dep in s.deps
        deps_inner(dep, deps)
    end
    nothing
end

goes_back(origin::SDG) = (dependent::SDG) -> goes_back(origin, dependent)
goes_back(origin::SDG, dependent::Sample) = true
    #origin == dependent.dist
goes_back(origin::SDG, dependent::Det) = true
    #in(origin, dependent.args)
goes_back(origin::SDG, dependent::GetIndex) =
    in(origin, dependent.args) || origin == origin.struct[map(value, args...)]

#Propose, accept/reject.
resample_inner(s::Sample) = begin
    old_val  = s.value
    deps = Sample[]
    deps_recurse(s, deps)
    #deps = deps_outer(s)
    #s.deps = deps
    old_logp = mapreduce(dep->dep.logp, +, deps)
    s.value  = rand(s.dist)
    new_dists = map(s -> value(s.det), deps)
    new_logps  = map((d, s) -> logp(d, s.value), new_dists, deps)
    new_logp = sum(new_logps)
    if !(exp(new_logp - old_logp) > rand())
        #Reject change
        s.value = old_val
    else
        #Accept changes
        for i = 1:length(deps)
            deps[i].dist = new_dists[i]
            deps[i].logp = new_logps[i]
        end
    end
    nothing
end
resample() = 
    if length(samples) != 0
        index = rand(1:length(samples))
        resample_inner(samples[index])
    end

#Run sampling + garbage collection in background.
background_sampler() = 
    @async while true
        tic()
        resample()
        sleep(0.001)
#        if rand() < toq()
#            gc_church()
#        end
    end

#Sugar
#Two alternatives to allow construction of Distributions.
#Use getindex directly on the datatype, breaks array construction syntax!
#getindex{D<:Distribution}(dist::Type{D}, args...) = Det(dist, args)
#Define (and export) a function for each distribution!
for dist in filter!(isleaftype, subtypes(Distribution))
    tsym = symbol(string(dist))
    fsym = symbol(lowercase(string(dist)))
    eval(quote 
        $fsym(args...) = $tsym(args...)
    end)
    eval(parse("export $(string(fsym))"))
end

#Operators construct a Det if called with a sample.
for op in [+, -, *, /, \, .*, ./, .\]
    op(a::SDG, b::SDG) = op[a, b]
    op(a::SDG, b) = op[a, b]
    op(a, b::SDG) = op[a, b]
end

#GC related functions.
weaken(as::Vector) = begin
    for i = 1:length(as)
        as[i] = WeakRef(as[i])
    end
    nothing
end
strengthen(as::Vector) = begin
    filter!(x -> x != WeakRef(), as)
    for i = 1:length(as)
        as[i] = as[i].value
    end
    nothing
end
#weaken_deps(as::Vector) = begin
#    for s in as
#        for i = 1:length(s.deps)
#            if !s.deps[i].conditioned
#                s.deps[i] = WeakRef(s.deps[i])
#            end
#        end
#    end
#    nothing
#end
#strengthen_deps(as::Vector) = begin
#    for s in as
#        filter!(x -> x != WeakRef(), s.deps)
#        for i = 1:length(s.deps)
#            if isa(s.deps[i], WeakRef)
#                s.deps[i] = s.deps[i].value
#            end
#        end
#    end
#    nothing
#end
#weaken_branches() = begin
#    for b in branches
#        if bool(value(b.cond))
#            b.val_false = WeakRef(b.val_false)
#        else
#            b.val_true = WeakRef(b.val_true)
#        end
#    end
#    nothing
#end
#strengthen_branches() = begin
#    for branch in branches
#        if bool(value(branch.cond))
#            branch.val_false = strengthen_if(branch.val_false)
#        else
#            branch.val_true = strengthen_if(branch.val_true)
#        end
#    end
#    nothing
#end
gc_church() = begin
    gc_disable()
#    weaken_branches()
#    weaken_deps(samples); weaken_deps(branches)
    weaken(samples)
    gc_enable()
    gc()
    gc_disable()
    strengthen(samples)
#    strengthen_deps(samples); strengthen_deps(branches)
#    strengthen_branches()
    gc_enable()
end
#strengthen_if(wr::WeakRef) =
#    wr == WeakRef() ? uneval : wr.value
n_samples() = length(samples)
end
