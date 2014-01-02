module Church
using Distributions

export sample, condition, value, resample, gc_church, background_sampler, n_samples, pdf

#Types
type Sample{V}
    det::Any
    dist::Distribution
    value::V
    logp::Float64
    deps::Vector
end
type Det{T <: Tuple}
    f::Function
    args::T
    deps::Vector
end
type GetIndex{S, T <: Tuple}
    struct::S
    args::T
    deps::Vector
end
SDG = Union(Sample, Det, GetIndex)
SD = Union(Sample, Det)
WSDG = Union(WeakRef, Sample, Det, GetIndex)

#Store declared samples/branches
const samples = Array(Union(Sample, WeakRef), 0)
const dets = Array(Union(Det, WeakRef), 0)
const getindexes = Array(Union(GetIndex, WeakRef), 0)
const conditions = Array(Union(Sample, WeakRef), 0)

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
    s = Sample(det, dist, val, lp, Array(WSDG, 0))
    push!(samples, s)
    add_dep(s, det)
    s
end
condition(det, val) = begin
    dist = value(det)
    lp = logp(dist, val)
    s = Sample(det, dist, val, lp, Array(WSDG, 0))
    add_dep(s, det)
    s
end
import Base.getindex
#Construct Det
getindex(f::Function, args...) = begin
    d = Det(f, args, Array(WSDG, 0))
    push!(dets, d)
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
getindex(s, a1::SDG) = 
    GI(s, a1)
GI(struct, args...) = begin
    g = GetIndex(struct, args, Array(WSDG, 0))
    push!(getindexes, g)
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

deps_inner(_, s::Sample, deps::Vector{Sample}) = begin
    if !in(s, deps)
        push!(deps, s)
    end
    nothing
end
deps_inner(_, s::Det, deps::Vector{Sample}) =
    deps_recurse(s, deps)
deps_inner(origin, s::GetIndex, deps::Vector{Sample}) =
    if in(origin, s.args) || (origin == s.struct[map(value, s.args)...])
        deps_recurse(s, deps)
    end
deps_recurse(s, deps::Vector{Sample}) = begin
    for dep in s.deps
        deps_inner(s, dep, deps)
    end
    nothing
end

#check_dep(origin::SDG) = (dependent) -> check_dep(origin, dependent)
#check_dep(origin::SDG, dependent::Union(Sample, Det)) = begin
#    @assert in(origin, depdendent.args)
#    true
#end
#check_dep(origin::SDG, dependent::GetIndex) =
#    in(origin, dependent.args) || (origin == dependent.struct[map(value, dependent.args)...])
remove_deps(origin::SDG) =
    filter!(x -> isdependent(origin, x), origin.deps)
remove_deps(as::Vector) = begin
    for a in as
        remove_deps(a)
    end
    nothing
end

#Propose, accept/reject.
resample_inner(s::Sample) = begin
    old_val  = s.value
    deps = Sample[]
    deps_recurse(s, deps)
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
        sleep(1E-6)
        if rand() < toq()
            gc_church()
        end
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
weaken_store(as::Vector) = begin
    for i = 1:length(as)
        as[i] = WeakRef(as[i])
    end
    nothing
end
strengthen_store(as::Vector) = begin
    filter!(x -> x != WeakRef(), as)
    for i = 1:length(as)
        as[i] = as[i].value
    end
    nothing
end
weaken_deps(sdgs::Vector) = begin
    for s in sdgs
        for i = 1:length(s.deps)
            s.deps[i] = WeakRef(s.deps[i])
        end
    end
    nothing
end
strengthen_deps(sdgs::Vector) = begin
    for s in sdgs
        filter!(x -> x != WeakRef(), s.deps)
        for i = 1:length(s.deps)
            if isa(s.deps[i], WeakRef)
                s.deps[i] = s.deps[i].value
            end
        end
    end
    nothing
end
weaken_getindexes() = begin
    for g in getindexes
        weaken(g.struct)
    end
    nothing
end
strengthen_getindexes() = begin
    for g in getindexes
        strengthen(g.struct)
    end
    nothing
end

wr_value(wr::WeakRef) = wr.value
wr_value(sr) = sr

isdependent(parent::SDG, child::Sample)   = parent == child.det
isdependent(parent::SDG, child::Det)      = in(parent, child.args)
isdependent(parent::SDG, child::GetIndex) = 
    in(parent, child.args) || (parent == child.struct[map(value, child.args)...])

args(s::Sample) = s.det
args(d::Det) = d.args
args(g::GetIndex) = 
    SDG[g.args, (parent == g.struct[map(value, g.args)...])]

condition_rec(parent::SDG, child::SDG) = begin
    @assert isdependent(parent, child.args)
    i = findfirst(x -> x==wr_value(child), deps)
    if isa(parent.deps[i], WeakRef)
        parent.deps[i] = parent.deps[i].value
        for g_parent in args(parent)
            condition_rec(g_parent, parent)
        end
    end
    nothing
end
condition_rec(cond::Sample) = begin
    for arg in cond.args
        condition_rec(arg, cond)
    end
    nothing
end
condition_rec() = begin
    for cond in conditions
        condition_rec(cond)
    end
    nothing
end

apply_to_sdg(f::Function) = begin
    f(samples)
    f(dets)
    f(getindexes)
    f(conditions)
end
gc_church() = begin
    #Remove redundant dependencies.
    apply_to_sdg(remove_deps)
    gc_disable()
    weaken_getindexes()
    apply_to_sdg(weaken_deps)
    condition_rec()
    apply_to_sdg(weaken_store)
    gc_enable()
    gc()
    gc_disable()
    apply_to_sdg(strengthen_store)
    apply_to_sdg(strengthen_deps)
    strengthen_getindexes()
    gc_enable()
end
n_samples() = length(samples)
include("datastructures.jl")
end
