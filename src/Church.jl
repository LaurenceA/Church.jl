module Church
using Distributions

export sample, condition, value, resample, gc_church, background_sampler, n_samples, pdf, SDG, @lift

#Types
type Sample{V}
    det::Any
    dist::Distribution
    value::V
    logp::Float64
    deps::Vector
end
type Det{F <: Union(Function, DataType), T <: Tuple}
    f::F
    args::T
    deps::Vector
end
type GetIndex{S, T <: Tuple}
    struct::S
    args::T
    deps::Vector
    #Value only used in GC phase.
    value::Any
end
SDG = Union(Sample, Det, GetIndex)
SD = Union(Sample, Det)
WSDG = Union(WeakRef, Sample, Det, GetIndex)
type NoCond
end
const nocond = NoCond()

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
sample(det, cond::NoCond) = begin
    dist = value(det)
    val = rand(dist)
    lp = logp(dist, val)
    s = Sample(det, dist, val, lp, Array(WSDG, 0))
    push!(samples, s)
    add_dep(s, det)
    s
end
sample(det, val) = begin
    dist = value(det)
    lp = logp(dist, val)
    s = Sample(det, dist, val, lp, Array(WSDG, 0))
    add_dep(s, det)
    s
end
import Base.getindex
#Construct Det
Det(f::Union(Function, DataType), args) = begin
    d = Det(f, args, Array(WSDG, 0))
    push!(dets, d)
    for arg in args
        add_dep(d, arg)
    end
    d
end

#getindex(f::Function, args...) = Det(f, args)

#Construct GetIndex
#getindex(s, a1, a2, a3, a4::SDG, args...) = 
#    GI(s, a1, a2, a3, a4, args...)
#getindex(s, a1, a2, a3::SDG, args...) = 
#    GI(s, a1, a2, a3, args...)
#getindex(s, a1, a2::SDG, args...) = 
#    GI(s, a1, a2, args...)
#getindex(s, a1::SDG, args...) = 
#    GI(s, a1, args...)
#getindex(s, a1::SDG) = 
#    GI(s, a1)
GetIndex(struct, args...) = begin
    g = GetIndex(struct, args, Array(WSDG, 0), nothing)
    push!(getindexes, g)
    add_dep(g, struct[map(value, args)...])
    for arg in args
        add_dep(g, arg)
    end
    g
end

#Make this code work!
a(i::Int) = symbol("a$i")
type_a(i::Int, b::Bool) =
    b ? :($(a(i))::SDG) : a(i)
lift(f::Symbol, issdg::Vector{Bool}) =
    Expr(:(=), Expr(:call, f, {type_a(i,  issdg[i]) for i = 1:length(issdg)}...),
               Expr(:call, :Det, f, Expr(:tuple, map(a, 1:length(issdg))...)))
lift_gi(typ, issdg::Vector{Bool}) =
    Expr(:(=), Expr(:call, :getindex, :(s::$typ),{type_a(i,  issdg[i]) for i = 1:length(issdg)}...),
               Expr(:call, :GetIndex, :s, map(a, 1:length(issdg))...))
lift(f, n::Int, inner_func::Function=lift) = begin
    defs = Expr[]
    for nsamps in n:(-1):1
        for c in combinations(1:n, nsamps)
            issdg = fill(false, n)
            for i in c
                issdg[i] = true
            end
            push!(defs, inner_func(f, issdg))
        end
    end
    Expr(:block, defs...)
end
macro lift_gi(typ, n)
    esc(lift(typ, n, lift_gi))
end
macro lift(f, n)
    esc(lift(f, n))
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
        $fsym(args...; condition=nocond) = begin
            det = any(issdg, args) ? Det($tsym, args) : $tsym(args...)
            sample(det, condition)
        end
    end)
    eval(parse("export $(string(fsym))"))
end
issdg(x::SDG) = true
issdg(x) = false

#Operators construct a Det if called with a sample.
for op in [+, -, *, /, \, .*, ./, .\]
    op(a::SDG, b::SDG) = Det(op, (a, b))
    op(a::SDG, b) = Det(op, (a, b))
    op(a, b::SDG) = Det(op, (a, b))
end
include("gc.jl")
include("datastructures.jl")
end
