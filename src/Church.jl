module Church
using Distributions

export value, resample, gc_church, background_sampler, n_samples, pdf, SDG, @lift, Det, nocond

#Types
type NoCond
end
const nocond = NoCond()
type NoSampler
end
const nosampler = NoSampler()
type Sample{V, S<:Union(NoSampler, Function)}
    det::Any
    dist::Distribution
    value::V
    logp::Float64
    deps::Vector
    sampler::S
end
type Det{F, T <: Tuple}
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
fallback_sampler(d::Distribution, v) = rand(d)
import Distributions.sample
sample(det, cond::NoCond, sampler::Union(NoSampler, Function)) = begin
    dist = value(det)
    val = rand(dist)
    lp = logp(dist, val)
    s = Sample(det, dist, val, lp, Array(WSDG, 0), sampler)
    push!(samples, s)
    add_dep(s, det)
    s
end
sample(det, val, sampler::Union(NoSampler, Function)) = begin
    dist = value(det)
    lp = logp(dist, val)
    s = Sample(det, dist, val, lp, Array(WSDG, 0), sampler)
    push!(conditions, s)
    add_dep(s, det)
    s
end
Det(f, args) = begin
    d = Det(f, args, Array(WSDG, 0))
    push!(dets, d)
    add_dep(d, f)
    for arg in args
        add_dep(d, arg)
    end
    d
end

GetIndex(struct, args...) = begin
    g = GetIndex(struct, args, Array(WSDG, 0), nothing)
    push!(getindexes, g)
    add_dep(g, struct[map(value, args)...])
    for arg in args
        add_dep(g, arg)
    end
    g
end
#Allow indexing into multidimensional samples.
getindex(s::Sample, args...) = Det((s_, args_...) -> s_[args_...], tuple(s, args...))

#Lifting functions.
a(i::Int) = symbol("a$i")
type_a(i::Int, b::Bool) =
    b ? :($(a(i))::SDG) : a(i)
lift(f::Symbol, issdg::Vector{Bool}) =
    Expr(:(=), Expr(:call, f, {type_a(i,  issdg[i]) for i = 1:length(issdg)}...),
               Expr(:call, :(Church.Det), f, Expr(:tuple, map(a, 1:length(issdg))...)))
lift_gi(typ, issdg::Vector{Bool}) =
    Expr(:(=), Expr(:call, :getindex, :(s::$typ),{type_a(i,  issdg[i]) for i = 1:length(issdg)}...),
               Expr(:call, :GetIndex, :s, map(a, 1:length(issdg))...))
lift(f, n::Int, inner_func::Function) = begin
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
    esc(quote
        if Main == Base.function_module($f)
            $(lift(f, n, lift))
        else
            eval(parse(string("import ", Base.function_module($f), ".", $f)))
            $(lift(f, n, lift))
        end
    end)
end

add_dep(s::SDG, arg::SDG) =
    if !in(s, arg.deps)
        push!(arg.deps, s)
    end
add_dep(s::SDG, arg) = nothing

#Get the value of an expression.
value(s::Sample) = s.value
value(det::Det) = 
    value(det.f)(map(value, det.args)...)
value(g::GetIndex) = begin
    res = g.struct[map(value, g.args)...]
    #Add g as a dependent of struct [arg]
    add_dep(g, res)
    value(res)
end
value(s) = s

#Don't care whether cdf or pdf is needed.
logp(d::ContinuousDistribution, x) = logpdf(d, x)
logp(d::DiscreteDistribution, x) = logpmf(d, x)

#Find dependents
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

#Remove redundant dependents.
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
    (s.value, old_logp, new_logp) = propose(s.dist, old_val, s.sampler)
    old_logp += mapreduce(dep->dep.logp, +, deps)
    new_dists = map(s -> value(s.det), deps)
    new_logps = map((d, s) -> logp(d, s.value), new_dists, deps)
    new_logp += sum(new_logps)
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
propose(d::Distribution, _, s::NoSampler) =
    (rand(d), 0., 0.)
propose(d::Distribution, old_val, s::Function) = begin
    prop_dist_forward = s(d, old_val)
    prop_val = rand(prop_dist_forward)
    prop_dist_backward = s(d, prop_val)
    (prop_val, 
     logp(d, old_val )-logp(prop_dist_backward, old_val ), 
     logp(d, prop_val)-logp(prop_dist_forward , prop_val))
end
resample() = 
    if length(samples) != 0
        index = rand(1:length(samples))
        resample_inner(samples[index])
    end
resample(iter::Int) =
    for i = 1:iter
        resample()
    end

#Run sampling + garbage collection in background.
background_sampler(gc=true) = 
    @async while true
        tic()
        resample()
        sleep(1E-6)
        if gc && rand() < toq()
            gc_church()
        end
    end

#Define distributions.
for dist in {:MvNormal, :MvTDist, filter!(isleaftype, subtypes(Distribution))...}
    ssym = match(r"[a-zA-Z]+", string(dist)).match
    tsym = symbol(ssym)
    fsym = symbol(lowercase(ssym))
    eval(quote 
        $fsym(args...; condition=nocond, sampler=nosampler) = begin
            det = any(issdg, args) ? Det($tsym, args) : $tsym(args...)
            sample(det, condition, sampler)
        end
    end)
    eval(parse("export $(string(fsym))"))
end
issdg(x::SDG) = true
issdg(x) = false

#Overload operators.
for op in [+, -, *, /, \, .*, ./, .\]
    op(a::SDG, b::SDG) = Det(op, (a, b))
    op(a::SDG, b) = Det(op, (a, b))
    op(a, b::SDG) = Det(op, (a, b))
    op(a::SDG) = Det(op, (a,))
end
include("gc.jl")
include("datastructures.jl")
end
