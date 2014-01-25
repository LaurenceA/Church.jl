module Church
using Distributions

export value, resample, gc_church, background_sampler, n_samples, pdf, RV, @lift, Det, nocond, Prior, MHMC

#Types
type NoCond
end
const nocond = NoCond()
abstract RV
type Sample{V} <: RV
    det::Any
    value::V
    deps::Vector
end
type Det{F, T <: Tuple} <: RV
    f::F
    args::T
    deps::Vector
end
type GetIndex{S, T <: Tuple} <: RV
    struct::S
    args::T
    deps::Vector
    #Value only used in GC phase.
    value::Any
end
WRV = Union(WeakRef, RV)
WSample = Union(WeakRef, Sample)

#Store declared samples/branches
include("Sampler.jl")
const samplers = Array(SamplerVars, 0)
const samples = Array(Union(Sample, WeakRef), 0)
const dets = Array(Union(Det, WeakRef), 0)
const getindexes = Array(Union(GetIndex, WeakRef), 0)
const conditions = Array(Union(Sample, WeakRef), 0)

#Show functions
import Base.show
show(io::IO, s::RV) = print(io, "RV(", value(s), ")")

#Constructors
import Distributions.sample
sample(det, cond::NoCond, sampler::Sampler) = begin
    s = sample(det, cond)
    push!(samplers, SamplerVars(sampler, s))
    s
end
sample(det, cond::NoCond, sampler::NoSampler) = 
    s = sample(det, cond)
sample(det, cond::NoCond) = begin
    dist = value(det)
    val = rand(dist)
    s = Sample(det, val, Array(WRV, 0))
    push!(samples, s)
    add_dep(s, det)
    s
end
sample(det, val, sampler::Union(Sampler, NoSampler)) = begin
    s = Sample(det, val, Array(WRV, 0))
    push!(conditions, s)
    add_dep(s, det)
    s
end
Det(f, args) = begin
    d = Det(f, args, Array(WRV, 0))
    push!(dets, d)
    add_dep(d, f)
    for arg in args
        add_dep(d, arg)
    end
    d
end

GetIndex(struct, args...) = begin
    g = GetIndex(struct, args, Array(WRV, 0), nothing)
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
    b ? :($(a(i))::RV) : a(i)
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

add_dep(s::RV, arg::RV) =
    if !in(s, arg.deps)
        push!(arg.deps, s)
    end
add_dep(s::RV, arg) = nothing

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
#Compute logp for samples
logp(s::Sample) = logp(value(s.det), s.value)
log_likelihood(s) = mapreduce(logp, +, deps_likelihood(s))
log_joint(s) = mapreduce(logp, +, deps_joint(s))
#log_likelihood(s::(Sample,), val) = begin
#    old_val = s[1].value
#    s[1].value = val[1]
#    result = log_likelihood(s)
#    s[1].value = old_val
#    result
#end
#log_joint(s::(Sample,), val) = begin
#    old_val = s.value
#    s.value = val
#    result = log_joint(s)
#    s.value = old_val
#    result
#end
log_likelihood(ss::(Sample...), vals::(Any...)) = begin
    old_vals = map(s -> s.value, ss)
    #for i = 1:length(ss)
    #    ss[i].value = vals[i]
    #end
    map((s, v) -> (s.value=v; nothing), ss, vals)
    result = log_likelihood(ss)
    #for i = 1:length(ss)
    #    ss[i].value = old_vals[i]
    #end
    map((s, v) -> (s.value=v; nothing), ss, old_vals)
    result
end
#log_joint(ss::Vector{WSample}, vals::Vector) = begin
#    old_vals = map(s -> s.value, ss)
#    for i = 1:length(ss)
#        ss[i].value = vals[i]
#    end
#    result = log_joint(ss)
#    for i = 1:length(ss)
#        ss[i].value = old_vals[i]
#    end
#    result
#end

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
deps(s::Sample) = begin
    result = Sample[]
    deps_recurse(s, result)
    result
end
#deps_joint(s::(Sample,)) = vcat(s[1], deps(s[1]))
#deps_joint(ss::(Sample...)) =
#    union(ss, map(s -> deps(s), ss)...)
deps_likelihood(s::(Sample,)) = deps(s[1])
#deps_likelihood(ss::(Sample...)) =
#    filter(s -> !in(s, ss), deps_joint(ss))


#Remove redundant dependents.
remove_deps(origin::RV) =
    filter!(x -> isdependent(origin, x), origin.deps)
remove_deps(as::Vector) = begin
    for a in as
        remove_deps(a)
    end
    nothing
end

resample() = 
    if length(samplers) != 0
        index = rand(1:length(samplers))
        resample(samplers[index])
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
        $fsym(args...; condition=nocond, sampler=prior) = begin
            det = any(issdg, args) ? Det($tsym, args) : $tsym(args...)
            sample(det, condition, sampler)
        end
    end)
    eval(parse("export $(string(fsym))"))
end
issdg(x::RV) = true
issdg(x) = false

#Overload operators.
for op in [+, -, *, /, \, .*, ./, .\]
    op(a::RV, b::RV) = Det(op, (a, b))
    op(a::RV, b) = Det(op, (a, b))
    op(a, b::RV) = Det(op, (a, b))
    op(a::RV) = Det(op, (a,))
end
include("gc.jl")
include("datastructures.jl")
end
