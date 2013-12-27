module Church
using Distributions

export sample, condition, @branch, value, resample, gc_church, background_sampler, Branch, n_samples

#Types
type Sample{V}
    det::Any
    dist::Distribution
    value::V
    logp::Float64
    deps::Vector{Union(WeakRef, Sample)}
    conditioned::Bool
end
type Branch
    cond
    f_true::Function
    f_false::Function
    deps::Vector{Union(WeakRef, Sample)}
    conditioned::Bool
    val_true
    val_false
end
type Det{F <: Union(Function, DataType), T <: Tuple}
    f::F
    args::T
end
type Uneval
end
const uneval = Uneval()
SDB = Union(Sample, Det, Branch)

#Store declared samples/branches
const samples = Array(Union(Sample, WeakRef), 0)
const branches = Array(Union(Branch, WeakRef), 0)

#Show functions
import Base.show
show(io::IO, s::Sample) = print(io, "Sample", tuple(s.dist, s.value))
show(io::IO, b::Branch) = print(io, "Branch")

#Constructors
import Distributions.sample
sample(det) = begin
    dist = value(det)
    val = rand(dist)
    lp = logp(dist, val)
    s = Sample(det, dist, val, lp, Array(Union(WeakRef, Sample), 0), false)
    push!(samples, s)
    add_dep(s, det)
    s
end
condition(det, val) = begin
    dist = value(det)
    lp = logp(dist, val)
    s = Sample(det, dist, val, lp, Array(Union(WeakRef, Sample), 0), true)
    #push!(samples, s)
    add_dep(s, det)
    rec_cond(det)
    s
end
Branch(cond, f_true::Function, f_false::Function) = begin
    b = Branch(cond, f_true, f_false, Array(Union(WeakRef, Sample), 0), false, uneval, uneval)
    push!(branches, b)
    b
end

#Set conditioned bit.  Required by condition.
rec_cond(s::Sample) = begin
    s.conditioned=true
    rec_cond(s.det)
end
rec_cond(det::Det) =
    for arg in det.args
        rec_cond(arg)
    end
rec_cond(b::Branch) = begin
    b.conditioned = true
    rec_cond(b.cond)
    rec_cond(b.val_true)
    rec_cond(b.val_false)
end
rec_cond(any) = nothing

#Make the sample that you're currently creating a dependent of previous samples.
#Required by sample and condition.
add_dep(s::Sample, arg::Sample) = push!(arg.deps, s)
add_dep(s::Sample, det::Det) = 
    for arg in det.args
        add_dep(s, arg)
    end
add_dep(s::Sample, b::Branch) = begin
    #Add s as a dependent of the currently existing nodes
    add_dep(s, b.cond)
    add_dep(s, b.val_true)
    add_dep(s, b.val_false)
    #Add s as a dependent of the Branch.
    push!(b.deps, s)
end
add_dep(s::Sample, arg) = nothing

#Get the value of an expression.
value(det::Det) = begin
    det.f(map(value, det.args)...)
end
value(s::Sample) = s.value
value(s) = s
value(b::Branch) =
    if bool(value(b.cond))
        if isa(b.val_true, Uneval)
            b.val_true = b.f_true()
            #Add dependents to new nodes
            map(s -> add_dep(s, b.val_true), b.deps)
            #Add conditioning if conditioned
            if b.conditioned
                rec_cond(b.val_true)
            end
        end
        value(b.val_true)
    else
        if isa(b.val_false, Uneval)
            b.val_false = b.f_false()
            #Add dependents to new nodes
            map(s -> add_dep(s, b.val_false), b.deps)
            #Adding conditioning if conditioned
            if b.conditioned
                rec_cond(b.val_false)
            end
        end
        value(b.val_false)
    end

#Don't care whether cdf or pdf is needed.
logp(d::ContinuousDistribution, x) = logpdf(d, x)
logp(d::DiscreteDistribution, x) = logpmf(d, x)
#old_logp(s::Sample) = s.logp

#Propose, accept/reject.
resample_inner(s::Sample) = begin
    old_val  = s.value
    deps = s.deps
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
        if rand() < toq()
            gc_church()
        end
    end

#Sugar
#Allow construction of Det by f[sample()] syntax.  Breaks Float64[]
import Base.getindex
getindex(f::Function, args...) = Det(f, args)
getindex(dist::DataType, args...) = Det(dist, args)

#Use a macro to construct a branch
macro branch(cond, val_true, val_false)
    esc(:(Branch($cond, () -> $val_true, () -> $val_false)))
end

#Operators construct a Det if called with a sample.
for op in [+, -, *, /, \, .*, ./, .\]
    op(a::SDB, b::SDB) = Det(op, (a, b))
    op(a::SDB, b) = Det(op, (a, b))
    op(a, b::SDB) = Det(op, (a, b))
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
weaken_deps(as::Vector) = begin
    for s in as
        for i = 1:length(s.deps)
            if !s.deps[i].conditioned
                s.deps[i] = WeakRef(s.deps[i])
            end
        end
    end
    nothing
end
strengthen_deps(as::Vector) = begin
    for s in as
        filter!(x -> x != WeakRef(), s.deps)
        for i = 1:length(s.deps)
            if isa(s.deps[i], WeakRef)
                s.deps[i] = s.deps[i].value
            end
        end
    end
    nothing
end
weaken_branches() = begin
    for b in branches
        if bool(value(b.cond))
            b.val_false = WeakRef(b.val_false)
        else
            b.val_true = WeakRef(b.val_true)
        end
    end
    nothing
end
strengthen_branches() = begin
    for branch in branches
        if bool(value(branch.cond))
            branch.val_false = strengthen_if(branch.val_false)
        else
            branch.val_true = strengthen_if(branch.val_true)
        end
    end
    nothing
end
gc_church() = begin
    gc_disable()
    weaken_branches()
    weaken_deps(samples); weaken_deps(branches)
    weaken(samples); weaken(branches)
    gc_enable()
    gc();gc();gc()
    gc_disable()
    strengthen(samples); strengthen(branches)
    strengthen_deps(samples); strengthen_deps(branches)
    strengthen_branches()
    gc_enable()
end
strengthen_if(wr::WeakRef) =
    wr == WeakRef() ? uneval : wr.value
n_samples() = length(samples)
end
