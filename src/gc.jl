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
        g.value = g.struct[map(value, g.args)...]
    end
    for g in getindexes
        if !isweak(g.struct)
            weaken(g.struct)
        end
    end
    nothing
end
strengthen_getindexes() = begin
    for g in getindexes
        if isweak(g.struct)
            strengthen(g.struct)
        end
    end
    for g in getindexes
        g.value = nothing
    end
    nothing
end

weaken(svv::Vector{SamplerVars}) = begin
    for sv in svv
        sv.var = WeakRef(sv.var)
    end
    nothing
end
strengthen(svv::Vector{SamplerVars}) = begin
    filter!(x -> x.var != WeakRef(), svv)
    for sv in svv
        @assert isa(sv.var, WeakRef)
        sv.var = sv.var.value
    end
    nothing
end

wr_value(wr::WeakRef) = wr.value
wr_value(sr) = sr

isdependent(parent::RV, child::Sample)   = parent == child.det
isdependent(parent::RV, child::Det)      = in(parent, child.args)
isdependent(parent::RV, child::GetIndex) = 
    in(parent, child.args) || (parent == child.struct[map(value, child.args)...])

args(s::Sample) = {s.det}
args(d::Det) = {d.args...}
args(g::GetIndex) = 
    {g.struct[map(value, g.args)...], g.args...}

condition_rec(parent, child::RV) = nothing
condition_rec(parent::RV, child::RV) = begin
    @assert isdependent(parent, child)
    i = findfirst(x -> x==wr_value(child), parent.deps)
    if isa(parent.deps[i], WeakRef)
        parent.deps[i] = parent.deps[i].value
        for g_parent in args(parent)
            condition_rec(g_parent, parent)
        end
    end
    nothing
end
condition_rec(cond::Sample) = begin
    condition_rec(cond.det, cond)
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
    #Weaken links in lazy structures.
    weaken_getindexes()
    #Weaken all dependent links.
    apply_to_sdg(weaken_deps)
    #Strengthen conditioned links.
    condition_rec()
    #Weaken links in samples, dets, etc.
    apply_to_sdg(weaken_store)
    #Weaken links from samplers
    weaken(samplers)
    gc_enable()
    gc()
    gc_disable()
    #Strengthen links in samples, dets, etc.
    strengthen(samplers)
    apply_to_sdg(strengthen_store)
    apply_to_sdg(strengthen_deps)
    strengthen_getindexes()
    gc_enable()
end
n_samples() = length(samples)
