#Must define isweak, weaken, strengthen and @lift_gi the datastructure.
import Base.getindex

#Don't weaken arrays or dicts
isweak(a) = false
weaken(a) = nothing
strengthen(a) = nothing
@lift_gi(Vector, 1)
@lift_gi(Matrix, 2)
@lift_gi(Dict, 1)

#If
export If, @If
type If
    f_true::Function
    f_false::Function
    val_true
    val_false
    If(f1::Function, f2::Function) = new(f1, f2, uneval, uneval)
end

type Uneval
end
const uneval = Uneval()

getindex(i::If, cond::Int) = getindex(i, bool(cond))
getindex(i::If, cond::Bool) =
    if cond
        if isa(i.val_true, Uneval)
            i.val_true = i.f_true()
        end
        i.val_true
    else
        if isa(i.val_false, Uneval)
            i.val_false = i.f_false()
        end
        i.val_false
    end

macro If(cond, val_true, val_false)
    esc(:(If(() -> $val_true, () -> $val_false)[$cond]))
end

isweak(i::If) = isa(i.val_true, WeakRef)
weaken(i::If) = begin
    i.val_true = WeakRef(i.val_true)
    i.val_false = WeakRef(i.val_false)
end
strengthen(wr::WeakRef) = (wr == WeakRef()) ? uneval : wr.value
strengthen(i::If) = begin
    i.val_true = strengthen(i.val_true)
    i.val_false = strengthen(i.val_false)
end
@lift_gi(If, 1)

export Mem#, InfiniteVector

#type InfiniteVector{T}
#    vec::Vector{T}
#    defined::Vector{Bool}
#end
#InfiniteVector{T}() = InfiniteVector(Array(T, 0), Array(Bool, 0))
#InfiniteVector() = InfiniteVector(Array(Any, 0), Array(Bool, 0))
#getindex(iv::InfiniteVector, i::Int) = begin
#    if length(iv.vec) < i || !iv.defined[i]
#        error("Index not defined in InfiniteVector")
#    end
#    iv.vec[i]
#end
#setindex!(iv::InfiniteVector, val, i::Int) = begin
#    while length(iv.vec) < i
#        old_vec = iv.vec
#        iv.vec = Array(eltype(iv.vec), length(old_vec))
#        append!(iv.vec, Array(eltype(iv.vec), i - length(iv.vec)))
#        append!(iv.defined, fill(false, i - length(iv.defined)))
#    end
#    iv.defined[i] = true
#    iv.vec[i] = val
#end
#import Base.haskey
#haskey(iv::InfiniteVector, i::Int) = begin
#    @assert length(iv.vec) == length(iv.defined)
#    length(iv.vec) >= i && iv.defined[i]
#end

type Mem
    f::Function
    dict::Dict{Any, WSDG}
end
Mem(f::Function) = Mem(f, Dict{Any, WSDG}())
getindex(m::Mem, args...) = begin
    if !haskey(m.dict, args)
        m.dict[args] = m.f(args...)
    end
    m.dict[args]
end
#Need to do something here!
#isweak(m::Mem) = isa(first(m.dict)[2], WeakRef)
isweak(m::Mem) = begin
    isa(first(m.dict)[2], WeakRef)
end
weaken(m::Mem) = begin
    for key in keys(m.dict)
        m.dict[key] = WeakRef(m.dict[key])
    end
    nothing
end
strengthen(m::Mem) = begin
    for key in collect(keys(m.dict))
        if m.dict[key] == WeakRef()
            delete!(m.dict, key)
        else
            m.dict[key] = m.dict[key].value
        end
    end
    nothing
end
@lift_gi(Mem, 1)
