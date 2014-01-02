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

import Base.getindex
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

weaken(i::If) = begin
    i.val_true = WeakRef(i.val_true)
    i.val_false = WeakRef(i.val_false)
end
strengthen(wr::WeakRef) = (wr == WeakRef()) ? uneval : wr.value
strengthen(i::If) = begin
    i.val_true = strengthen(i.val_true)
    i.val_false = strengthen(i.val_false)
end
