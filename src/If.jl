export If
type If
    f_true::Function
    f_false::Function
    val_true
    val_false
    If(f1::Function, f2::Function) = new(f1, f2)
end

import Base.getindex
getindex(i::If, cond::Int) = getindex(i, bool(cond))
getindex(i::If, cond::Bool) =
    if cond
        if !isdefined(i, :val_true)
            i.val_true = i.f_true()
        end
        i.val_true
    else
        if !isdefined(i, :val_false)
            i.val_false = i.f_false()
        end
        i.val_false
    end


