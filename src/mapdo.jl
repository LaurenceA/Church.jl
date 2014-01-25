import Base.Callable
# 0 argument function
mapdo(f::Callable) = f()
# 1 argument function
mapdo(f::Callable, t::()) = nothing
mapdo(f::Callable, t::(Any,)) = (f(t[1]); nothing)
mapdo(f::Callable, t::(Any, Any)) = (f(t[1]); f(t[2]); nothing)
mapdo(f::Callable, t::(Any, Any, Any)) = (f(t[1]); f(t[2]); f(t[3]); nothing)
mapdo(f::Callable, t::(Any, Any, Any, Any)) = (f(t[1]); f(t[2]); f(t[3]); f(t[4]); nothing)
mapdo(f::Callable, t::Tuple) = begin
    for i = 1:length(t)
        f(t[i])
    end
    nothing
end
# 2 argument function
mapdo(f::Callable, t::(), s::()) = ()
mapdo(f::Callable, t::(Any,), s::(Any,)) = (f(t[1],s[1]); nothing)
mapdo(f::Callable, t::(Any,Any), s::(Any,Any)) = (f(t[1],s[1]); f(t[2],s[2]); nothing)
mapdo(f::Callable, t::(Any,Any,Any), s::(Any,Any,Any)) =
    (f(t[1],s[1]); f(t[2],s[2]); f(t[3],s[3]); nothing)
mapdo(f::Callable, t::(Any,Any,Any,Any), s::(Any,Any,Any,Any)) =
    (f(t[1],s[1]); f(t[2],s[2]); f(t[3],s[3]); f(t[4],s[4]); nothing)
        # n argument function
mapdo(f::Callable, ts::Tuple...) = begin
    for i = 1:length(ts[1])
        f(map(t -> t[i], ts)...)
    end
    nothing
end
