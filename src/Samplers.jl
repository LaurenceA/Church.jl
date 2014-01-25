abstract Sampler
abstract Prior <: Sampler

immutable Prior1 <: Prior
    s::WSample
end
immutable
resample(p::Prior1) = begin
    s = p.var
    old_val  = s.value
    old_logp = log_likelihood(s)
    proposal = rand(value(s.det))
    new_logp = log_likelihood(s, proposal)
    if rand() < exp(new_logp - old_logp)
        #Accept change
        s.value = proposal
    end
    nothing
end
