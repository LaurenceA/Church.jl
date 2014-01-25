abstract Sampler

type SamplerVars
    sampler::Sampler
    vars::Union((Sample...), (WeakRef...))
end
SamplerVars(s::Sampler, vars::Sample...) = SamplerVars(s, vars)
resample(sv::SamplerVars) = resample(sv.sampler, sv.vars)

type NoSampler
end
const nosampler = NoSampler()

immutable Prior <: Sampler
end
const prior = Prior()


resample(::Prior, ss::(Sample...)) = begin
    old_logp = log_likelihood(ss)
    proposals = map(s -> rand(value(s.det)), ss)
    new_logp = log_likelihood(ss, proposals)
    if rand() < exp(new_logp - old_logp)
        #Accept change
        mapdo((s, p) -> s.value=p, ss, proposals)
    end
    nothing
end

immutable MHMC <: Sampler
    proposal::Function
end
resample(sampler::MHMC, ss::(Sample...)) = begin
    old_logp = log_joint(ss)
    old_vals = map(value, ss)
    forward_proposal_dists = tuplify(sampler.proposal(old_vals...))
    proposals = map(rand, forward_proposal_dists)
    forward_logp = sum(map(logp, forward_proposal_dists, proposals))
    backward_logp = sum(map(logp, tuplify(sampler.proposal(proposals...)), old_vals))
    new_logp = log_joint(ss, proposals)
    if rand() < exp(new_logp - old_logp + backward_logp - forward_logp)
        #Accept change
        mapdo((s, p) -> s.value=p, ss, proposals)
    end
    nothing
end

tuplify(x::(Any...)) = x
tuplify(x) = (x,)
