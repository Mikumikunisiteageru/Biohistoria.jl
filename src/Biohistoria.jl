# src/Biohistoria.jl

module Biohistoria

using Entropics

export lineage, LAR, integrate

struct Lineage
	distribution::Entropics.Distribution
	properties::Dict{Symbol,Any}
	Lineage(distribution; kwargs...) = new(distribution, kwargs)
end

function lineage(data; 
		old = convert(float(eltype(data)), Inf), now = zero(float(eltype(data))), kwargs...)
	distribution = bound(smooth(sample(data)), now, old)
	return Lineage(distribution; kwargs...)
end

struct LAR <: Function
	lineages::Vector{Lineage}
	properties::Dict{Symbol,Any}
	LAR(lineages; kwargs...) = new(lineages, kwargs)
end

function integrate(lar::LAR, near::Real=0.0, far::Real=Inf)
	near > far && ((near, far) = (far, near))
	scdf(x) = sum(map.(cdf.(getfield.(lar.lineages, :distribution)), (x,)))
	return scdf(far) - scdf(near)
end
integrate(lar::LAR, limits::NTuple{2,Real}) = integrate(lar, limits...)

(lar::LAR)(x) = sum(map.(pdf.(getfield.(lar.lineages, :distribution)), (x,)))

end # module Biohistoria
