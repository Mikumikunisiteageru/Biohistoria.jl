# example/levatai/resampler.jl

module Resampler

export resample

gaussian(t) = sqrt(2*pi) \ exp(-2 \ t^2)

getstep(v::StepRange) = step(v)
getstep(v) = (maximum(v) - minimum(v)) / (length(v) - 1)

function rspgaussian(xf, yf, xc; sigma=getstep(xc)/2)
	function y(x)
		w = gaussian.((x .- xf) ./ sigma)
		return sum(w .* yf) / sum(w)
	end
	return y.(xc)
end

function rspuniform(xf, yf, xc; radius=getstep(xc)/2)
	function y(x)
		diff = abs.(x .- xf)
		idlt = diff .< radius
		ideq = diff .== radius
		return ((2*sum(yf[idlt]) + sum(yf[ideq])) / (2*count(idlt) + count(ideq)))
	end
	return y.(xc)
end

function resample(xf, yf, xc; h=getstep(xc)/2, kernel=:gaussian)
	if kernel == :gaussian
		return rspgaussian(xf, yf, xc; sigma=h)
	elseif kernel == :uniform
		return rspuniform(xf, yf, xc; radius=h)
	else
		throw(ArgumentError("kernel not supported"))
	end
end

end # module Resampler
