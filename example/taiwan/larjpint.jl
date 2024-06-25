# example/taiwan/larjpint.jl

using Biohistoria
using DelimitedFiles
using JointPoint
using ZipFile

const x = 0:0.01:10
const kk = 0:3
const intlimits = [(0.0, 0.7), (0.7, 1.2), (0.0, 1.2), (1.2, Inf), (0, Inf)]

function calclar(dir, lar)
	y = lar(x)
	writedlm(dir * "_lar.txt", hcat(x, y))
	return y
end

function calcjp(dir, y)
	jpresults = []
	for k = kk
		xj, yj = findjoint(x, y, k)
		for (xji, yji) = zip(xj, yj)
			push!(jpresults, (k, xji, yji))
		end
	end
	writedlm(dir * "_jp.txt", permutedims(hcat(collect.(jpresults)...)))
end

function calcint(dir, lar)
	intresults = integrate.((lar,), intlimits)
	writedlm(dir * "_int.txt", hcat(intlimits, intresults))
end

function work(dir::AbstractString)
	@info(dir)
	files = readdir(dir; join=true)
	datas = vec.(readdlm.(files))
	lineages = lineage.(datas)
	lar = LAR(lineages)
	@time y = calclar(dir, lar)
	@time calcjp(dir, y)
	@time calcint(dir, lar)
end

function workdirs(root::AbstractString)
	w = ZipFile.Writer(root * "_larjpint.zip")
	for (dir, subdirs, _) = walkdir(root)
		isempty(subdirs) || continue
		work(dir)
		reldir = strip(chop(dir; head=length(root), tail=0), ['/', '\\'])
		for suffix = ["_lar.txt", "_jp.txt", "_int.txt"]
			t = ZipFile.addfile(w, reldir * suffix)
			write(t, read(dir * suffix))
			close(t)
		end
	end
	close(w)
end

# workdirs(root)
