# example/levatai/fig4.jl

using Dates          # 1.10.2
using DelimitedFiles # 1.10.2
using Distributions  # 0.25.109
using GeoInterface   # 1.3.5
using LibGEOS        # 0.9.2
using Printf         # 1.10.2
using PyPlot         # 2.11.2, Python 3.12.3, MatPlotLib 3.8.4

file(pathcomp...) = joinpath(ENV["LEVATAI_PATH"], pathcomp...)

function readint(filename::AbstractString)
	lines = readlines(file("Ages", "TaiwanEndemism_divergence_time", filename))
	S07, S7d, S0d, _, S0i = parse.(Float64, last.(split.(lines, '\t')))
	SvS = S0d / S0i
	H07 = S07 / (0.7 - 0.0)
	H7d = S7d / (1.2 - 0.7)
	HvH = H07 / H7d
	return S0i, S0d, SvS, H7d, H07, HvH
end
table = Matrix{Float64}(undef, 4, 6)
table[1, :] .= readint("3k-4k_altitude_int.txt")
table[2, :] .= readint("2k-3k_altitude_int.txt")
table[3, :] .= readint("1k-2k_altitude_int.txt")
table[4, :] .= readint("0-1k_altitude_int.txt")

matrix = readdlm(file("Section.tsv"))
mx = matrix[:, 1]
mx = (mx .- minimum(mx)) ./ (maximum(mx) - minimum(mx))
realize(t) = isa(t, Real) ? t : NaN
my = realize.(matrix[:, 2])
tw = Polygon([[collect.(zip(mx, my)); [[mx[1], my[1]]]]])

# original = ["#2a6544", "#fee5a1", "#ce7e41", "#984426"]
palette = ["#9bf1c1", "#fee5a1", "#f1bf98", "#dd9b84"]
heights = 0:1000:4000

function geom2mat(p::Polygon)
	coords = Vector{Float64}[]
	for ls = getgeom(p)
		for pt = getgeom(ls)
			push!(coords, collect(getcoord(pt)))
		end
		push!(coords, [NaN, NaN])
	end
	return coords[1:end-1, :]
end
function geom2mat(mp::MultiPolygon)
	coords = Vector{Float64}[]
	for p = getgeom(mp)
		for ls = getgeom(p)
			for pt = getgeom(ls)
				push!(coords, collect(getcoord(pt)))
			end
			push!(coords, [NaN, NaN])
		end
		push!(coords, [NaN, NaN])
	end
	return coords[1:end-1, :]
end

function cutfrom(l, u)
	band = Polygon([[[0.0,l], [0.0,u], [1.0,u], [1.0,l], [0.0,l]]])
	is = intersection(tw, band)
	mat = geom2mat(is)
	return first.(mat), last.(mat)
end

close()
figure(figsize=(6.4, 3.2))

try axt.remove() catch; end
axt = PyPlot.axes([0.11, 0.03, 0.88, 0.40])
for i = 1:4
	l, u = heights[i:i+1]
	c = palette[i]
	x, y = cutfrom(l, u)
	axt.fill([-1,-1,2,2], [l,u,u,l], c; lw=0, alpha=0.3)
	axt.fill(x, y, c; lw=0)
end
# axt.plot(mx, my, "k"; lw=0.2)
axt.set_xlim(-0.01, 1.005)
# axt.set_xticks(LinRange(0, 1, 13))
axt.set_xticks([])
axt.set_ylim(0, 4000)
axt.set_ylabel("Elevation (m)"; labelpad=3)
for j = 1:6
	dx = maximum(table[:, j]) > 10 ? 0.041 : 0.048
	for i = 1:4
		x = j/6 - dx
		y = i*1000 - 520
		t = @sprintf("%.3f", table[5-i, j])
		axt.text(x, y, t; ha="right", va="center")
	end
end
axt.text(1/12, 4650, raw"$S_{[\infty, 0]}$"; ha="center", va="center")
axt.text(3/12, 4650, raw"$S_{[1.2, 0]}$"; ha="center", va="center")
axt.text(5/12, 4650, raw"$\frac{S_{[1.2, 0]}}{S_{[\infty, 0]}}$"; ha="center", va="center", fontsize=14)
axt.text(7/12, 4650, raw"$H_{[1.2, 0.7]}$"; ha="center", va="center")
axt.text(9/12, 4650, raw"$H_{[0.7, 0]}$"; ha="center", va="center")
axt.text(11/12, 4650, raw"$\frac{H_{[0.7, 0]}}{H_{[1.2, 0.7]}}$"; ha="center", va="center", fontsize=14)

normal = Normal()
f(x) = pdf(normal, x)
x = -3:0.03:3
y = f.(x)
sx = -1.62:0.03:0.6
sy = f.(sx)
s = cdf(normal, 0.6) - cdf(normal, -1.62)
h = s / (0.6 + 1.62)

try ax1.remove() catch; end
ax1 = PyPlot.axes([0.02, 0.58, 0.3, 0.4])
ax1.plot(x, y, "-k"; lw=1.2)
ax1.plot([-3, 3], [0, 0], "-k"; lw=1.2)
ax1.fill_between(sx, 0, sy; color="#66ccff", lw=0, alpha=0.5)
ax1.fill_between(sx, 0, sy; color="#c584df", lw=0, alpha=0.5)
ax1.plot([-1.62, -1.62], [-0.02, f(-1.62)], "-k"; lw=0.8)
ax1.plot([0.6, 0.6], [-0.02, f(0.6)], "-k"; lw=0.8)
ax1.text(-1.62, -0.06, raw"$a$"; ha="center", va="center")
ax1.text(0.6, -0.06, raw"$b$"; ha="center", va="center")
ax1.set_xlim(-3, 3)
ax1.set_ylim(-0.07, 0.43)
ax1.text(-0.35, 0.119, raw"$S_{[a, b]}$"; ha="center", va="center")
ax1.text(1.63, 0.248, "LAR"; ha="center", va="center")
ax1.set_axis_off()

try ax2.remove() catch; end
ax2 = PyPlot.axes([0.34, 0.58, 0.3, 0.4])
ax2.plot(x, y, "-k"; lw=1.2)
ax2.plot([-3, 3], [0, 0], "-k"; lw=1.2)
ax2.fill_between(sx, 0, sy; color="#66ccff", lw=0, alpha=0.5)
ax2.fill_between([-1.62, 0.6], 0, [h, h]; color="#c584df", lw=0, alpha=0.5)
ax2.plot([-1.62, -1.62], [-0.02, f(-1.62)], "-k"; lw=0.8)
ax2.plot([0.6, 0.6], [-0.02, f(0.6)], "-k"; lw=0.8)
ax2.text(-1.62, -0.06, raw"$a$"; ha="center", va="center")
ax2.text(0.6, -0.06, raw"$b$"; ha="center", va="center")
ax2.set_xlim(-3, 3)
ax2.set_ylim(-0.07, 0.43)
ax2.text(-0.35, 0.119, raw"$S_{[a, b]}$"; ha="center", va="center")
ax2.plot([-1.62, -1.62, 0.6], [f(-1.62), h, h], "k"; lw=0.8)
ax2.plot([0.6, 3], [h, h], "k"; lw=0.8)
ax2.plot([2.6, 2.6, NaN, 2.6, 2.6], [0, 0.35, NaN, 0.65, 1] .* h, "k"; lw=0.8)
ax2.text(2.6, 0.5h, raw"$H_{[a, b]}$"; ha="center", va="center")
ax2.text(-0.96, 0.39, "LAR"; ha="center", va="center")
ax2.set_axis_off()

try ax0.remove() catch; end
ax0 = PyPlot.axes([0, 0, 1, 1])
ax0.set_facecolor("None")
ax0.set_axis_off()
ax0.set_xlim(0, 1)
ax0.set_ylim(0, 1)
y0 = 0.913
dy = 0.066
ax0.text(0.677, y0 - 0dy, raw"$a$: Starting time (Ma)"; ha="left", va="center")
ax0.text(0.677, y0 - 1dy, raw"$b$: Stopping time (Ma)"; ha="left", va="center")
ax0.text(0.677, y0 - 2dy, raw"$S_{[a, b]}$: Area under LAR curve"; ha="left", va="center")
ax0.text(0.677, y0 - 3dy, raw"$H_{[a, b]}$: Average of LAR (Ma$^{-1}$)"; ha="left", va="center")
ax0.text(0.820, 0.619, raw"$S_{[a, b]} = H_{[a, b]} \cdot |a-b|$"; ha="center", va="center")

date = Dates.format(today(), "YYYYmmdd")
savefig(file("Fig4_$date.pdf"))
