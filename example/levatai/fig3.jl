# example/levatai/fig3.jl

using Colors          # 0.12.11
using Dates           # 1.10.2
using DelimitedFiles  # 1.10.2
using GeologicTime    # 0.1.0
using HypothesisTests # 0.11.0
using Printf          # 1.10.2
using PyPlot          # 2.11.2, Python 3.12.3, MatPlotLib 3.8.4
using Statistics      # 1.10.2

file(pathcomp...) = joinpath(ENV["LEVATAI_PATH"], pathcomp...)

carex = "#f6ecce" # "#f0d78d"

function readxy(path::AbstractString; maxx=nextfloat(2.6f0))
	x, y = eachcol(readdlm(path))
	id = x .<= maxx
	return x[id], y[id]
end

time, sea = readxy(file("Factors", "SeaElev_2M6.tsv"))
time, mag = readxy(file("Factors", "Magnet_2M6.tsv"))
time, temp = readxy(file("Factors", "Temp_2M6.tsv"))
time, prec = readxy(file("Factors", "Prec_2M6.tsv"))
time, dp = readxy(file("Ages", "TaiwanEndemism_dispersal_node_time_lar.txt"))
time, edm = readxy(file("Ages", "TaiwanEndemism_divergence_time", "All_endemism_lar.txt"))
time, edm34 = readxy(file("Ages", "TaiwanEndemism_divergence_time", "3k-4k_altitude_lar.txt"))
time, edm23 = readxy(file("Ages", "TaiwanEndemism_divergence_time", "2k-3k_altitude_lar.txt"))
time, edm12 = readxy(file("Ages", "TaiwanEndemism_divergence_time", "1k-2k_altitude_lar.txt"))
time, edm01 = readxy(file("Ages", "TaiwanEndemism_divergence_time", "0-1k_altitude_lar.txt"))
qq = [sea, mag, temp, prec, dp, edm, edm34, edm23, edm12, edm01]

k = 8 / 6.4
close()
figure(figsize=(6.4, 6.4k))

function drawroute(ax, x, y; c=carex, sx=0.02, hx=0.03, rx=0.01, sy=0.02/k, hy=hx/k, ry=rx/k, kwargs...)
	tt = LinRange(0, pi/2, 51)
	xx = [sx; sx; 
		x-hx-rx .+ (rx+2hx) .* cos.(reverse(tt)); 
		x+hx; x-hx; 
		x-hx-rx .+ rx .* cos.(tt); 
		sx; sx]
	yy = [y; y+hy; 
		y-hy-ry .+ (ry+2hy) .* sin.(reverse(tt)); 
		sy; sy; 
		y-hy-ry .+ ry .* sin.(tt); 
		y-hy; y]
	ax.fill(xx, yy; c=c, lw=0, kwargs...)
end
function drawroute(ax, y; c=carex, sx=0.02, tx=1-sx, hy=0.03/k, kwargs...)
	xx = [sx, sx, tx, tx]
	yy = [y-hy, y+hy, y+hy, y-hy]
	ax.fill(xx, yy; c=c, lw=0, kwargs...)
end
dg = 0.11
x = 0.54 .+ (0:3) .* dg
y = [0.09+0dg/k, 0.09+dg/k, 0.09+2dg/k, 0.09+3dg/k, 
	 0.09+4dg/k, 2\(1+dg/k), 
	 0.91-3dg/k, 0.91-2dg/k, 0.91-1dg/k, 0.91-0dg/k] |> reverse
xa, xb = 0.25, 0.47
dp = 0.05 + dg/k - 0.13

try ax0.remove() catch; end
ax0 = PyPlot.axes([0, 0, 1, 1])
ax0.set_facecolor("None")
ax0.set_xlim(0, 1)
ax0.set_ylim(0, 1)
ax0.set_axis_off()

try
	for i = 1:10
		axx[i].remove()
	end
catch
	global axx = Vector(undef, 10)
end
for i = 1:10
	if i <= 4
		drawroute(ax0, x[i], y[i])
	else
		drawroute(ax0, y[i])
	end
	# ax0.fill([xa, xb, xb, xa, xa], y[i] .+ [-0.04, -0.04, +0.04, +0.04, -0.04]; c=carex, ec="k", lw=0.0)
	axx[i] = PyPlot.axes([xa, y[i]-0.04, xb-xa, 0.08])
	axx[i].plot(time, qq[i], "k"; lw = i <= 4 ? 0.5 : 0.8)
	axx[i].set_xlim(2.6, 0)
	axx[i].set_facecolor(carex)
	# axx[i].set_axis_off()
	axx[i].set_xticks([])
	axx[i].set_yticks([])
end

try axr.remove() catch; end
axr = PyPlot.axes([xa, 0.02, xb-xa, 0.03-dp])
axr.fill([0,0,1,1], [0,1,1,0]; c="k", lw=0)
axr.fill([2,2,3,3], [0,1,1,0]; c="k", lw=0)
axr.set_xlim([2.6, 0])
axr.set_ylim([0, 1])
axr.set_xticks([])
axr.set_yticks([])

try axg.remove() catch; end
axg = PyPlot.axes([xa, 0.95+dp, xb-xa, 0.03-dp])
drawtimescale(axg, 2.6, 0, [4]; facealpha=1, texts = Dict("Pleistocene" => "Pleistocene"))

function maketext(ax, i, line1)
	ax.text(2\(0.02+xa), y[i], line1; ha="center", va="center")
end
function maketext(ax, i, line1, line2)
	ax.text(0.033, y[i] + 0.01, line1; ha="left", va="center")
	ax.text(xa-0.013, y[i] - 0.01, line2; ha="right", va="center")
end
maketext(ax0, 1, "Sea level")
maketext(ax0, 2, "Soil magnetic", "susceptibility")
maketext(ax0, 3, "Temperature")
maketext(ax0, 4, "Precipitation")
maketext(ax0, 5, "LAR, dispersed")
maketext(ax0, 6, "LAR, endemic", "all-elevation")
maketext(ax0, 7, "LAR, endemic", "3000–4000 m") # en dash "–"
maketext(ax0, 8, "LAR, endemic", "2000–3000 m")
maketext(ax0, 9, "LAR, endemic", "1000–2000 m")
maketext(ax0, 10, "LAR, endemic", "0–1000 m")
ax0.text(2\(0.02+xa), 0.965+2\dp, "Geologic Time"; ha="center", va="center")
ax0.text(2\(0.02+xa), 0.035-2\dp, "Time (Ma)"; ha="center", va="center")
for (i, t) = zip(1:4, ["SL", "SMS", "T", "P"])
	ax0.text(x[i], 0.035-2\dp, t; ha="center", va="center")
	ax0.text(x[i] .- 0.01, y[i] .- 0.01/k, t; ha="center", va="center")
end
for (i, t) = zip(5:10, ["D", "E", "E\$_{\\,34}\$", "E\$_{\\,23}\$", "E\$_{\\,12}\$", "E\$_{\\,01}\$"])
	ax0.text(0.95, y[i], t; ha="center", va="center")
end

cm = ColorMap("Spectral_r")
# cm = ColorMap("RdBu_r")
function tcwt(c)
	ctest = RGBA(c...)
	black = colorant"black"
	white = colorant"white"
	if colordiff(ctest, black) > colordiff(ctest, white)
		return "k", "normal"
	else
		return "w", "bold" # not working for formula, whatever
	end
end

corrs = Matrix(undef, 10, 4)
for i = 2:10
	for j = 1:4
		j >= i && continue
		r = cor(qq[i], qq[j])
		p = pvalue(CorrelationTest(qq[i], qq[j]))
		corrs[i, j] = (r, p)
		@assert p < 0.001
		c = cm(2 \ (1 + r))
		tc, wt = tcwt(c)
		ax0.fill(x[j] .+ [-0.045, -0.045, 0.045, 0.045], y[i] .+ [-0.045, 0.045, 0.045, -0.045] ./ k; lw=0.8, ec="k", c=c)
		ax0.text(x[j], y[i], @sprintf("\$%+.3f\$", r); c=tc, ha="center", va="center", weight=wt)
	end
end

xl1 = 0.6
xl3 = 0.975
yl3 = 0.98
yl2 = 0.95 + dp
xl2 = xl3 - (yl3-yl2)*k
yl1 = yl3 - (xl3-xl1)/k
ax0.plot([xl2, xl3, xl3, xl2, xl2, xl1, xl1, xl2], [yl3, yl3, yl1, yl1, yl2, yl2, yl3, yl3], "k"; lw=0.8)
yl4 = 0.951
xl4 = xl2 - (yl2-yl4)*k
xl5 = 0.9
yl5 = yl2 - (xl2-xl5)/k
for r = 0:0.2:1
	xl = r * xl1 + (1-r) * xl5
	yl = r * yl1 + (1-r) * yl5
	ax0.plot([xl4, xl2], [yl, yl], "k"; lw=0.8)
	ax0.plot([xl, xl], [yl4, yl2], "k"; lw=0.8)
	ax0.text(xl, 0.935, @sprintf("\$%+.1f\$", r); ha="center", va="center")
	ax0.text(0.94-(0.98-xl3), yl, @sprintf("\$%-.1f\$", -r); ha="right", va="center")
end

n = 51
larm = Array{Float64}(undef, 1, n, 3)
barm = Array{Float64}(undef, n, 1, 3)
for (i, r) = enumerate(LinRange(0, 1, n))
	larm[1, i, :] .= cm(2\(1+r))[1:3]
	barm[i, 1, :] .= cm(2\(1-r))[1:3]
end
ax0.imshow(larm; extent = [xl5, xl1, yl2, yl3], aspect="auto", interpolation="none", zorder=1)
ax0.imshow(barm; extent = [xl2, xl3, yl1, yl5], aspect="auto", interpolation="none", zorder=1)
ax0.fill([xl5-0.05, xl5-0.05, xl3, xl3, xl2, xl2], [yl2, yl3, yl3, yl5-0.05, yl5-0.05, yl2]; c=cm(0.5), lw=0, zorder=-2)

ax0.text(0.797, 0.895, "Pearson's"; ha="center", va="center")
ax0.text(0.797, 0.875, "correlation"; ha="center", va="center")
ax0.text(0.797, 0.855, "coefficient \$r\$"; ha="center", va="center")

ax0.text(0.797, 0.797, "\$p<0.001\$"; ha="center", va="center")
ax0.text(0.797, 0.777, "for all"; ha="center", va="center")

date = Dates.format(today(), "YYYYmmdd")
savefig(file("Fig3_$date.pdf"))
