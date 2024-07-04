# example/levatai/fig2.jl

using Dates
using DelimitedFiles
using GeologicTime    # 0.1.0
using HypothesisTests # 0.11.0
using JointPoint      # 0.0.5
using Loess           # 0.6.3
using Printf
using PyPlot          # 2.11.2
using Statistics

file(pathcomp...) = joinpath(ENV["LEVATAI_PATH"], pathcomp...)

function readxy(path::AbstractString; maxx=Inf)
	x, y = eachcol(readdlm(path))
	id = x .<= maxx
	return x[id], y[id]
end

function readjpxy(path::AbstractString; k=3)
	njp, x, y = eachcol(readdlm(path))
	id = njp .== k
	return x[id], y[id]
end

xc = LinRange(0, 9, 101)
yc = sin.((xc ./ 9) .^ 2 .* 9) ./ exp.(xc ./ 10)
xj, yj = findjoint(xc, yc, 3)

x, ae = readxy(file("Ages", "TaiwanEndemism_divergence_time", "All_endemism_lar.txt"))
jpxae, jpyae = readjpxy(file("Ages", "TaiwanEndemism_divergence_time", "All_endemism_jp.txt"))

x, dp = readxy(file("Ages", "TaiwanEndemism_dispersal_node_time_lar.txt"))
jpxdp, jpydp = readjpxy(file("Ages", "TaiwanEndemism_dispersal_node_time_jp.txt"))

x, se = readxy(file("Factors", "SeaElev_10M.tsv"))
model = loess(x, se; span=0.1)
sel = predict(model, x)

dy = 30
y3toy4(y) = y + dy
y4toy3(y) = y - dy

function drawmt(ax, csvpath; xspan=1, yspan=1, xstart=0, ystart=0, kwargs...)
	mtx, mty = eachcol(readdlm(file(csvpath), ','))
	x = xstart .- xspan .* mtx
	y = ystart .+ yspan .* mty
	ax.plot(x, y; lw=0.8, c="k", kwargs...)
	i = findfirst(isnan.(x))
	ax.fill_between(x, 0, y; facecolor="#93c69c")
	return ax
end

function reg(ax, x, y, cx, cy, dcx, dcy, nx, ny)
	mx = sum(x) / length(x)
	my = sum(y) / length(y)
	b = sum((x .- mx) .* (y .- my)) / sum((x .- mx) .^ 2)
	a = my - b * mx
	h(x) = a + b * x
	ax.plot(x, y, "."; c=dcy, ms=2)
	l, r = ax.get_xlim()
	d, u = ax.get_ylim()
	ax.plot([l, r], h.([l, r]), "--k"; lw=0.8)
	ax.plot([l, l], [d, u]; c=cy, lw=2.4)
	ax.plot([l, r], [u, u]; c=cx, lw=2.4)
	ax.plot([r, r], [d, u]; c=cy, lw=2.4)
	ax.plot([l, r], [d, d]; c=cx, lw=2.4)
	ax.plot([l, l], [d, 2\(d+u)]; c=cy, lw=2.4)
	ax.fill([l, l, r, r], [u, d, d, u], "#ffffff"; lw=0, zorder=-10)
	ax.set_xlim(l, r)
	ax.set_ylim(d, u)
	ax.set_axis_off()
	corr = cor(x, y)
	p = pvalue(CorrelationTest(x, y))
	ax.text(1.47r - 0.47l, 0.9u + 0.9d, nx; ha="center", va="center", c=dcx)
	ax.text(1.47r - 0.47l, 0.72u + 0.28d, "vs."; ha="center", va="center")
	ax.text(1.47r - 0.47l, 0.56u + 0.44d, ny; ha="center", va="center", c=dcy)
	ax.text(1.07r - 0.07l, 0.28u + 0.72d, "\$r\$ ="; ha="left", va="center")
	ax.text(1.87r - 0.87l, 0.28u + 0.72d, @sprintf("\$%.4f\$", corr); ha="right", va="center")
	if p < 0.0001
		ax.text(1.07r - 0.07l, 0.1u + 0.9d, "\$p <\$"; ha="left", va="center")
		ax.text(1.87r - 0.87l, 0.1u + 0.9d, "0.0001"; ha="right", va="center")
	end
	return ax, corr, p
end

close()
figure(figsize=(6.4, 6.4))

try axb.remove() catch; end
axb = PyPlot.axes([0, 0, 1, 1])
axb.set_facecolor("None")
axb.set_axis_off()
axb.set_xlim(0, 1)
axb.set_ylim(0, 1)
axb.fill([0.942, 0.942, 0.989, 0.989], [0.636, 0.791, 0.791, 0.636], "#f26a31"; lw=0, alpha=0.2)
axb.fill([0.942, 0.942, 0.989, 0.989], [0.636, 0.450, 0.450, 0.636], "#6630b2"; lw=0, alpha=0.2)
axb.fill([0.013, 0.013, 0.079, 0.079], [0.455, 0.054, 0.054, 0.455], "#66a7e6"; lw=0, alpha=0.3)

try ax1.remove() catch; end
ax1 = PyPlot.axes(([0.09, 0.85, 0.84, 0.13]))
ax1.set_xlim(10, 0)
ax1.set_ylim(0, 0.8)
drawmt(ax1, "XY_LowMt.csv"; xstart=7.4, xspan=1.5, ystart=-0.2, yspan=0.4)
drawmt(ax1, "XY_LowMt.csv"; xstart=5.6, xspan=1.5, ystart=-0.03, yspan=0.4)
drawmt(ax1, "XY_HiMt.csv"; xstart=3.6, xspan=1.7, ystart=-0.1, yspan=0.9)
drawmt(ax1, "XY_FullMt.csv"; xstart=1.7, xspan=1.7, ystart=-0.08, yspan=0.9)
ax1.plot([10, 0], [0, 0], "#66a7e6"; lw=8)
ax1.set_xticks([])
ax1.set_yticks([])

try ax2.remove() catch; end
ax2 = PyPlot.axes([0.09, 0.80, 0.84, 0.04])
drawtimescale(ax2, 10, 0, [4]; facealpha=0.2, 
	texts = Dict("Miocene" => "Miocene", "Pliocene" => "Pliocene", "Pleistocene" => "Pleistocene"))

try ax4.remove() catch; end
ax4 = PyPlot.axes([0.09, 0.07, 0.84, 0.72])
ax4.plot(x, y3toy4.(ae), "#f26a31"; label="LAR, species endemic to Taiwan", zorder=5, lw=3)
ax4.plot(jpxae, y3toy4.(jpyae), "*:"; label="          Joint points of the LAR above", c="#b45229", zorder=9, ms=10)
ax4.plot(x, y3toy4.(dp), "#6630b2"; label="LAR, lineages dispersed into Taiwan ", zorder=7) # four-per-em space " "
ax4.plot(jpxdp, y3toy4.(jpydp), "*:"; label="          Joint points of the LAR above", c="#512b86", zorder=8, ms=10)
ax4.fill([7.1, 7.1, 2.5, 2.5], [-200, 200, 200, -200], "#f1dc83"; lw=-10, zorder=0, label="Uplift of Central Mountains", alpha=0.3)
ax4.fill([1.2, 1.2, 0.7, 0.7], [-200, 200, 200, -200], "#c3f8f9"; lw=-10, zorder=0, label="Mid-Pleistocene Transition", alpha=0.5)
ax4.plot(x, se; c="#66a7e6", lw=0.8, label="Sea level, 10 ka resolution")
ax4.plot(x, sel; ls="--", c="#386b9c", label="Sea level, after LOESS")
ax4.set_ylim(-100, 135)
ax4.set_yticks(-100:10:20)
ax4.set_xlim(10, 0)
ax4.set_xticks(10:-1:0)
ax4.set_xlabel("Time (Ma)"; labelpad=0)
ax4.set_ylabel("          Sea level (m)", labelpad=-25)
ax4.legend(loc="upper left", handlelength=2.15, labelspacing=0.3)

try ax3.remove() catch; end
ax3 = ax4.twinx()
ax3.plot([10, 0], [0, 0], "--k"; lw=0.8)
ax3.set_ylim(y4toy3.(ax4.get_ylim()))
ax3.set_yticks(0:10:100, ["  0", " 10", " 20", " 30", " 40", " 50", " 60", " 70", " 80", " 90", "100"]) # figure space " "
ax3.set_ylabel("Lineage accumulation rate (Ma\$^{-1}\$)                ", labelpad=-15) # em space " "

try axr1.remove() catch; end
axr1 = PyPlot.axes([0.103, 0.083, 0.16, 0.16])
reg(axr1, se, ae, "#66a7e6", "#f26a31", "#386b9c", "#b45229", "Sea level", "LAR\$_{\\sf{endemic}}\$")

try axr2.remove() catch; end
axr2 = PyPlot.axes([0.42, 0.083, 0.16, 0.16])
reg(axr2, se, dp, "#66a7e6", "#6630b2", "#386b9c", "#512b86", "Sea level", "LAR\$_{\\sf{dispersed}}\$")

try axjp.remove() catch; end
axjp = PyPlot.axes([0.589, 0.586, 0.847-0.589, 0.777-0.586])
axjp.plot(xc, yc; c="#aaa", lw=3, label="Original")
axjp.plot(xj, yj, "*:"; c="k", ms=10, label="Joint points")
axjp.set_xlim(-1, 10)
axjp.set_ylim(-0.9, 1)
axjp.set_xticks([])
axjp.set_yticks([])
axjp.text(-0.26, -0.52, "Joint points", ha="left")
axjp.text(-0.26, -0.73, "minimizing RSS", ha="left", fontsize=8)

try ax0.remove() catch; end
ax0 = PyPlot.axes([0, 0, 1, 1])
ax0.set_facecolor("None")
ax0.set_axis_off()
ax0.set_xlim(0, 1)
ax0.set_ylim(0, 1)
ax0.text(0.103, 0.956, "(a)"; ha="left", va="center")
ax0.text(0.103, 0.526, "(b)"; ha="left", va="center")
ax0.text(0.602, 0.753, "(c)"; ha="left", va="center")
ax0.text(0.103, 0.07+(0.72/235*130)-0.024, "(d)"; ha="left", va="center")
ax0.text(0.250, 0.219, "(e)"; ha="right", va="center")
ax0.text(0.567, 0.219, "(f)"; ha="right", va="center")

date = Dates.format(today(), "YYYYmmdd")
savefig(file("Fig2_$date.pdf"))
