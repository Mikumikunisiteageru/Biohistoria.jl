# example/levatai/figsublar.jl

using Dates           # 1.10.2
using DelimitedFiles  # 1.10.2
using GeologicTime    # 0.1.0
using PyPlot          # 2.11.2, Python 3.12.3, MatPlotLib 3.8.4

file(pathcomp...) = joinpath(ENV["LEVATAI_PATH"], pathcomp...)

function readxy(path::AbstractString; maxx=Inf)
	x, y = eachcol(readdlm(path))
	id = x .<= maxx
	return x[id], y[id]
end

x, y = readxy(file("Ages", "TaiwanEndemism_divergence_time", "All_endemism_lar.txt"))
x, allopatric = readxy(file("Ages", "TaiwanEndemism_divergence_time", "endemism_allopatric_lar.txt"))
x, sympatric = readxy(file("Ages", "TaiwanEndemism_divergence_time", "endemism_sympatric_lar.txt"))
x, temp = readxy(file("Ages", "TaiwanEndemism_divergence_time", "temperate_endemism_lar.txt"))
x, trop = readxy(file("Ages", "TaiwanEndemism_divergence_time", "tropic_endemism_lar.txt"))

close()
plot(x, y)

close()
plot(x, allopatric)
plot(x, sympatric)

# plot(x, allopatric + sympatric)
close()
plot(x, temp)
plot(x, trop)
# plot(x, temp + trop)

plot(x, allopatric ./ y)
plot(x, temp ./ y)

figure();
plot(y, temp + trop, "#6cf"; lw=3)
plot([0, 100], [0, 100], "k"; lw=0.8)

m1 = maximum(y)
m2 = maximum(max.(allopatric, sympatric))
m3 = maximum(max.(temp, trop))

close()
figure(figsize=(6.4, 7.2))

try axg.remove() catch; end
axg = PyPlot.axes(([0.02, 0.95, 0.84, 0.03]))
drawtimescale(axg, 10, 0, [4]; facealpha=0.2, 
	texts = Dict("Miocene" => "Miocene", "Pliocene" => "Pliocene", "Pleistocene" => "Pleistocene"))

try ax1.remove() catch; end
ax1 = PyPlot.axes(([0.02, 0.06 + 0.88 * (m2 + m3) / (m1 + m2 + m3), 0.84, 0.88 * m1 / (m1 + m2 + m3)]))
ax1.tick_params(left=false, labelleft=false, bottom=false, labelbottom=false, right=true, labelright=true)
ax1.plot(x, y, "#f26a31"; label="All", lw=3)
ax1.set_xlim(10, 0)
ax1.set_ylim(0, m1 * 1.1)
ax1.grid(axis="x", color="#ccc")
ax1.legend(loc="upper left")

try ax2.remove() catch; end
ax2 = PyPlot.axes(([0.02, 0.05 + 0.88 * m3 / (m1 + m2 + m3), 0.84, 0.88 * m2 / (m1 + m2 + m3)]))	
ax2.tick_params(left=false, labelleft=false, bottom=false, labelbottom=false, right=true, labelright=true)
ax2.plot(x, allopatric, "goldenrod"; label="Allopatric", lw=3.5)
ax2.plot(x, sympatric, "darkkhaki"; label="Sympatric", lw=2)
ax2.set_xlim(10, 0)
ax2.set_ylim(0, m2 * 1.1)
ax2.set_yticks(0:20:40)
ax2.grid(axis="x", color="#ccc")
ax2.legend(framealpha=1)

try ax3.remove() catch; end
ax3 = PyPlot.axes(([0.02, 0.04, 0.84, 0.88 * m3 / (m1 + m2 + m3)]))
ax3.tick_params(left=false, labelleft=false, right=true, labelright=true)
ax3.plot(x, temp, "lightgreen"; label="Temperate", lw=3.5)
ax3.plot(x, trop, "deepskyblue"; label="Tropical", lw=2)
ax3.set_xlim(10, 0)
ax3.set_ylim(0, m3 * 1.1)
ax3.set_yticks(0:20:60)
ax3.grid(axis="x", color="#ccc")
ax3.legend()

try ax0.remove() catch; end
ax0 = PyPlot.axes([0, 0, 1, 1])
ax0.set_facecolor("None")
ax0.set_axis_off()
ax0.set_xlim(0, 1)
ax0.set_ylim(0, 1)
ax0.text(0.939, 0.015, "Time (Ma)"; ha="center", va="center")
ax0.text(0.945, 0.490, "Lineage accumulation rate (Ma\$^{-1}\$)"; ha="center", va="center", rotation=90)

date = Dates.format(today(), "YYYYmmdd")
savefig(file("Fig_LAR_Sub_$date.pdf"))
