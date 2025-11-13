## Dependencies
using PlotlyKaleido: restart as restart_kaleido
restart_kaleido(plotly_version = "2.35.2", mathjax = true) 
using Latexify, LaTeXStrings, PlotlyJS, LaTeXTabulars, Statistics, Colors
Latexify.set_default(fmt = "%.4g")
using MAT: matread
include("tableFormatter.jl")
##
metricName = "matVec"
rel_tol_kkt = 1e-8
abs_tol_kkt = 1e-8
nMin = 100 
nMax = 150
fname = "amphi.lattice.n_5.p_8.cDist_2.5.warmStart"
matdic = matread("/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/$fname.mat")
_results = matdic["results"]
mcGood = Int.(matdic["mcGood"])[:]
algoNames = _results["name"][mcGood[1], :]
results = Dict()
for (ixAlgo, name) in enumerate(algoNames)
    if name == "CVX"
        continue
    end
    results[name] = Dict()
    for (k, v) in _results 
        if k == "name" || k == "algo"
            continue 
        end
        results[name][k] = v[mcGood,ixAlgo]
    end
end
##
PLOTLYJS_COLORS = [
    colorant"#1f77b4",  # muted blue
    colorant"#ff7f0e",  # safety orange
    colorant"#2ca02c",  # cooked asparagus green
    colorant"#d62728",  # brick red
    colorant"#9467bd",  # muted purple
    colorant"#8c564b",  # chestnut brown
    colorant"#e377c2",  # raspberry yogurt pink
    colorant"#7f7f7f",  # middle gray
    colorant"#bcbd22",  # curry yellow-green
    colorant"#17becf"   # blue-teal
]
trs = AbstractTrace[]
for (kk,name) in enumerate(algoNames)
    println(name)
    metric = []
    for errHist in results[name]["errHist"]
        errHist[1,2] = Inf
        ix = findfirst((errHist[:,1] .< rel_tol_kkt) .|| (errHist[:,2] .< abs_tol_kkt))
        if !isnothing(ix)
            # number of iterations is the index found
            # number of matVecs is the third column of the errHist
            push!(
                metric, 
                metricName == "matVec" ? errHist[ix,3] : ix
            )
        else
             push!(
                metric, 
                length(errHist[:,1]), 
            )
        end
    end
    results[name][metricName] = metric
    if contains(name, "\\kappa")
        prts = split(name, "\\kappa")
        name = "\$\\text{"*prts[1]*"}\\kappa"*prts[end]*"\$"
    end
    println(name)
    push!(
        trs, 
        box(
            name=name,
            x=metric,
            marker_color=PLOTLYJS_COLORS[kk]
        )
    )
end
problem_size_string = "\$"#"\$n \\in [$(nMin),$(nMax)), "
num_problems = length(mcGood)

p = plot(
    trs,
    Layout(
        title=attr(
            text="Number of $( metricName == "matVec" ? "Matrix Vector Products" : "Iterations")",
            x=0.5,
            xanchor="center",
        ),
        xaxis_type="log",
        showlegend=false,
        annotations=[attr(
            text="$(problem_size_string)\\varepsilon_\\text{rel} = $(rel_tol_kkt), \\varepsilon_\\text{abs} = $(abs_tol_kkt)\$",
            font=attr(
                size= 13, # Adjust font size as needed
                color= "rgb(116, 101, 130)" # Set subtitle color
            ),
            showarrow= false, # Hide the arrow associated with annotations
            align= "center", # Align the subtitle horizontally
            x= 0.35, # Center horizontally (0.5 for paper reference)
            y= 1.05, # Position at the top (1 for paper reference)
            xref= "paper", # Reference coordinates to the plot paper
            yref= "paper" # Reference coordinates to the plot paper
        )]
    )
)
display(p)
nothing
##
save_dir = "/Users/niru8088/scratch/Spheroidal3D-collisions/docs/fig"
# open(joinpath(save_dir, "$(fname)_caption.tex"), "w") do f
#     s = """
#     All the LCP's with a problem of size of $problem_size_string were ran by all algorithms. There were a total of $num_problems problems. The tolerance of 1E-6 was used for both the absolute and relative kkt condition.
# """
#     write(f, s)
# end
PlotlyJS.savefig(
    p,
    joinpath(save_dir, "$(fname)_boxPlot.pdf"),
    height=600,
    width=800
)
##
rows = Any[]
push!(rows, ["", "Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum"])
    push!(rows, Rule(:top))
for name in algoNames
    metric = results[name][metricName]
    if contains(name, "\\kappa")
        prts = split(name, "\\kappa")
        name = string("\\text{"*prts[1]*"}\\kappa"*prts[end])
        println(name )
        name = replace(name, 
            "\\kappa"=>raw"\kappa", "\\eta"=>raw"\eta", "\\tau"=>raw"\tau", "\\_"=>raw"\_")
        name = LaTeXString("\$"*name*"\$")
    end
    push!(rows, [name, minimum(metric), quantile(metric,0.25), quantile(metric,0.5), quantile(metric, 0.75), maximum(metric)])
end
push!(rows, Rule(:bottom))

table_file = joinpath(save_dir, "$(fname)_table.tex")
@info "Saving table to $table_file"
latex_tabular(
    table_file, 
    Tabular("lccccc"), 
    rows; formatter=myFormatter
)
