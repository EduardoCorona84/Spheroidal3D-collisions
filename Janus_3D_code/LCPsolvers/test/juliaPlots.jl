## Dependencies
using PlotlyKaleido: restart as restart_kaleido
restart_kaleido(plotly_version = "2.35.2", mathjax = true) 
using Latexify, LaTeXStrings, PlotlyJS, LaTeXTabulars, Statistics
Latexify.set_default(fmt = "%.4g")
using MAT: matread
include("tableFormatter.jl")
##
metricName = "matVec"
rel_tol_kkt = 1e-8
abs_tol_kkt = 1e-8
nMin = 100 
nMax = 150
fname = "randProblems_diagDom_n_$(nMin)_$(nMax)"
matdic = matread("/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/results_$(fname).mat")
_results = matdic["results"]
mcGood = Int.(matdic["mcGood"])[:]
algoNames = _results["name"]
results = Dict()
for (ix, name) in enumerate(algoNames)
    if name == "CVX"
        continue
    end
    results[name] = Dict()
    for (k, v) in _results 
        if k == "name" || k == "algo"
            continue 
        end
        # try
        results[name][k] = v[ix][:][mcGood]
        # catch 
        #     @show k
        # end
    end
end
##
trs = AbstractTrace[]
for name in algoNames
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
        end
    end
    results[name][metricName] = metric
    if contains(name, "\\kappa")
        prts = split(name, "\\kappa")
        name = "\$\\text{"*prts[1]*"}\\kappa"*prts[end]*"\$"
    end
    push!(
        trs, 
        box(
            name=name,
            x=metric
        )
    )
end
problem_size_string = "\$n \\in [$(nMin),$(nMax))\$"
num_problems = length(mcGood)
save_dir = "/Users/niru8088/scratch/Spheroidal3D-collisions/docs/fig"
open(joinpath(save_dir, "$(fname)_caption.tex"), "w") do f
    s = """
    All the LCP's with a problem of size of $problem_size_string were ran by all algorithms. There were a total of $num_problems problems. The tolerance of 1E-6 was used for both the absolute and relative kkt condition.
"""
    write(f, s)
end

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
            text="$(problem_size_string[1:end-1]), \\varepsilon_\\text{rel} = $(rel_tol_kkt), \\varepsilon_\\text{abs} = $(abs_tol_kkt)\$",
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
savefig(
    p,
    joinpath(save_dir, "$(fname)_boxPlot.pdf"),
    height=600,
    width=800
)
display(p)
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
