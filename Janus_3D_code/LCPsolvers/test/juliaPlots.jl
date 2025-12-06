## Dependencies
using Printf
using PlotlyKaleido: restart as restart_kaleido
restart_kaleido(plotly_version = "2.35.2", mathjax = true) 
using Latexify, LaTeXStrings, PlotlyJS, LaTeXTabulars, Statistics, Colors
Latexify.set_default(fmt = "%.4g")
using MAT: matread
include("tableFormatter.jl")
##
name2Display = Dict(
    "PGD"=>L"\mathrm{PGD}",
    "Accelerated PGD"=>L"\mathrm{A-PGD}",
    "zeroSR1"=>L"\mathrm{zeroSR1}",
    "L-BFGS-B"=>L"\mathrm{L-BFGS-B}",
    "PQN"=>L"\mathrm{PQN}",
    "Min-Map Newton"=>L"\mathrm{Min-Map}\;\mathrm{Newton}",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)"=>L"\mathrm{B-PQN}\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>L"\mathrm{B-PQN}\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>L"\mathrm{B-PQN}\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)",
) 

name2Color = Dict(
    "PGD"=>colorant"#1f77b4",  # muted blue
    "Accelerated PGD"=>colorant"#ff7f0e",  # safety orange
    "zeroSR1"=>colorant"#2ca02c",  # cooked asparagus green
    "L-BFGS-B"=>colorant"#d62728",  # brick red
    "PQN"=>colorant"#9467bd",  # muted purple
    "Min-Map Newton"=>colorant"#8c564b",  # chestnut brown
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)"=>colorant"#e377c2",  # raspberry yogurt pink
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>colorant"#7f7f7f",  # middle gray
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>colorant"#bcbd22",  # curry yellow-green
)
# colorant"#17becf"   # blue-teal
##
function createBoxPlotAndTable(
    metricName,
    prefix;
    save_dir = "/Users/niru8088/scratch/Spheroidal3D-collisions/docs/fig",
    data_dir = "/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data",
    rel_tol_kkt = 1e-8,
    abs_tol_kkt = 1e-8,
    nMin = 100 ,
    nMax = 150,
)
    matdic = matread(joinpath(data_dir, "$prefix.mat"))
    _results = matdic["results"]
    mcGood = findall(isa.(_results["matVecs"][:,1], Real) .&& (0 .< _results["matVecs"][:,1]))
    # Int.(matdic["mcGood"])[:]
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
    trs = AbstractTrace[]
    for (kk,name) in enumerate(algoNames)
        println(name)
        # metric = []
        # for errHist in results[name]["errHist"]
        #     errHist[1,2] = Inf
        #     ix = findfirst((errHist[:,1] .< rel_tol_kkt) .|| (errHist[:,2] .< abs_tol_kkt))
        #     if !isnothing(ix)
        #         # number of iterations is the index found
        #         # number of matVecs is the third column of the errHist
        #         push!(
        #             metric, 
        #             metricName == "matVec" ? errHist[ix,3] : ix
        #         )
        #     else
        #          push!(
        #             metric, 
        #             length(errHist[:,1]), 
        #         )
        #     end
        # end
        # results[name][metricName] = metric
        # if contains(name, "\\kappa")
        #     prts = split(name, "\\kappa")
        #     name = "\$\\text{"*prts[1]*"}\\kappa"*prts[end]*"\$"
        # end
        println(name)
        push!(
            trs, 
            box(
                name=name2Display[name],
                x=results[name][metricName],
                marker_color=name2Color[name]
            )
        )
    end
    # problem_size_string = "\$"#"\$n \\in [$(nMin),$(nMax)), "
    # num_problems = length(mcGood)
    (title_text, xaxis_title_text,xaxis_range,xaxis_tickvals,xaxis_type) = if metricName == "matVecs" 
        "Number of MVPs", "MVPs", log10.([.9,155]),nothing,"log"
    elseif metricName == "eMatVecs"
        "Number of Effective MVPs", "Effective MVPs", log10.([.9,20]),nothing,"log"
    elseif metricName == "iters" 
        "Number of Iterations", "Iterations"
    elseif metricName == "estimTime"
        "Estimated Wall Time", "Seconds",log10.([1e2,1e5]),nothing,"log"
    end
    font_size =30
    p = plot(
        trs,
        Layout(
            title=attr(
                text=title_text,
                font_size=40,
                x=0.5,
                xanchor="center",
            ),
            yaxis_tickfont_size = font_size,
            xaxis_tickfont_size = font_size,
            xaxis=attr(
                title=attr(
                    text=xaxis_title_text,
                    font_size=font_size
                ),
                dtick=1,
                range=xaxis_range,
                tickvals=xaxis_tickvals,
                type=xaxis_type
            ),
            showlegend=false,
            # annotations=[attr(
            #     text="$(problem_size_string)\\varepsilon_{\\text{kkt}_\\text{rel}} = $(rel_tol_kkt), \\varepsilon_{\\text{kkt}_\\text{abs}} = $(abs_tol_kkt)\$",
            #     font=attr(
            #         size=font_size, # Adjust font size as needed
            #         color= "rgb(116, 101, 130)" # Set subtitle color
            #     ),
            #     showarrow= false, # Hide the arrow associated with annotations
            #     align= "center", # Align the subtitle horizontally
            #     x= 0.35, # Center horizontally (0.5 for paper reference)
            #     y= 1.05, # Position at the top (1 for paper reference)
            #     xref= "paper", # Reference coordinates to the plot paper
            #     yref= "paper" # Reference coordinates to the plot paper
            # )],
            # uniformtext_minsize=font_size,
        )
    )
    display(p)
    # open(joinpath(save_dir, "$(prefix)_caption.tex"), "w") do f
    #     s = """
    #     All the LCP's with a problem of size of $problem_size_string were ran by all algorithms. There were a total of $num_problems problems. The tolerance of 1E-6 was used for both the absolute and relative kkt condition.
    # """
    #     write(f, s)
    # end
    boxPlotFile =  joinpath(save_dir, "$(prefix)_$(metricName)_boxPlot.pdf")
    @info "Saveing box plot to $boxPlotFile"
    PlotlyJS.savefig(
        p,
        boxPlotFile,
        height=650,
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
        push!(rows, [LaTeXString(name),  
            @sprintf("%.4g",minimum(metric)), 
            @sprintf("%.4g",quantile(metric,0.25)), 
            @sprintf("%.4g",median(metric)), 
            @sprintf("%.4g",quantile(metric, 0.75)),  
            @sprintf("%.4g",maximum(metric))])
    end
    push!(rows, Rule(:bottom))

    table_file = joinpath(save_dir, "$(prefix)_$(metricName)_table.tex")
    @info "Saving table to $table_file"
    latex_tabular(
        table_file, 
        Tabular("lccccc"), 
        rows; formatter=myFormatter
    )
end
## Estimated Time
createBoxPlotAndTable("estimTime", "amphi.lattice.n_5.monoFidelity" #=prefix=#)
createBoxPlotAndTable("estimTime", "amphi.lattice.n_5.warmStart" #=prefix=#)
createBoxPlotAndTable("estimTime", "amphi.lattice.n_5.bifi" #=prefix=#)
createBoxPlotAndTable("estimTime", "amphi.lattice.n_5.warmStart.bifi" #=prefix=#)
## MatVecs
createBoxPlotAndTable("matVecs", "amphi.lattice.n_5.monoFidelity" #=prefix=#)
createBoxPlotAndTable("matVecs", "amphi.lattice.n_5.warmStart" #=prefix=#)
createBoxPlotAndTable("eMatVecs", "amphi.lattice.n_5.bifi" #=prefix=#)
createBoxPlotAndTable("eMatVecs", "amphi.lattice.n_5.warmStart.bifi" #=prefix=#)