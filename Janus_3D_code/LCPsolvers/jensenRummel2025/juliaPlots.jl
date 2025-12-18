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
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)"=>"B-PQN: p=3, ϵ=1e-5",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN: p=4, ϵ=1e-6",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN: p=6, ϵ=1e-6",
) 
name2LatexDisplay = Dict(
    "PGD"=>L"\mathrm{PGD}",
    "Accelerated PGD"=>L"\mathrm{A-PGD}",
    "zeroSR1"=>L"\mathrm{zeroSR1}",
    "L-BFGS-B"=>L"\mathrm{L-BFGS-B}",
    "PQN"=>L"\mathrm{PQN}",
    "Min-Map Newton"=>L"\mathrm{Min-Map}\;\mathrm{Newton}",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)"=>"B-PQN: p=3, "*L"\epsilon_\mathrm{gmres}=1e-5",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN: p=4, "*L"\epsilon_\mathrm{gmres}=1e-6",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN: p=6, "*L"\epsilon_\mathrm{gmres}=1e-6",
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
    algoNames,
    metricName,
    prefix,
    root;
    fig_dir = "/Users/niru8088/scratch/Spheroidal3D-collisions/docs/fig",
    rel_tol_kkt = 1e-8,
    abs_tol_kkt = 1e-8,
    nMin = 100 ,
    nMax = 150,
)
    matdic = matread(joinpath(root, "$prefix.mat"))
    _results = matdic["results"]
    mcGood = findall(isa.(_results["matVecs"][:,1], Real) .&& (0 .< _results["matVecs"][:,1]))
    # Int.(matdic["mcGood"])[:]
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
    # display(p)
    # open(joinpath(fig_dir, "$(prefix)_caption.tex"), "w") do f
    #     s = """
    #     All the LCP's with a problem of size of $problem_size_string were ran by all algorithms. There were a total of $num_problems problems. The tolerance of 1E-6 was used for both the absolute and relative kkt condition.
    # """
    #     write(f, s)
    # end
    boxPlotFile =  joinpath(fig_dir, "$(prefix)_$(metricName)_boxPlot.pdf")
    @info "Saveing box plot to $boxPlotFile"
    PlotlyJS.savefig(
        p,
        boxPlotFile,
        height=650,
        width=800
    )
    ##
    rows = Any[]
    push!(rows, ["", "Minimum", "Lower Quartile", "Median", "Mean", "Upper Quartile", "Maximum"])
        push!(rows, Rule(:top))
    list_minimum = zeros(length(results))
    list_low_quantile = zeros(length(results))
    list_median = zeros(length(results))
    list_mean = zeros(length(results))
    list_up_quantile = zeros(length(results))
    list_maximum = zeros(length(results))
    for (i,name) in enumerate(algoNames)
        metric = results[name][metricName]
        list_minimum[i] = minimum(metric) 
        list_low_quantile[i] = quantile(metric, 0.25) 
        list_median[i] = median(metric) 
        list_mean[i] = mean(metric) 
        list_up_quantile[i] = quantile(metric, 0.75) 
        list_maximum[i] = maximum(metric) 
    end
    for (i,name) in enumerate(algoNames)
        metric = results[name][metricName]
        if contains(name, "\\kappa")
            prts = split(name, "\\kappa")
            name = string("\\text{"*prts[1]*"}\\kappa"*prts[end])
            println(name )
            name = replace(name, 
                "\\kappa"=>raw"\kappa", "\\eta"=>raw"\eta", "\\tau"=>raw"\tau", "\\_"=>raw"\_")
            name = LaTeXString("\$"*name*"\$")
        end
        row = [LaTeXString(name2LatexDisplay[name])]
        push!(row, list_minimum[i] == minimum(list_minimum) ? 
            LaTeXString(@sprintf("\\textbf{%.4g}",minimum(metric))) : 
            @sprintf("%.4g",list_minimum[i]))
        push!(row, list_low_quantile[i] == minimum(list_low_quantile) ? 
            LaTeXString(@sprintf("\\textbf{%.4g}",list_low_quantile[i])) :
            @sprintf("%.4g",list_low_quantile[i]))
        push!(row, list_median[i] == minimum(list_median) ? 
            LaTeXString(@sprintf("\\textbf{%.4g}",list_median[i])) :
            @sprintf("%.4g",list_median[i]))
        push!(row, list_mean[i] == minimum(list_mean) ? 
            LaTeXString(@sprintf("\\textbf{%.4g}",list_mean[i])) :
            @sprintf("%.4g",list_mean[i]))
        push!(row, list_up_quantile[i] == minimum(list_up_quantile) ? 
            LaTeXString(@sprintf("\\textbf{%.4g}",list_up_quantile[i])) :
            @sprintf("%.4g",list_up_quantile[i]))
        push!(row, list_maximum[i] == minimum(list_maximum) ? 
            LaTeXString(@sprintf("\\textbf{%.4g}",list_maximum[i])) :
            @sprintf("%.4g",list_maximum[i]))
        push!(rows, row)
    end
    push!(rows, Rule(:bottom))

    table_file = joinpath(fig_dir, "$(prefix)_$(metricName)_table.tex")
    @info "Saving table to $table_file"
    latex_tabular(
        table_file, 
        Tabular("lcccccc"), 
        rows; formatter=myFormatter
    )
end
## Estimated Time
# monodisperse (old)
root = "/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/goodData";
prefix = "amphi.lcp.lattice.n_5.p_8.cDist_2.5";
# polydisperse
# root = "/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord"
# prefix = "amphi.lcp.lattice.n_5.p_8.cDist_3.lcpSlvr_proxquasinewton.polyDisperseRatio_0.2"
#
monoNames = [
    "PGD",
    "Accelerated PGD",
    "zeroSR1",
    "L-BFGS-B",
    "PQN",
    "Min-Map Newton",
]
bifiNames = [
    "PGD",
    "PQN",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"
]
createBoxPlotAndTable(monoNames,"estimTime", "$prefix.mono", root)
createBoxPlotAndTable(monoNames,"estimTime", "$prefix.warmStart", root)
createBoxPlotAndTable(bifiNames, "estimTime", "$prefix.bifi", root)
createBoxPlotAndTable(bifiNames, "estimTime", "$prefix.warmStart.bifi", root)
## MatVecs
createBoxPlotAndTable(monoNames,"matVecs", "$prefix.mono", root)
createBoxPlotAndTable(monoNames,"matVecs", "$prefix.warmStart", root)
createBoxPlotAndTable(bifiNames, "eMatVecs", "$prefix.bifi", root)
createBoxPlotAndTable(bifiNames, "eMatVecs", "$prefix.warmStart.bifi", root)