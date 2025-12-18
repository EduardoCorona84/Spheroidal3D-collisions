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
    "PGD"=>"PGD",
    "Accelerated PGD"=>"A-PGD",
    "zeroSR1"=>"zeroSR1",
    "L-BFGS-B"=>"L-BFGS-B",
    "PQN"=>"PQN",
    "Min-Map Newton"=>"Min-Map Newton",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)"=>"B-PQN: p=3, ϵ=1e-5",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN: p=4, ϵ=1e-6",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN: p=6, ϵ=1e-6",
) 
name2LatexDisplay = Dict(
    "PGD"=>"PGD",
    "Accelerated PGD"=>"A-PGD",
    "zeroSR1"=>"zeroSR1",
    "L-BFGS-B"=>"L-BFGS-B",
    "PQN"=>"PQN",
    "Min-Map Newton"=>"Min-Map Newton",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)"=>"  p=3,"*L"\epsilon_\mathrm{gmres}"*"=1e-5",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"  p=4,"*L"\epsilon_\mathrm{gmres}"*"=1e-6",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"  p=6,"*L"\epsilon_\mathrm{gmres}"*"=1e-6",
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
font_size =35
# colorant"#17becf"   # blue-teal
##
function _createBoxPlot(p,r,c,algoNames,metricName, results)
    
    for name in reverse(algoNames)
        add_trace!(p,
            box(
                name=name2Display[name],
                x=results[name][metricName],
                marker_color=name2Color[name],
                showlegend=false
            ), row=r, col=c
        )
    end
  
   return results
end

function createBoxPlotAndTable(algoNames,metrics,prefix,root,fig_dir, title_text; rel_tol_kkt = 1e-8, abs_tol_kkt = 1e-8)
    ## Get results from mat file
    matdic = matread(joinpath(root, "$prefix.mat"))
    _results = matdic["results"]
    mcGood = Int.(matdic["mcGood"])[:]
    results = Dict()
    for name in algoNames
        if name == "CVX"
            continue
        end
        results[name] = Dict()
        ixAlgo = findfirst(name .== _results["name"][mcGood[1],:])
        for (k, v) in _results 
            if k == "name" || k == "algo"
                continue 
            end
            results[name][k] = v[mcGood,ixAlgo]
        end
        ## Fix the matVec so that we have tota
        for (i,errHist) in enumerate(_results["errHist"][mcGood,ixAlgo])
            errHist[1,2] = Inf
            ix = findfirst((errHist[:,1] .< rel_tol_kkt) .|| (errHist[:,2] .< abs_tol_kkt))
            results[name]["matVecs"][i] = errHist[ix,3]
            results[name]["estimTime"][i] = errHist[ix,end]
        end
    end
    ## Make Plots
    p = make_subplots(rows=1, cols=2, shared_yaxes=true,)# subplot_titles=[metrics[1] == "eMatVecs" ? "Number of Effective MVPs" : "Number of MVPs" "Estimated Time"])
    results = _createBoxPlot(p,1,1,algoNames,metrics[1],results)
    _createBoxPlot(p,1,2,algoNames,metrics[2],results)
    relayout!(p, 
        title=attr(
            text=title_text,
            font_size=40,
            x=0.6,
            xanchor="center",
        ),
        yaxis_tickfont_size = font_size,
        xaxis=attr(
            title=attr(
                text="MVPs",
                font_size=font_size
            ),
            tickfont_size = font_size,
            dtick=1,
            range=log10.([2,50]),
            type="log"
        ),  
        xaxis2=attr(
            title=attr(
                text= "Seconds",
                font_size=font_size
            ),
            tickfont_size = font_size,
            dtick=1,
            range=log10.([1e2,1e4]),
            type="log"
        ),
        uniformtext_minsize=font_size,
        uniformtext_mode="show"
    )
    ## Save Plot 
    boxPlotFile = joinpath(fig_dir, "$(prefix)_boxPlot.pdf")
    @info "Saveing box plot to $boxPlotFile"
    PlotlyJS.savefig(
        p,
        boxPlotFile,
        height=800,
        width=1400
    )

    ## LaTex Table
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
        metric = results[name][metrics[1]]
        list_minimum[i] = minimum(metric) 
        list_low_quantile[i] = quantile(metric, 0.25) 
        list_median[i] = median(metric) 
        list_mean[i] = mean(metric) 
        list_up_quantile[i] = quantile(metric, 0.75) 
        list_maximum[i] = maximum(metric) 
    end
    seenBPQN = false
    for (i,name) in enumerate(algoNames)
        metric = results[name][metrics[1]]
        if !seenBPQN && contains(name,"B-PQN")
            seenBPQN = true
            push!(rows, ["B-PQN", "", "", "", "", "", ""])
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

    table_file = joinpath(fig_dir, "$(prefix)_$(metrics[1])_table.tex")
    @info "Saving table to $table_file"
    latex_tabular(
        table_file, 
        Tabular("lcccccc"), 
        rows; formatter=myFormatter
    )
    return p
end


## 
fig_dir = "/Users/niru8088/scratch/Spheroidal3D-collisions/docs/fig"
# monodisperse (old)
root = "/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/goodData";
prefix = "amphi.lcp.lattice.n_5.p_8.cDist_2.5";
# polydisperse
# root = "/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord"
# prefix = "amphi.lcp.lattice.n_5.p_8.cDist_3.lcpSlvr_proxquasinewton.polyDisperseRatio_0.2"
#
monoNames = ["PGD", "Accelerated PGD", "zeroSR1", "L-BFGS-B", "PQN", "Min-Map Newton",
]
bifiNames = ["PGD", "PQN", "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)", "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)", "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"
]
## Monofidelity
p1 = createBoxPlotAndTable(monoNames,["matVecs", "estimTime"], "$prefix.mono", root, fig_dir,"Monofidelity Comparison")
## warmStart
p2 = createBoxPlotAndTable(monoNames,["matVecs", "estimTime"], "$prefix.warmStart", root, fig_dir,"Warm Starting Comparison")
## bifi 
p3 = createBoxPlotAndTable(bifiNames,["eMatVecs", "estimTime"], "$prefix.bifi", root, fig_dir,"Bifidelity Comparison")
## warmStart.bifi
p4 = createBoxPlotAndTable(bifiNames,["eMatVecs", "estimTime"], "$prefix.warmStart.bifi", root, fig_dir,"Warm Starting and Bifidelity Comparison")
nothing