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
    "BB-PGD"=>"BB-PGD",
    "A-PGD"=>"A-PGD",
    "zeroSR1"=>"zeroSR1",
    "L-BFGS-B"=>"L-BFGS-B",
    "Monofidelity PQN"=>"Mono-PQN",
    "Min-Map Newton"=>"Min-Map Newton",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN: p=3, ϵ=1e-5",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN: p=4, ϵ=1e-6",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN: p=6, ϵ=1e-6",
) 
name2LatexDisplay = Dict(
    "BB-PGD"=>"BB-PGD",
    "A-PGD"=>"A-PGD",
    "zeroSR1"=>"zeroSR1",
    "L-BFGS-B"=>"L-BFGS-B",
    "Min-Map Newton"=>"Min-Map Newton",
    "Monofidelity PQN"=>"Mono-PQN",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN : "*L"p=3, \epsilon_\mathrm{gmres}=10^{-6}",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN : "*L"p=4, \epsilon_\mathrm{gmres}=10^{-6}",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>"B-PQN : "*L"p=6, \epsilon_\mathrm{gmres}=10^{-6}",
) 

name2Color = Dict(
    "BB-PGD"=>colorant"#1f77b4",  # muted blue
    "A-PGD"=>colorant"#ff7f0e",  # safety orange
    "zeroSR1"=>colorant"#2ca02c",  # cooked asparagus green
    "L-BFGS-B"=>colorant"#d62728",  # brick red
    "Monofidelity PQN"=>colorant"#9467bd",  # muted purple
    "Min-Map Newton"=>colorant"#8c564b",  # chestnut brown
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>colorant"#e377c2",  # raspberry yogurt pink
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>colorant"#7f7f7f",  # middle gray
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"=>colorant"#bcbd22",  # curry yellow-green
)
font_size = 35
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
        ## Fix the matVec, mvps & emvps so that we have the convergence at the proper
        for (i,errHist) in enumerate(_results["errHist"][mcGood,ixAlgo])
            errHist[1,2] = Inf
            ix = findfirst((errHist[:,1] .< rel_tol_kkt) .|| (errHist[:,2] .< abs_tol_kkt))
            results[name]["matVecs"][i] = errHist[ix,3]
            results[name]["eMatVecs"][i] = errHist[ix,4]
            results[name]["estimTime"][i] = errHist[ix,end]
        end
    end
    ## Make Plots
    p = make_subplots(rows=1, cols=2, shared_yaxes=true,) # subplot_titles=[metrics[1] == "eMatVecs" ? "Number of Effective MVPs" : "Number of MVPs" "Estimated Time"])
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
                text="Effective MVPs",
                font_size=font_size
            ),
            tickfont_size = font_size,
            dtick=1,
            range=log10.([2,200]),
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
    push!(rows, Rule(:top))
    push!(rows, ["", L"\text{\large Minimum}",
    #  "Lower Quartile", 
     L"\text{\large Median}", L"\text{\large Mean}", 
    #  "Upper Quartile",
     L"\text{\large Maximum}"])
    push!(rows, Rule(:mid))
    push!(rows, [L"\text{\emph{Baseline Methods}}", "", "", "", "", "", ""])
    emvps_minimum = zeros(length(results))
    emvps_low_quantile = zeros(length(results))
    emvps_median = zeros(length(results))
    emvps_mean = zeros(length(results))
    emvps_up_quantile = zeros(length(results))
    emvps_maximum = zeros(length(results))
    mvps_minimum = zeros(length(results))
    mvps_low_quantile = zeros(length(results))
    mvps_median = zeros(length(results))
    mvps_mean = zeros(length(results))
    mvps_up_quantile = zeros(length(results))
    mvps_maximum = zeros(length(results))
    for (i,name) in enumerate(algoNames)
        emvps = results[name]["eMatVecs"]
        emvps_minimum[i] = minimum(emvps) 
        emvps_low_quantile[i] = quantile(emvps, 0.25) 
        emvps_median[i] = median(emvps) 
        emvps_mean[i] = mean(emvps) 
        emvps_up_quantile[i] = quantile(emvps, 0.75) 
        emvps_maximum[i] = maximum(emvps) 
        mvps = results[name]["matVecs"]
        mvps_minimum[i] = minimum(mvps) 
        mvps_low_quantile[i] = quantile(mvps, 0.25) 
        mvps_median[i] = median(mvps) 
        mvps_mean[i] = mean(mvps) 
        mvps_up_quantile[i] = quantile(mvps, 0.75) 
        mvps_maximum[i] = maximum(mvps) 
    end
    seenPQN = false
    for (i,name) in enumerate(algoNames)
        metric = results[name][metrics[1]]
        if !seenPQN && contains(name,"PQN")
            seenPQN = true
            push!(rows, Rule(:mid))
            push!(rows, [L"\text{\emph{Proposed Methods}}", "", "", "", "", "", ""])
        end
        row = [LaTeXString(name2LatexDisplay[name])]
        if emvps_minimum[i] == mvps_minimum[i]
            push!(row, @sprintf("%.2g",emvps_minimum[i]))
            # push!(row, @sprintf("%.2g",emvps_low_quantile[i]))
            push!(row, @sprintf("%.2g",emvps_median[i]))
            push!(row, @sprintf("%.2g",emvps_mean[i]))
            # push!(row, @sprintf("%.2g",emvps_up_quantile[i]))
            push!(row, @sprintf("%.2g",emvps_maximum[i]))
        else
            push!(row, @sprintf("%.2g (%.2g)",emvps_minimum[i], mvps_minimum[i]))
            # push!(row, @sprintf("%.2g (%.2g)",emvps_low_quantile[i], mvps_low_quantile[i]))
            push!(row, @sprintf("%.2g (%.2g)",emvps_median[i], mvps_median[i]))
            push!(row, @sprintf("%.2g (%.2g)",emvps_mean[i], mvps_mean[i]))
            # push!(row, @sprintf("%.2g (%.2g)",emvps_up_quantile[i], mvps_up_quantile[i]))
            push!(row, @sprintf("%.2g (%.2g)",emvps_maximum[i], mvps_maximum[i]))
        end
        if contains(name,"B-PQN") && contains(name, "p=4")
            for i = 1:length(row)
                row[i] = "\\textbf{$(row[i])}"
            end
        end
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
# Paths to offline results 
root = "/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/offlineResults";
prefix = "amphi.lcp.lattice.n_5.p_8.cDist_2.5";
postfix = "final"
#
algoNames = [
    "BB-PGD",
    "A-PGD",
    "zeroSR1",
    "L-BFGS-B",
    "Min-Map Newton",
    "Monofidelity PQN",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)",
    "B-PQN"*L"\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)"
]
## final 
p = createBoxPlotAndTable(algoNames,["eMatVecs", "estimTime"], "$prefix.$postfix", root, fig_dir, "LCP Solver Comparison")
nothing