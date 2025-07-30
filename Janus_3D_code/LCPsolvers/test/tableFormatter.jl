import LaTeXTabulars.SimpleCellFormatter
myFormatter = SimpleCellFormatter()
function (formatter::SimpleCellFormatter)(x)
    if isnothing(x)
        out = formatter.nothing
    elseif ismissing(x)
        out = formatter.missing
    elseif x isa Real
        if x == Inf
            out = formatter.inf
        elseif x == -Inf
            out = formatter.minus_inf
        elseif isnan(x)
            out = formatter.NaN
        else
            out = string(x)
        end
    else
        out = x
    end
    # @show x, out
    out 
end
import LaTeXTabulars.latex_cell
function latex_cell(io::IO, x, formatter)
    print(io, formatter(x))
end