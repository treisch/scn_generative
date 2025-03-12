module temp_vat_analytics
# import Pkg; Pkg.add("Plots")
using Plots

export my_pdf, my_ccdf, my_ccdf2, calc_ccdf



function my_pdf(dic1, dic2, label1="vec 1", label2="vec 2", xlab="x lab")    
    x = collect(keys(dic1));
    y = getindex.(Ref(dic1), x);
    plot(x.+1,y, seriestype=:scatter, xaxis=:log10, yaxis=:log10, label=label1, xlabel=xlab, ylabel="counts")

    x2 = collect(keys(dic2));
    y2 = getindex.(Ref(dic2), x2);
    plot!(x2.+1, y2, seriestype=:scatter, xaxis=:log10, yaxis=:log10, label=label2, xlabel=xlab, ylabel="counts")
end

function calc_ccdf(dic)
    x = sort(collect(keys(dic)))
    y = getindex.(Ref(dic), x);
    y = cumsum(y)
    ccdf = 1 .- (y ./ y[end])
    return x,ccdf
end

function my_ccdf2(dic1, dic2, label1="vec 1", label2="vec 2", xlab="x lab")    
    x1,ccdf1 = calc_ccdf(dic1)
    plot(x1[1:end-1] .+1, ccdf1[1:end-1], xaxis=:log10, yaxis=:log10, label=label1, xlabel=xlab, ylabel="CCDF, p(X > x)")
    
    x2,ccdf2 = calc_ccdf(dic2)
    plot!(x2[1:end-1] .+1, ccdf2[1:end-1], xaxis=:log10, yaxis=:log10, label=label2)
end

function my_ccdf(dic1, dic2, label1="vec 1", label2="vec 2", xlab="x lab")    
    x1,ccdf1 = calc_ccdf(dic1)
    pushfirst!(ccdf1, 1)
    plot(x1 .+1, ccdf1[1:end-1], xaxis=:log10, yaxis=:log10, label=label1, xlabel=xlab, ylabel="CCDF, p(X ≥ x)")
    
    x2,ccdf2 = calc_ccdf(dic2)
    pushfirst!(ccdf2,1)
    plot!(x2 .+1, ccdf2[1:end-1], xaxis=:log10, yaxis=:log10, label=label2)
end


end