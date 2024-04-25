## --- Set up 
    # Relationship between d13C and H/C ratios can be modeled with a Rayleigh-style model 
    # of fractionation  

    # Packages
    using StatGeochem, Plots
    using LsqFit

    include("utilities/Utilities.jl")

    # Load data 
    data = importdataset("data/compilation.csv", ',', importas=:Tuple)

    
## --- [PLOT] Relationship between H/C and d13c 
    # Get data and exclude some weird outliers, or don't
    t = @. !isnan(data.d13c_org) & !isnan(data.hc);
    # t .&= data.d13c_org .> -40;

    # Fit Rayleigh fractionation model 
    params = fit_rayleigh(data.d13c_org[t], data.hc[t])
    x,y = rayleigh_curve(params, data.hc[t])

    # Plot data
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        xlabel="Kerogen H/C ratio", ylabel="Organic δ¹³C [‰]",
        fg_color_legend=:white,
        size=(600,500),
        right_margin=(25,:px),
        ylims=(-62,2),
    )
    plot!(h, data.hc[t], data.d13c_org[t], label="", 
        zcolor=data.age[t], color=:thermal, msc=:auto,
        seriestype=:scatter,
        colorbartitle="\nSample Age [Ma.]"
    )
    plot!(h, x, y, label="Rayleigh Fractionation", 
        color=:red,
        linewidth=2,
    )
    display(h)
    savefig(h, "figures/rayleigh.pdf")


## --- End of file