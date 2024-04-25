## --- Set up
    # Try to replicate the carbon isotope curves plotted in Des Marais et al., 1992

    # Packages 
    using StatGeochem, Plots  # LsqFit

    include("utilities/Utilities.jl")

    # Load data
    data = importdataset("data/compilation.csv", ',', importas=:Tuple)

    # Restrict to only data from Strauss et al., 1992 (data available to Des Marais)
    strauss = data.data_reference .== "Schopf and Klein 1992";


## --- Plot uncorrected data over time [Fig. 1]
    xmin, xmax = 500, 2700
    nbins = Int((xmax - xmin) / 100)
    
    # Plot 
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        size=(400,600),
        xlabel="Age [Ma.]", ylabel="δ¹³C [‰]",
        xlims=(500,3000),
        xflip=true,
        # ylims=(-60,10)
        fg_color_legend=:white
    );

    # Raw data
    plot!(data.age[strauss], data.d13c_carb[strauss], label="Carbonate",
        seriestype=:scatter,
        markersize=2, markershape=:+,
        color=:blue, msc=:auto, alpha=0.5
    )
    plot!(data.age[strauss], data.d13c_org[strauss], label="Organic",
        seriestype=:scatter,
        markersize=2,
        color=:green, msc=:auto, alpha=0.5
    )

    # Moving means 
    c,m,e = binmeans(data.age[strauss], data.d13c_carb[strauss], xmin, xmax, nbins, relbinwidth=2)
    plot!(c, m, label="", color=:blue, linewidth=2)

    c,m,e = binmeans(data.age[strauss], data.d13c_org[strauss], xmin, xmax, nbins, relbinwidth=2)
    plot!(c, m, label="", color=:green, linewidth=2)

    display(h)
    

## --- Correct organic carbon for alteration using H/C ratios [Fig. 2]
    # Figure 2 contains data only for kerogen with ash contents less than 25% and H/C values 
    # greater than 0.1
    s = data.hc[strauss] .>= 0.1;

    # Shift in δ¹³C (org) using the equation from Des Marias 
    Δδ(hc::Number) = 4.05 - 3.05*hc + 0.785/hc + 0.0165/hc^2 - (8.79E-4)/hc^3
    corrected_org = data.d13c_org[strauss][s] .- Δδ.(data.hc[strauss][s])

    # Plot 
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        size=(400,600),
        xlabel="Age [Ma.]", ylabel="δ¹³C [‰]",
        xlims=(500,3000),
        xflip=true,
        fg_color_legend=:white
    )
    plot!(data.age[strauss], data.d13c_carb[strauss], label="Carbonate",
        seriestype=:scatter,
        markersize=2, markershape=:+,
        color=:blue, msc=:auto, alpha=0.5
    )

    plot!(data.age[strauss][s], corrected_org, label="Organic",
        seriestype=:scatter,
        markersize=2,
        color=:green, msc=:auto, alpha=0.5
    )

    # Moving means 
    c,m,e = binmeans(data.age[strauss], data.d13c_carb[strauss], xmin, xmax, nbins, relbinwidth=2)
    plot!(c, m, label="", color=:blue, linewidth=2)

    c,m,e = binmeans(data.age[strauss][s], corrected_org, xmin, xmax, nbins, relbinwidth=2)
    plot!(c, m, label="", color=:green, linewidth=2)

    display(h)


## --- Fraction of carbon buried as organic [Fig. 3]
    # 100 Ma bins of organic and carbonate data between 2500 and 600 Ma.
    xmin = 600
    xmax = 2500
    nbins = Int((xmax - xmin) / 100)
    q = @. xmin <= data.age[strauss] <= xmax; 

    s = data.hc[strauss] .>= 0.1;
    corrected_org = data.d13c_org[strauss][s .& q] .- Δδ.(data.hc[strauss][s .& q])

    c_carb, m_carb, e_carb = binmeans(data.age[strauss][q], data.d13c_carb[strauss][q], xmin,xmax,nbins)
    c_org, m_org, e_org = binmeans(data.age[strauss][s .& q], corrected_org, xmin,xmax,nbins)

    # Fraction of carbon buried as organic 
    frog_curve = f_org.(-5, 0, m_org)

    # Original data
    des_marais = (;
        age=[2650.0, 2495.5516180173463, 2047.2862072181294, 1949.6316329385436, 1853.1940688240234, 
            1747.3141844633037, 1646.8618856663252, 1553.2220460691974, 1451.8744754266536, 1350.582859274457, 
            1251.9162470079896, 1051.5664666811736, 957.5075512657118, 850.7471608277119, 756.0325287284863, 
            656.6550224336061],
        forg=[0.08958593677142582, 0.10020889676396522, 0.1494628368926606, 0.19049706238925668, 0.17004010071808257, 
            0.11977338431409115, 0.12975286766763028, 0.14035064813951315, 0.14039261400727404, 0.14105567471789607, 
            0.16009238967121547, 0.1497514947102282, 0.19076697804304346, 0.2210044219590288, 0.21133867232428727, 
            0.14991501911778413]
    )

    # Plot 
    t = @. !isnan(frog_curve);
    h = plot(
        xlabel="Age [Ma.]", ylabel="f (org)",
        framestyle=:box,
        fontfamily=:Helvetica,
        xlims=(500,3000),
        xflip=true,
        ylims=(0,0.3),
        size=(400,500),
        legend=:topleft
    )
    plot!(des_marais.age, des_marais.forg, label="Des Marais",
        markershape=:circle, linestyle=:dash,
        color=:blue, msc=:auto,
    )
    plot!(c_carb[t], frog_curve[t], label="Replicate",
        markershape=:circle,
        color=:black, msc=:auto,
    )
    display(h)
    

## --- End of file 