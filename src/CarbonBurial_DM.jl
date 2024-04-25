## --- Set up
    # Test organic carbon burial over geologic time using methods modified from 
    # Des Marais et al., 1992 (10.1038/359605a0)

    # Set up 
    using StatGeochem, Plots
    using LsqFit

    include("utilities/Utilities.jl")

    # Load data 
    data = importdataset("data/compilation.csv", ',', importas=:Tuple)

    
## --- [PLOT] Raw / uncorrected isotope data 
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        size=(400,600),
        xlabel="Age [Ma.]", ylabel="δ¹³C [‰]",
        xlims=(0,3800),
        # xflip=true,
        # ylims=(-60,10)
        fg_color_legend=:white,
        legend=:bottomleft,
    );

    # Raw data 
    plot!(data.age, data.d13c_carb, label="Carbonate",
        seriestype=:scatter,
        markersize=2, markershape=:+,
        color=colors.carb_light, msc=:auto,
    )
    plot!(data.age, data.d13c_org, label="Organic",
        seriestype=:scatter,
        markersize=2,
        color=colors.org_light, msc=:auto,
    )

    # Moving means of 100 Ma averaged over 200 Ma moving window
    c,m,e = binmeans(data.age, data.d13c_carb, 0, 3800, 38, relbinwidth=3)
    plot!(c, m, label="", color=colors.carb_dark, linewidth=2)

    c,m,e = binmeans(data.age, data.d13c_org, 0, 3800, 38, relbinwidth=3)
    plot!(c, m, label="", color=colors.org_dark, linewidth=2,)

    # Show plot 
    display(h)
    savefig(h, "figures/isotope_raw.pdf")


## --- Correct organic carbon for post-depositional alteration with Des Marais curve 
    Δδ(hc::Number) = 4.05 - 3.05*hc + 0.785/hc + 0.0165/hc^2 - (8.79E-4)/hc^3

    t = @. !isnan(data.d13c_org) & !isnan(data.hc) .& (data.hc >= 0.1)
    d13c_org_corrected = data.d13c_org[t] .- Δδ.(data.hc[t])

    # Plot data 
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        size=(400,600),
        xlabel="Age [Ma.]", ylabel="δ¹³C [‰]",
        xlims=(0,3800),
        # xflip=true,
        # ylims=(-60,10)
        fg_color_legend=:white,
        legend=:bottomleft,
    );

    # Carbonate
    plot!(data.age, data.d13c_carb, label="Carbonate",
        seriestype=:scatter,
        markersize=2, markershape=:+,
        color=colors.carb_light, msc=:auto,
    )
    c,m,e = binmeans(data.age, data.d13c_carb, 0, 3800, 38, relbinwidth=3)
    plot!(c, m, label="", color=colors.carb_dark, linewidth=2)

    # Uncorrected organic 
    plot!(data.age[t], data.d13c_org[t], label="Uncorrected Organic",
        seriestype=:scatter,
        markersize=2,
        color=colors.org_light, msc=:auto,
    )
    c,m,e = binmeans(data.age[t], data.d13c_org[t], 0, 3800, 38, relbinwidth=3)
    s = .!isnan.(m)
    plot!(c[s], m[s], label="",linewidth=1, 
        color=colors.org_dark, msc=colors.org_dark
    )
    plot!(c[s], m[s], yerror=e[s], label="",
        linewidth=2,
        seriestype=:scatter,
        markershape=:circle, markersize=3,
        color=colors.org_dark, msc=colors.org_dark
    )

    # Corrected organic
    plot!(data.age[t], d13c_org_corrected, label="Corrected Organic",
        seriestype=:scatter,
        markersize=2,
        color=:palegreen, msc=:auto,
    )
    c,m,e = binmeans(data.age[t], d13c_org_corrected, 0, 3800, 38, relbinwidth=3)
    s = .!isnan.(m)
    plot!(c[s], m[s], label="",linewidth=1, 
        color=:forestgreen, msc=:forestgreen
    )
    plot!(c[s], m[s], yerror=e[s], label="",
        linewidth=2,
        seriestype=:scatter,
        markershape=:circle, markersize=3,
        color=:forestgreen, msc=:forestgreen
    )

    display(h)
    savefig(h, "figures/isotope_corrected.pdf")


## Calculate fraction of carbon buried as organic 
    # Use running mean of carbonate values to smooth CIEs
    mantle = -5.5
    c,carbonate,e = binmeans(data.age, data.d13c_carb, 0, 3800, 38, relbinwidth=3)

    # Plot base
    h = plot(
        xlabel="Age [Ma.]", ylabel="Fraction Buried as Organic",
        framestyle=:box,
        fontfamily=:Helvetica,
        # xlims=(500,3000),
        # xflip=true,
        ylims=(0,0.3),
        size=(400,500),
        legend=:bottomleft,
        fg_color_legend=:white
    );

    # All data without correction 
    c,m,e = binmeans(data.age, data.d13c_org, 0, 3800, 38, relbinwidth=3)
    frog = (mantle .- carbonate) ./ (m .- carbonate)
    s = .!isnan.(frog)
    plot!(c[s], frog[s], label="All data n = $(count(.!isnan.(data.d13c_org)))",
        color=:springgreen, msc=:auto,
        markershape=:circle,
    )

    # # Uncorrected data with known H/C ratios
    # c,m,e = binmeans(data.age[t], data.d13c_org[t], 0, 3800, 38, relbinwidth=3)
    # frog = (mantle .- carbonate) ./ (m .- carbonate)
    # s = .!isnan.(frog)
    # plot!(c[s], frog[s], label="Uncorrected n = $(count(t))",
    #     color=:navy, msc=:auto,
    #     markershape=:circle,
    # )

    # Corrected data given known H/C ratios
    c,m,e = binmeans(data.age[t], d13c_org_corrected, 0, 3800, 38, relbinwidth=3)
    frog = (mantle .- carbonate) ./ (m .- carbonate)
    s = .!isnan.(frog)
    plot!(c[s], frog[s], label="Corrected n = $(count(t))",
        color=:darkgreen, msc=:auto,
        markershape=:circle,
    )

    # Des Marais curve 
    des_marais = (;
        age=[2650.0, 2495.5516180173463, 2047.2862072181294, 1949.6316329385436, 
            1853.1940688240234, 1747.3141844633037, 1646.8618856663252, 1553.2220460691974, 
            1451.8744754266536, 1350.582859274457, 1251.9162470079896, 1051.5664666811736, 
            957.5075512657118, 850.7471608277119, 756.0325287284863, 656.6550224336061],
        forg=[0.08958593677142582, 0.10020889676396522, 0.1494628368926606, 0.19049706238925668, 
            0.17004010071808257, 0.11977338431409115, 0.12975286766763028, 0.14035064813951315, 
            0.14039261400727404, 0.14105567471789607, 0.16009238967121547, 0.1497514947102282, 
            0.19076697804304346, 0.2210044219590288, 0.21133867232428727, 0.14991501911778413]
    )
    plot!(des_marais.age, des_marais.forg, label="Des Marais",
        markershape=:circle, linestyle=:dash,
        color=:darkorange, msc=:auto,
    )

    display(h)
    savefig(h, "figures/fraction_organic.pdf")

## --- End of file