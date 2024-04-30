## --- Set up 
    # Attempt to H/C correct all samples to give ourselves more data  

    # What if we average the d13C values for each formation? And have a standard 
    # deviation... that would maybe reduce the number of CIEs that show up

    # Packages
    using StatGeochem, Plots, Measurements, LsqFit
    using ProgressMeter
    using Isoplot: yorkfit
    using LoopVectorization: @turbo
    using LogExpFunctions: logsumexp
    using KernelDensity: kde

    include("utilities/Utilities.jl")

    # Load data 
    data = importdataset("data/compilation.csv", ',', importas=:Tuple)


## --- Create a resampled dataset with an H/C value for each organic d13c value 
    # Preallocate
    nsims = Int(1e6)
    simout = (;
        carb = Array{Float64}(undef, nsims, 2),    # d13c, age 
        org = Array{Float64}(undef, nsims, 5),     # d13c, H/C, age, lat, lon
    )

    # Isotope data without uncertainty is assigned an uncertainty of 0.02 ‰
    data.d13c_carb_uncert[isnan.(data.d13c_carb_uncert)] .= 0.02
    data.d13c_org_uncert[isnan.(data.d13c_org_uncert)] .= 0.02

    # Age data is assigned a minimum uncertainty of 10%
    for i in eachindex(data.age_uncert)
        data.age_uncert[i] = nanmax(data.age_uncert[i], data.age[i]*0.05)
    end

    # Resample carbonate isotope data
    t = .!isnan.(data.d13c_carb)
    k = invweight(data.lat[t], data.lon[t], data.age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    simout.carb .= bsresample([data.d13c_carb[t] data.age[t]],
        [data.d13c_carb_uncert[t] data.age_uncert[t]],
        nsims,p
    )

    # Organic and H/C 
    t = .!isnan.(data.d13c_org)
    k = invweight(data.lat[t], data.lon[t], data.age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    simout.org .= bsresample([data.d13c_org[t] data.hc[t] data.age[t] data.lat[t] data.lon[t]],
        [data.d13c_org_uncert[t] fill(0.01, count(t)) data.age_uncert[t] data.loc_uncert[t] data.loc_uncert[t]],
        nsims,p
    )


## --- Try again!! But this time use an age - H/C value model 
    # Ages without uncertainty assigned an uncertainty of 5%
    ageuncert = Array{Float64}(undef, length(data.age), 1)
    for i in eachindex(ageuncert)
        ageuncert[i] = ifelse(isnan(data.age_uncert[i]), data.age[i]*0.05, data.age_uncert[i])
    end

    # Fit model, forcing a y intercept of 1.5
    t = @. !isnan(data.d13c_org) & !isnan(data.hc) & (data.hc > 0);
    t .&= data.std_fm_name .!= "Onverwacht Gp";

    x = [fill(0, 500); data.age[t]]
    y = [log.(fill(1.5, 500)); log.(data.hc[t])]

    x_sigma = [fill(1e-8, 500); ageuncert[t]]
    y_sigma = [fill(1e-8, 500); log.(fill(0.1, count(t)))]

    fobj = yorkfit(x, x_sigma, y, y_sigma)
    println("$(fobj.intercept)")
    hc_age(age) = exp(age * (fobj.slope) + (fobj.intercept))

    x = 50:50:3800
    y = Measurements.value.(hc_age.(x))
    e = Measurements.uncertainty.(hc_age.(x))
    h = plot(data.age[t], data.hc[t],
        label="Observed", 
        seriestype=:scatter, markersize=3,
        color=:darkturquoise, msc=:auto,
        xlabel="Age [Ma.]", ylabel="[LOG] H/C Ratio",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        yaxis=:log10
    )
    plot!(x, y, ribbon=2*e, 
        label="Modeled ± 2 s.d.", 
        linewidth=2, color=:teal
    )
    
    display(h)
    savefig(h, "figures/hc_age.pdf")

    # For any organic carbon value without a H/C ratio: 
        # Pick age randomly from a Gaussian distribution: mean = age, std = age uncert
        # Pick H/C randomly from a Gaussian distribution: mean = modeled HC, std = uncert 
    hc_assigned = Array{Float64}(undef, length(data.hc), 1)
    for i in eachindex(data.hc)
        hc = hc_age(randn() * ageuncert[i] + data.age[i])
        hc_assigned[i] = ifelse(!isnan(data.hc[i]), data.hc[i], randn() * hc.err + hc.val)
    end

    # Sanity check: re-plot data for Rayleigh curve 
    t = @. !isnan(data.d13c_org) 
    h = plot(hc_assigned[t], data.d13c_org[t], 
        seriestype=:scatter, label="Assigned", 
        color=:lightblue, msc=:auto,
        markersize=2,
        framestyle=:box, 
        ylabel="d13c organic", xlabel="H/C ratio"
    )
    plot!(data.hc[t], data.d13c_org[t], 
        seriestype=:scatter, label="Observed", 
        color=:red, msc=:auto,
        markersize=2
    )

    t = @. !isnan(data.d13c_org) & !isnan(data.hc)
    params = fit_rayleigh(data.d13c_org[t], data.hc[t])
    x,y = rayleigh_curve(params, data.hc[t])
    plot!(x,y, label="Rayleigh Model", color=:black, linewidth=2)
    
    display(h)
    savefig(h, "figures/rayleigh_agemodel.pdf")


## --- Correct organic carbon data for post-depositional alteration with Rayleigh curve
    # Model with existing H/C data and correct with assigned H/C data
    t = @. !isnan(data.d13c_org) & !isnan(data.hc)
    params = fit_rayleigh(data.d13c_org[t], data.hc[t])
    d13c_org_corr = data.d13c_org .- (vec(r₀(hc_assigned, params)) .-  r₀(1.5, params))
    
    # Plot data 
    t = @. !isnan(data.d13c_org);
    h = plot(
        framestyle=:box,
        xlabel="Age [Ma.]", ylabel="d13c organic",
    )

    c,m,e = binmeans(data.age[t], data.d13c_org[t], 0,3800,38)
    s = .!isnan.(m)
    plot!(c[s], m[s], yerror=2e, label="Observed", 
        color=:darkorange, lcolor=:darkorange, msc=:auto, 
        # seriestype=:scatter,
        markershape=:circle,
    )

    c,m,e = binmeans(data.age[t], d13c_org_corr[t], 0,3800,38)
    s = .!isnan.(m)
    plot!(c[s], m[s], yerror=2e, label="Corrected", 
        color=:seagreen, lcolor=:seagreen, msc=:auto, 
        # seriestype=:scatter,
        markershape=:circle,
    )

    display(h)
    savefig(h, "figures/isotope_corrected_v2.pdf")


## --- Fraction buried as organic? Something is horribly and terribly wrong here 
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
        # ylims=(0,0.3),
        ylims=(0,0.4),
        size=(400,500),
        legend=:bottomleft,
        fg_color_legend=:white
    );

    # Uncorrected data
    t = .!isnan.(data.d13c_org);
    c,m,e = binmeans(data.age[t], data.d13c_org[t], 0, 3800, 38, relbinwidth=3)
    frog = (mantle .- carbonate) ./ (m .- carbonate)
    s = .!isnan.(frog)
    plot!(c[s], frog[s], label="Uncorrected n = $(count(t))",
        color=:darkorange, msc=:auto,
        markershape=:circle,
    )

    # Corrected data given estimated H/C ratios
    c,m,e = binmeans(data.age[t], d13c_org_corr[t], 0, 3800, 38, relbinwidth=3)
    frog = (mantle .- carbonate) ./ (m .- carbonate)
    s = .!isnan.(frog)
    plot!(c[s], frog[s], label="Corrected (est. H/C) n = $(count(t))",
        color=:seagreen, msc=:auto,
        markershape=:circle,
    )

    # Corrected data, but only those with recorded H/C analyses
    t = .!isnan.(data.d13c_org) .& .!isnan.(data.hc);
    c,m,e = binmeans(data.age[t], d13c_org_corr[t], 0, 3800, 38, relbinwidth=3)
    frog = (mantle .- carbonate) ./ (m .- carbonate)
    s = .!isnan.(frog)
    plot!(c[s], frog[s], label="Corrected (obs. H/C) n = $(count(t))",
        color=:royalblue, msc=:auto,
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
        color=:black, msc=:auto,
    )

    display(h)
    savefig(h, "figures/fraction_organic_v2.pdf")


## --- End of file