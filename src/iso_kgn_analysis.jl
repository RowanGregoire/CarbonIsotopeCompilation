## -- Load the packages, functions, and data we'll be using
using Plots
using StatGeochem
using Statistics
using Chron
using LsqFit
using StatsPlots
using LoopVectorization
using ProgressMeter

include("intake.jl")

data = importdataset("data/compilation.csv", ',', importas=:Tuple)

## -- Plot d13C org and H/C relationship and calibrated fractionation curve
t = @. !isnan(data.d13c_org) && !isnan(data.hc)
plotmatch = plot(data.hc[t], data.d13c_org[t], seriestype=:scatter, msc=:auto,
    framestyle=:box, label="", xlabel="H/C ratio", ylabel="δC13 org", 
    title="H/C and δ13C Relationship"
)

params = fit_rayleigh(data.d13c_org[t], data.hc[t])
x,y = rayleigh_curve(params, data.hc[t])
plot!(plotmatch, x, y, label="Rayleigh model")

# don't know why this is erroring here but the x and y vectors aren't the same length
savefig(plotmatch, "output/org_kgn_matches.pdf")


## -- Visualize non-resampled data
isotoperaw = plot(data.age, data.d13c_org, label="organic",
    seriestype=:scatter, msc=:auto, alpha=0.25, framestyle=:box
)
plot!(data.age, data.d13c_carb, label="carbonate", seriestype=:scatter,
    msc=:auto, alpha=0.25, framestyle=:box, xlabel="Age [Ma]", ylabel="d13C [PDB OR VPDB]",
    title="Raw Isotope Data"
)
savefig(isotoperaw, "output/raw_data.pdf")


## -- Estimate uncertainties for samples without known uncertainties
# Make a new data set so we can reference back to the old one if needed
newdata = deepcopy(data)

# Uncertainties
uage = 0.05     # multiply by 5%
uiso = 0.02     # ± 0.02‰

add_uncert!(newdata.age, newdata.age_uncertainty, err=uage, method=:percent, addnote=true, 
    flag=newdata.flag, comment=newdata.comments, note="est. 5% age uncertainty"
)
add_uncert!(newdata.d13c_carb, newdata.ccarb_uncertainty, err=uiso, method=:plusminus, addnote=true, 
    flag=newdata.flag, comment=newdata.comments, note="est. ± 0.02‰ δ13C carb uncertainty"
)
add_uncert!(newdata.d13c_org, newdata.corg_uncertainty, err=uiso, method=:plusminus, addnote=true, 
    flag=newdata.flag, comment=newdata.comments, note="est. ± 0.02‰ δ13C org uncertainty"
)


## -- Resample d13C and H/C values
org_rs = NamedTuple{(:c, :m, :e)}(bin_bsr_means(newdata.age, newdata.d13c_org, 0, 3800, 38,
    x_sigma=newdata.age_uncertainty, y_sigma=newdata.corg_uncertainty, sem=:sigma
))

carb_rs = NamedTuple{(:c, :m, :e)}(bin_bsr_means(newdata.age, newdata.d13c_carb, 0, 3800, 38,
    x_sigma=newdata.age_uncertainty, y_sigma=newdata.ccarb_uncertainty, sem=:sigma
))

hc_rs = NamedTuple{(:c, :m, :e)}(bin_bsr_means(newdata.age, newdata.hc, 0, 3800, 38,
    x_sigma=newdata.age_uncertainty, y_sigma=newdata.hc*0.01, sem=:sigma
))


## -- Bootstrap resample data in an array using a sample-per-row method
# data_head = ["lat" "long" "age" "d13c_org" "d13c_carb" "hc"]
# sigma_head = ["lat_sigma" "long_sigma" "age_sigma" "d13c_org_sigma" "d13c_carb_sigma" "hc_sigma"]

# Format data for bsresample()
fmtdata = unelementify(newdata, (:lat, :long, :age, :d13c_org, :d13c_carb, :hc), floatout=true)
fmtsigma = unelementify(newdata, (:coordinate_uncertainty, :coordinate_uncertainty, :age_uncertainty, :corg_uncertainty, :ccarb_uncertainty,), floatout=true)
fmtsigma = hcat(fmtsigma, data.hc*0.01)   # Add 1% uncertainty for H/C ratio

nrows = 100000  # Gets mad at 1e6
k = invweight(newdata.lat, newdata.long, newdata.age)   # Calculate inverse weights
p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)             # Keep rougly one-fith of the data in each resampling
bsr_data = bsresample(fmtdata, fmtsigma, nrows, p)


## --- Sanitize data and get H/C ratios for all samples
t = @. bsr_data[:,3] > 0    # Only data where the age is greater than 0

sbsr = (lat = bsr_data[:,1][t], long = bsr_data[:,2][t], age = bsr_data[:,3][t], 
    org = bsr_data[:,4][t], carb = bsr_data[:,5][t], hc = bsr_data[:,6][t]
)

# Estimate H/C ratios for samples that don't have them
bound = 25  # ± 25 Ma range to grab an H/C ratio. Fairly arbitrary
p = Progress(length(sbsr.age) ÷ 10 , desc = "Finding H/C Ratios: ")

#=
    check - this would not be a good place to use @turbo becuase I need to do things before
    assigning a value to sbsr.hc
=#
for i in eachindex(sbsr.hc)
    if isnan(sbsr.hc[i])
        # Get all samples in the age range
        target_age = sbsr.age[i]
        t = @. (target_age - bound) < sbsr.age < (target_age + bound)

        # Make sure the bound has data
        s = @. !isnan(sbsr.hc[t])
        newbound = bound
        while isempty(sbsr.hc[t][s])
            # L + ratio + no data in your subset
            # Is there a way to make this recursive? seems to throw me into a infinite loop
            # Is that even useful?
            newbound += 5
            t = @. (target_age - newbound) < sbsr.age < (target_age + newbound)
            s = @. !isnan(sbsr.hc[t])
        end

        # Randomly pick one and assign it to be the H/C value for that sample
        sbsr.hc[i] = rand(sbsr.hc[t][s])
    end
    (i % 10 == 0) && next!(p)
end


## --- Compute H/C adjustment using Rayleigh and Des Mariais style corrections
org_initial_rf = sbsr.org .- Δδ(sbsr.hc, params)

# Exclude values that are really small because it makes the Des Marais correction sad 
t = @. sbsr.hc > 0.015  # Anything smaller than this hits some weird asymptote
org_initial_dm = sbsr.org[t] .- Δδ(sbsr.hc[t], :DesMarais)

# Bin data and find mean and error
org_bin = NamedTuple{(:c, :m, :e)}(bin_means(sbsr.age, sbsr.org, 0:100:3800))       # Uncorrected
carb_bin = NamedTuple{(:c, :m, :e)}(bin_means(sbsr.age, sbsr.carb, 0:100:3800))     # No relevant correction
dm_bin = NamedTuple{(:c, :m, :e)}(bin_means(sbsr.age[t], org_initial_dm, 0:100:3800))
rf_bin = NamedTuple{(:c, :m, :e)}(bin_means(sbsr.age, org_initial_rf, 0:100:3800))


## --- Plot results
# Set things up for the f_org calculation
# forg = (δmantle - δcarb) / (δorg - δcarb)
d13c_mantle = -5.5
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

# Compare uncorrected vs Des Marais corrected data and compute f_org
pdm = plot(org_bin.c, org_bin.m, yerror=2*org_bin.e, seriestype=:scatter, msc=:auto, 
    label="uncorrected", xlabel="Age [Ma]", ylabel="δ13C (org) [‰]", framestyle=:box, 
    legend=:bottomleft, ylims=(-70,10)
)
plot!(pdm, dm_bin.c, dm_bin.m, yerror=2*dm_bin.e, seriestype=:scatter, msc=:auto, 
    label="Des Marais correction", xlabel="Age [Ma]", ylabel="δ13C (org) [‰]", framestyle=:box
)
display(pdm)
savefig(pdm, "output/corrected_desmarais.pdf")

forg_dm = @. (d13c_mantle .- carb_bin.m) / (dm_bin.m - carb_bin.m)
plotforg = plot(rf_bin.c, forg, marker=:circle, msc=:auto, label="calculated f(org)")
plot!(plotforg, des_marais.age, des_marais.forg, marker=:circle, msc=:auto, 
    label="Des Marais f(org)", xlabel="Age [Ma]", ylabel="Fraction of carbon buried as organic", 
    title="Preliminary f(org) - Des Marais Correction", framestyle=:box   
)
display(plotforg)
savefig(plotforg, "output/forg_desmarais.pdf")

# Compare uncorrected vs. Rayleigh fractionation corrected data and compute f_org
prf = plot(org_bin.c, org_bin.m, yerror=2*org_bin.e, seriestype=:scatter, msc=:auto, 
    label="uncorrected", xlabel="Age [Ma]", ylabel="δ13C (org) [‰]", framestyle=:box
)
plot!(prf, rf_bin.c, rf_bin.m, yerror=2*rf_bin.e, seriestype=:scatter, msc=:auto, 
    label="Rayleigh correction", xlabel="Age [Ma]", ylabel="δ13C (org) [‰]", framestyle=:box, 
    legend=:bottomleft, ylims=(-70,10)
)
display(prf)
savefig(prf, "output/corrected_rayleigh.pdf")

forg_rf = @. (d13c_mantle .- carb_bin.m) / (rf_bin.m - carb_bin.m)
plotforg = plot(rf_bin.c, forg, marker=:circle, msc=:auto, label="calculated f(org)")
plot!(plotforg, des_marais.age, des_marais.forg, marker=:circle, msc=:auto, 
    label="Des Marais f(org)", xlabel="Age [Ma]", ylabel="Fraction of carbon buried as organic", 
    title="Preliminary f(org) - Rayleigh Correction", framestyle=:box   
)
display(plotforg)
savefig(plotforg, "output/forg_rayleigh.pdf")


#=
TO DO: This has not been corrected for the new resampled data because this data does not
propagate uncertainties. Unsure what to put the uncertainty needed for this function as


## --- Calculate δ13C MCMC minimum for each bin using the Rayleigh-corrected data
nsteps = 10000
dist = ones(10)
binedges = 0:100:3800
bincntrs = cntr(binedges)

mcmc_d13c = (iso = similar(bincntrs), uncert = similar(bincntrs))

for i in 1:lastindex(bincntrs)
    # Only look at the data that is in the bin and is not NaN
    t = @. (binedges[i] < sbsr.age < binedges[i+1]) && !isnan(sbsr.age) && !isnan(org_initial_rf)

    # Only do the metropolis calculation if there are more than 2 data points
    if count(t) > 2
        d13c_t = org_initial_rf[t]
        sig_t = newdata.corg_uncertainty[t]

        mindist = metropolis_min(nsteps, dist, d13c_t, sig_t, burnin=5000)
        mcmc_d13c.iso[i] = nanmean(mindist)
        mcmc_d13c.uncert[i] = nanstd(mindist)
    else
        mcmc_d13c.iso[i] = mcmc_d13c.uncert[i] = NaN
    end
end

## Replot the organics plot - it's far enough up in the code that it makes sense to redo it
organics = plot(newdata.age, newdata.d13c_org, label="organic (raw data)", seriestype=:scatter,
    msc=:auto, alpha=0.25, framestyle=:box, xlabel="Age [Ma]", ylabel="d13C [‰]"
)
plot!(organics, org_rs.c, org_rs.m, label="resampled mean", yerror=org_rs.e, 
    seriestype=:scatter, msc=:auto
)
plot!(organics, bincntrs, mcmc_d13c.iso, yerror=mcmc_d13c.uncert, label="MCMC minimum",
    seriestype=:scatter, msc=:auto, legend=:bottomleft, size=(600,1200)
)
display(organics)
savefig(organics, "output/mcmc_d13c.pdf")
=#