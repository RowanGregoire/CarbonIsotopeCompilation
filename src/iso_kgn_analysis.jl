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
    framestyle=:box, label="", xlabel="H/C ratio", ylabel="δC13 org" 
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
    msc=:auto, alpha=0.25, framestyle=:box, xlabel="Age [Ma]", ylabel="d13C [PDB OR VPDB]"
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
# Format data for bsresample()
ndata = length(newdata.key)
fmtdata = unelementify(newdata, (:lat, :long, :age, :d13c_org, :corg_uncertainty, :d13c_carb, :ccarb_uncertainty, :hc), floatout=true)
fmtsigma = unelementify(newdata, (:coordinate_uncertainty, :coordinate_uncertainty, :age_uncertainty, :corg_uncertainty), floatout=true)

# Add 1% uncertainty for H/C ratio and 0 uncertainty for resampled uncertainties
fmtsigma = hcat(fmtsigma, zeros(ndata), newdata.ccarb_uncertainty, zeros(ndata), data.hc*0.01)

# Resample
nrows = 100000  # Gets mad at 1e6
k = invweight(newdata.lat, newdata.long, newdata.age)   # Calculate inverse weights
p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)             # Keep rougly one-fith of the data in each resampling
bsr_data = bsresample(fmtdata, fmtsigma, nrows, p)


## --- Sanitize data and get H/C ratios for all samples
t = @. bsr_data[:,3] > 0    # Only data where the age is greater than 0
sbsr = (lat = bsr_data[:,1][t], long = bsr_data[:,2][t], age = bsr_data[:,3][t], 
    org = bsr_data[:,4][t], org_sigma = bsr_data[:,5][t], 
    carb = bsr_data[:,6][t], carb_sigma = bsr_data[:,7][t],
    hc = bsr_data[:,8][t]
)

estimate_hc!(sbsr.age, sbsr.hc, 5)


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
savefig(pdm, "output/corrected_desmarais.pdf")

forg_dm = @. (d13c_mantle .- carb_bin.m) / (dm_bin.m - carb_bin.m)
plotforg = plot(rf_bin.c, forg_dm, marker=:circle, msc=:auto, label="calculated f(org)")
plot!(plotforg, des_marais.age, des_marais.forg, marker=:circle, msc=:auto, 
    label="Des Marais f(org)", xlabel="Age [Ma]", ylabel="Fraction of carbon buried as organic", 
     framestyle=:box, legend=:topright    
)
savefig(plotforg, "output/forg_desmarais.pdf")

# Compare uncorrected vs. Rayleigh fractionation corrected data and compute f_org
prf = plot(org_bin.c, org_bin.m, yerror=2*org_bin.e, seriestype=:scatter, msc=:auto, 
    label="uncorrected", xlabel="Age [Ma]", ylabel="δ13C (org) [‰]", framestyle=:box
)
plot!(prf, rf_bin.c, rf_bin.m, yerror=2*rf_bin.e, seriestype=:scatter, msc=:auto, 
    label="Rayleigh correction", xlabel="Age [Ma]", ylabel="δ13C (org) [‰]", framestyle=:box, 
    legend=:bottomleft, ylims=(-70,10)
)
savefig(prf, "output/corrected_rayleigh.pdf")

forg_rf = @. (d13c_mantle .- carb_bin.m) / (rf_bin.m - carb_bin.m)
plotforg = plot(rf_bin.c, forg_rf, marker=:circle, msc=:auto, label="Rayleigh Corrected f(org)")
plot!(plotforg, des_marais.age, des_marais.forg, marker=:circle, msc=:auto, 
    label="Des Marais f(org)", xlabel="Age [Ma]", ylabel="Fraction of carbon buried as organic", 
    framestyle=:box, legend=:topright   
)
savefig(plotforg, "output/forg_rayleigh.pdf")


## --- Using the Rayleigh corrected data, make plots of the raw and resampled data 
# Carbonate 
plot(data.age, data.d13c_carb, seriestype=:scatter, label="raw data", alpha = 0.25, msc=:auto)
plot!(carb_rs.c, carb_rs.m, yerror=carb_rs.e, seriestype=:scatter, label="resampled mean", msc=:auto,
    framestyle=:box, legend=:topright, xlabel="Age [Ma]", ylabel="δ13C [‰]"
)
savefig("output/rs_carbonates.pdf")

# Organic 
plot(data.age, data.d13c_org, seriestype=:scatter, label="raw data", alpha = 0.25, msc=:auto)
plot!(org_rs.c, org_rs.m, yerror=org_rs.e, seriestype=:scatter, label="resampled mean", msc=:auto,
    framestyle=:box, legend=:topright, xlabel="Age [Ma]", ylabel="δ13C [‰]"
)
plot!(rf_bin.c, rf_bin.m, yerror=org_bin.e, seriestype=:scatter, label="Rayleigh-corrected", msc=:auto)
savefig("output/rs_organics.pdf")

# H/C
plot(data.age, data.hc, seriestype=:scatter, label="raw data", alpha = 0.25, msc=:auto)
plot!(hc_rs.c, hc_rs.m, yerror=hc_rs.e, seriestype=:scatter, label="resampled mean", msc=:auto,
    framestyle=:box, legend=:topright, xlabel="Age [Ma]", ylabel="H/C"
)
savefig("output/rs_kerogen.pdf")


## -- Plot the Phanerozoic in greater resolution to see what's going on
t = @. newdata.age <= 541
phan = plot(newdata.age[t], newdata.d13c_org[t], seriestype=:scatter, label="organic", msc=:auto, alpha = 0.25)
plot!(phan, newdata.age[t], newdata.d13c_carb[t], seriestype=:scatter, label="carbonate", msc=:auto, alpha = 0.25, 
    framestyle=:box, legend=:topleft, xlabel="Age [Ma]", ylabel="δ13C [‰]", ylims=(-40,20), size=(600,600)
)
savefig("output/phanerozoic.pdf")


## -- Calculate MCMC minimum for the new resampled dataset
# Preallocate
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
        sig_t = sbsr.carb_sigma[t]

        mindist = metropolis_min(nsteps, dist, d13c_t, sig_t, burnin=5000)
        mcmc_d13c.iso[i] = nanmean(mindist)
        mcmc_d13c.uncert[i] = nanstd(mindist)
    else
        mcmc_d13c.iso[i] = mcmc_d13c.uncert[i] = NaN
    end
end

plot(data.age, data.d13c_org, seriestype=:scatter, label="raw data", alpha = 0.25, msc=:auto)
plot!(rf_bin.c, rf_bin.m, yerror=org_bin.e, seriestype=:scatter, label="Rayleigh-corrected", msc=:auto, 
    framestyle=:box, legend=:topright, xlabel="age", ylabel="δ13C [‰]"
)
plot!(bincntrs, mcmc_d13c.iso, yerror=mcmc_d13c.uncert, label="MCMC minimum", msc=:auto, seriestype=:scatter)
savefig("output/mcmc_d13c.pdf")



## -- Try replicating the Des Marais plots 
# Restrict the data to the Schopf and Klein 1992 data
dm_sources = [
    "Abell et al. 1985 \"Archean\""
    "Schidlowski et al. 1979"
    "Schopf and Klein 1992"
    "Thode and Goodwin 1983"
    "Veizer and Hoefs 1976"
    "McKirdy and Powell 1974"
]

# Since we'll be using it frequently, make a new data set that's just the data we want
t = falses(length(newdata.key))      # Initialize BitVector
for i in eachindex(newdata.data_reference)
    for j in dm_sources
        if newdata.data_reference[i] == j
            t[i] = true
            break
        end
    end
end

dm_data = (
    key = newdata.key[t], 
    lat = newdata.lat[t],
    long = newdata.long[t],
    coordinate_uncertainty = newdata.coordinate_uncertainty[t],
    age = newdata.age[t],
    age_uncertainty = newdata.age_uncertainty[t],
    d13c_org = newdata.d13c_org[t], 
    corg_uncertainty = newdata.corg_uncertainty[t],
    d13c_carb = newdata.d13c_carb[t],
    ccarb_uncertainty = newdata.ccarb_uncertainty[t],
    hc = newdata.hc[t]
)

# Apply Des Marais correction and calculate forg
δ₀ = dm_data.d13c_org .- Δδ(dm_data.hc, :DesMarais)
obin = NamedTuple{(:c, :m, :e)}(bin_means(dm_data.age, dm_data.d13c_org, 0:100:3800))
cbin = NamedTuple{(:c, :m, :e)}(bin_means(dm_data.age, dm_data.d13c_carb, 0:100:3800))
forg = @. (d13c_mantle .- cbin.m) / (obin.m - cbin.m)

plotforg = plot(obin.c, forg, marker=:circle, msc=:auto, label="calculated f(org)")
plot!(plotforg, des_marais.age, des_marais.forg, marker=:circle, msc=:auto, 
    label="Des Marais f(org)", xlabel="Age [Ma]", ylabel="Fraction of carbon buried as organic", 
    title="Des Marias f(org) Replicate", framestyle=:box, legend=:topright   
)
display(plotforg)
# hmmm 🤔😟😨😞😢