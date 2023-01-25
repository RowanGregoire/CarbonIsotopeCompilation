#= 
    Match H/C ratios and δ13C org values from the same formations and uses
        this plot to fit a Rayleigh fractionation curve to the H/C and δ13C org relationship

    Resample δ13C and H/C values

    Correct δ13C org values for alteration using the Rayleigh fractionation curve

    Calculate fraction of carbon buried as organic matter
=#

## -- Load the packages, functions, and data we'll be using
using Plots
using StatGeochem
using Statistics
using Chron
using LsqFit
using StatsPlots

include("intake.jl")

isotopes = importdataset("data/isotopes.csv", ',', importas=:Tuple)
kerogen = importdataset("data/kerogen.csv", ',', importas=:Tuple)
#=
    Note that isotopes has samples calibrated to both VPDB and PDB - this should be fixed ASAP
=#

## -- VPDB to PDB conversion
# also adding uncertanties for age and isotopes if they dont have them?
# or is it better to only use the ones that have uncertainties?
n = length(isotopes.key)
isopdb = (key = isotopes.key, d13c_org = zeros(n), d13c_carb = zeros(n))

for i in 1:lastindex(isotopes.key)
    if isotopes.c_standard[i] == "PDB"
        #= 
            Based on the sync_vectors function I can probably modify the arrays
            that are in the isotopes tuple... and there's no reason to save the old 
            VPDB standard data here
        =#
    else
        convert_vpdb!()
    end
end


## -- Plot data and calibrated fractionation curve
matches = build_match(kerogen, isotopes)
plotmatch = plot(matches.hc, matches.d13c_org, seriestype=:scatter, msc=:auto,
    framestyle=:box, label="", xlabel="H/C ratio", ylabel="δC13 org", 
    title="H/C and δ13C Relationship"
)
savefig(plotmatch, "output/org_kgn_matches.png")

# hopefully this will work when I get the domain error to work?
(x, y) = fit_rayleigh(matches.d13c_org, matches.hc)

plot!(plotmatch, x, y, label="Rayleigh model")




## -- Visualize non-resampled data
isotoperaw = plot(isotopes.age, isotopes.d13c_org, label="organic",
    seriestype=:scatter, msc=:auto, alpha=0.25, framestyle=:box
)
plot!(isotopes.age, isotopes.d13c_carb, label="carbonate", seriestype=:scatter,
    msc=:auto, alpha=0.25, framestyle=:box, xlabel="Age [Ma]", ylabel="d13C [PDB OR VPDB]",
    title="Raw Isotope Data"
)
savefig(isotoperaw, "output/raw_data.png")

## -- Only use the data from isotopes and kerogen that has age and isotope uncertainties
naniso_age = NamedTuple{(:age, :uncert)}(sync_vectors(isotopes.age, isotopes.age_uncertainty))
nanorganic = NamedTuple{(:iso, :uncert)}(sync_vectors(isotopes.d13c_org, isotopes.corg_uncertainty))
nancarbonate = NamedTuple{(:iso, :uncert)}(sync_vectors(isotopes.d13c_carb, isotopes.ccarb_uncertainty))
nankgn_age = NamedTuple{(:age, :uncert)}(sync_vectors(kerogen.age, kerogen.age_uncertainty))


## -- Visualize the data loss
# Count the data in the filtered and unfiltered datasets and put the counts into matrices
edges = 0:100:3800
(countall_org, all_org_ctr) = count_data(isotopes.age, isotopes.d13c_org, edges)
(countorg, org_ctr) = count_data(isotopes.age, nanorganic.iso, edges)
knitorg = knit_vectors([countall_org, countorg])

(countall_carb, all_carb_ctr) = count_data(isotopes.age, isotopes.d13c_carb, edges)
(countcarb, carb_ctr) = count_data(isotopes.age, nancarbonate.iso, edges)
knitcarb = knit_vectors([countall_carb, countcarb])

(countall_hc, all_hc_ctr) = count_data(kerogen.age, kerogen.hc, edges)
(counthc, hc_ctr) = count_data(nankgn_age.age, nankgn_age.age, edges)
knithc = knit_vectors([countall_hc, counthc])

# Plot
gbo = groupedbar(knitorg, bar_width=0.7, xlabel="100 Ma bin", ylabel="count",
    title="δ13C organic", xlims=(0,38), label="", framestyle=:box
)
gbc = groupedbar(knitcarb, bar_width=0.7, xlabel="100 Ma bin", ylabel="count",
    title="δ13C carbonate", xlims=(0,38), label="", framestyle=:box
)
gbh = groupedbar(knithc, bar_width=0.7, xlabel="100 Ma bin", ylabel="count",
    title="H/C ratio (ages without uncertainties removed)", xlims=(0,38), label="",
    framestyle=:box
)
savefig(gbo, "output/loss_organic.png")
savefig(gbc, "output/loss_carbonate.png")
savefig(gbh, "output/loss_hc.png")


## -- Resample d13C and H/C values using only age and isotope values that have uncertainty
org_rs = NamedTuple{(:c, :m, :e)}(bin_bsr_means(naniso_age.age, nanorganic.iso, 0, 3800, 38,
    x_sigma=naniso_age.uncert, y_sigma=nanorganic.uncert, sem=:sigma
))

carb_rs = NamedTuple{(:c, :m, :e)}(bin_bsr_means(naniso_age.age, nancarbonate.iso, 0, 3800, 38,
    x_sigma=naniso_age.uncert, y_sigma=nancarbonate.uncert, sem=:sigma
))

hc_rs = NamedTuple{(:c, :m, :e)}(bin_bsr_means(nankgn_age.age, kerogen.hc, 0, 3800, 38,
    x_sigma=nankgn_age.uncert, y_sigma=kerogen.hc*0.01, sem=:sigma
))


## -- Use the Rayleigh fractionation curve to correct d13C values
#= 
    Parameters in as 0 for now because the fractionation function doesn't work :(

    It's so odd that this looks like it corrects the values evenly across the board??
    with the WILD exception of the last two values WHAT is going on??

    Maybe because the last two bins have such low values? There's not a lot of data there...
    It could also be a product of the equation used - switching to the Rayleigh equation
    could help 
=#
org_initial = org_rs.m .- Δδ.(hc_rs.m, 0)


## -- Plot resampled values with raw data for comparison
hydrocarbons = plot(kerogen.age, kerogen.hc, label="H/C ratio (raw data)", seriestype=:scatter,
    msc=:auto, alpha=0.25, framestyle=:box, xlabel="Age [Ma]", ylabel="H/C ratio"
)
plot!(hydrocarbons, hc_rs.c, hc_rs.m, label="resampled mean", yerror=hc_rs.e, 
    seriestype=:scatter, msc=:auto, title="Resampled H/C values"
)
#display(hydrocarbons)
savefig(hydrocarbons, "output/rs_kerogen.png")

organics = plot(isotopes.age, isotopes.d13c_org, label="organic (raw data)", seriestype=:scatter,
    msc=:auto, alpha=0.25, framestyle=:box, xlabel="Age [Ma]", ylabel="d13C [‰]"
)
plot!(organics, org_rs.c, org_rs.m, label="resampled mean", yerror=org_rs.e, 
    seriestype=:scatter, msc=:auto
)

#=
plot!(organics, org_rs.c, org_initial, label="corrected resampled mean", yerror=org_rs.e, 
    seriestype=:scatter, msc=:auto, legend=:bottomleft, title="Resampled δ13C organic"
)
=#
display(organics)
savefig(organics, "output/rs_organics.png")

carbonates = plot(isotopes.age, isotopes.d13c_carb, label="carbonate", seriestype=:scatter,
    msc=:auto, alpha=0.25, framestyle=:box, xlabel="Age [Ma]", ylabel="d13C [‰]"
)
plot!(carbonates, carb_rs.c, carb_rs.m, label="resampled mean", yerror=carb_rs.e, 
    seriestype=:scatter, msc=:auto, title="Resampled δ13C carbonate"
)
savefig(carbonates, "output/rs_carbonates.png")


## -- Calculate forg with resampled means
d13c_mantle = -5.5

# forg = (δmantle - δcarb) / (δorg - δcarb)
forg = @. (-1 * carb_rs.m + d13c_mantle) / (org_initial - carb_rs.m)

# Plot calculated forg and compare to the Des Marais plot and GOE
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

plotforg = plot(org_rs.c, forg, marker=:circle, msc=:auto, label="calculated f(org)")
plot!(plotforg, des_marais.age, des_marais.forg, marker=:circle, msc=:auto, 
    label="Des Marais f(org)", xlabel="Age [Ma]", ylabel="Fraction of carbon buried as organic", 
    title="Preliminary f(org)"   
)
plot!(plotforg, [2500], seriestype=:vline, label="GOE [2500 Ma]")
savefig(plotforg, "output/prelim_forg")


## -- Calculate δ13C MCMC minimum for each bin
# I feel like I'm not a huge fan of doing this for not-resampled data...
nsteps = 10000
dist = ones(10)
binedges = 0:100:3800
bincntrs = cntr(binedges)

mcmc_d13c = (iso = similar(bincntrs), uncert = similar(bincntrs))

for i in 1:lastindex(bincntrs)
    # Only look at the data that is in the bin and is not NaN
    # This is where it's actually important to only use the data that has uncertainties! Or something
    t = @. (binedges[i] < naniso_age.age < binedges[i+1]) && !isnan(naniso_age.age) && !isnan(nanorganic.iso)

    # Now that we've removed most of the data, only do the metropolis calculation if there are more than 2 data points
    if count(t) > 2
        d13c_t = nanorganic.iso[t]
        sig_t = nanorganic.uncert[t]

        mindist = metropolis_min(nsteps, dist, d13c_t, sig_t, burnin=5000)
        mcmc_d13c.iso[i] = nanmean(mindist)
        mcmc_d13c.uncert[i] = nanstd(mindist)
    else
        mcmc_d13c.iso[i] = mcmc_d13c.uncert[i] = NaN
    end

end

## Replot the organics plot - it's far enough up in the code that it makes sense to redo it
organics = plot(isotopes.age, isotopes.d13c_org, label="organic (raw data)", seriestype=:scatter,
    msc=:auto, alpha=0.25, framestyle=:box, xlabel="Age [Ma]", ylabel="d13C [‰]"
)
plot!(organics, org_rs.c, org_rs.m, label="resampled mean", yerror=org_rs.e, 
    seriestype=:scatter, msc=:auto
)
plot!(organics, bincntrs, mcmc_d13c.iso, yerror=mcmc_d13c.uncert, label="MCMC minimum",
    seriestype=:scatter, msc=:auto, legend=:bottomleft, size=(600,1200)
)
display(organics)
savefig(organics, "output/mcmc_d13c.png")
