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
            This won't work because you've already initialized the tuple
            The algorithm needs to initialize the tuple at the same time it puts 
            data into it

            But based on the sync_vectors function I can probably modify the arrays
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

#= hopefully this will work when I get the domain error to work?
(x, y) = fit_rayleigh(matches.d13c_org, matches.hc)

plot!(plotmatch, x, y, label="Rayleigh model")

=#


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
display(hydrocarbons)
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
#display(organics)
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
forg = @. (-1 * carb_rs.m + d13c_mantle) / (org_intial - carb_rs.m)

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



#=
## -- Find minimum using metropolis and subtract from each data point
# base metropolis code is from calculate_forg.jl, added HC calculation here
nsteps = 10000
dist = ones(10)

binedges = 0:100:3800
bincenters = cntr(binedges)
d13C_min = similar(bincenters)          #δ13C
d13C_min_sigma = similar(bincenters)
HC_min = similar(bincenters)            # H/C
HC_min_sigma = similar(bincenters)

# Iterate through each bin center
for i = 1:length(bincenters)
    # Get indexes for data that falls in the bin and is not NaN
    local t = @. (binedges[i] < isotopes.Age < binedges[i+1]) & !isnan(isotopes.d13C_org)
    local s = @. (binedges[i] < kerogen.Age < binedges[i+1]) & !isnan(kerogen.HC_Ratio)

    # δ13C: only do metropolis calculation if there are more than 2 data points: δ13C
    if count(t) > 2
        d13C_t = isotopes.d13C_org[t]
        # Error is generally 0.2‰. Create a vector the same length as d13C_t, fill with 0.2
        sigma_t = similar(d13C_t)
        for i in eachindex(sigma_t)
            sigma_t[i] = d13C_error
        end

        mindist = metropolis_min(nsteps, dist, d13C_t, sigma_t, burnin=5000)
        d13C_min[i] = nanmean(mindist)
        d13C_min_sigma[i] = nanstd(mindist)
    else
        d13C_min[i] = d13C_min_sigma[i] = NaN
    end

    # HC: only do metropolis calculation if there are more than 2 data points
    if count(s) > 2
        HC_s = kerogen.HC_Ratio[s]
        sigma_s = abs.(HC_s .* 0.01)  # Error at 1% of value

        mindist = metropolis_min(nsteps, dist, HC_s, sigma_s, burnin=5000)
        HC_min[i] = nanmean(mindist)
        HC_min_sigma[i] = nanstd(mindist)
    else
        HC_min[i] = HC_min_sigma[i] = NaN
    end

end

=#


















## old code old code old code!!
#=
# Global variables
global d13C_error = 0.2       # Strauss et al. (1992): uncertainty generally < 0.2‰


## -- Import data
isotopes_old = importdataset("data/Compiled_CarbonIsotopes_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)
kerogen_old = importdataset("data/Compiled_hcRatio_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)
schopf = importdataset("data/schopf_5-5_data.csv", ',', importas=:Tuple)


## -- Identify formations that have both H/C and δ13C_org data
# Set of formations that appear in both datasets
uniques = intersect(Set(isotopes_old.Std_Fm_Name), Set(kerogen_old.Std_Fm_Name))

# Create a new NamedTuple to store values
matches_old = (
    Std_Fm_Name = String[],
    Age = Float64[],
    d13C_org = Float64[],
    d13C_carb = Float64[],
    HC = Float64[],
    TOC = Float64[]
)

# For each element in match, get the average H/C ratio, δC13_org, and δC13_carb value
for i in uniques
    i_iso = findall(x -> x == i, isotopes_old.Std_Fm_Name)  # Returns an array of indices
    i_kgn = findall(x -> x == i, kerogen_old.Std_Fm_Name)

    # Average δC13_org, δC13_carb, and H/C values and add to NamedTuple
    push!(matches_old.Std_Fm_Name, i)
    push!(matches_old.Age, nanmean(isotopes_old.Age[i_iso]))
    push!(matches_old.d13C_org, nanmean(isotopes_old.d13C_org[i_iso]))
    push!(matches_old.d13C_carb, nanmean(isotopes_old.d13C_carb[i_iso]))
    push!(matches_old.HC, nanmean(kerogen_old.HC_Ratio[i_kgn]))
    push!(matches_old.TOC, nanmean(isotopes_old.TOC[i_iso]))
end

# Only values younger than 1.6Ga
t = (matches_old.Age .< 1600)


## -- Plot results
rel = plot(xlabel="H/C Ratio", ylabel="δC13", title="Known Relationship Between H/C ratio and δ13C",
    xlims=(0, 1.4), ylims=(-55, 15), framestyle=:box, size=(750,750),
    legend=:bottomright
)

# Show data from Schopf (1983) for comparison
s = (schopf.Age .< 1600)
plot!(rel, schopf.HC[.!s], schopf.d13C[.!s],  # Greater than 1.6 Ga
    color=:"#ffac1c", seriestype=:scatter, label="Schopf 1983 > 1.6 Ga",
    msc=:auto, shape=:utriangle
)

plot!(rel, schopf.HC[s], schopf.d13C[s],
    color=:"#ff5733", seriestype=:scatter, label="Schopf 1983 < 1.6 Ga",
    msc=:auto, shape=:utriangle
)

# Plot δ13C (org) data
plot!(rel, matches_old.HC[.!t], matches_old.d13C_org[.!t], # Greater than 1.6 Ga
    color=:"#038cfc", seriestype=:scatter, label="δ13C (org) > 1.6Ga", msc=:auto
)

plot!(rel, matches_old.HC[t], matches_old.d13C_org[t],     # Less than 1.6 Ga
    color=:"#003763", seriestype=:scatter, label="δ13C (org) < 1.6Ga", msc=:auto
)

# Plot δC13 (carb) data for comparison
plot!(rel, matches_old.HC, matches_old.d13C_carb,
    color=:"#3dbd46", seriestype=:scatter, label="δ13C (carb)", msc=:auto
)


## -- Fit Rayliegh-fractionation style equation to δ13C-vs-HC
using LsqFit

s = (schopf.Age .< 1600)

@. r₀(HC, p) = p[1]/(HC/p[2] + p[3])^(p[4]-1) + p[5]
p₀ = Float64[1, 2, 0, 1.3, -30]
lb = Float64[0, 0, 0, 1, -50]
ub = Float64[Inf, 10, 1, 10, 0]

fitray = curve_fit(r₀, schopf.HC[s], schopf.d13C[s], p₀, lower=lb, upper=ub)
display(fitray.param)

x = 0:0.01:1.5
plot!(rel, x, r₀(x, fitray.param), label="Rayleigh model")
display(rel)
#savefig("HC ratio and δ13C.pdf")


## -- Normalize data by subtracting out the average value for that 100Ma bin
# Resample d13C-org, HC values
org_rs = NamedTuple{(:c, :m, :e)}(bin_bsr_means(isotopes.Age, isotopes.d13C_org, 0, 3700, 37,
    x_sigma=isotopes.Uncertainty, y_sigma=isotopes.d13C_org * d13C_error, sem=:sigma
))
hc_rs = NamedTuple{(:c, :m, :e)}(bin_bsr_means(kerogen.Age, kerogen.HC_Ratio, 0, 3700, 37,
    x_sigma=kerogen.Uncertainty, y_sigma=kerogen.HC_Ratio*0.01, sem=:sigma
))

# Create a NamedTuple to store normalized data
# Each element in the tuple stores the value and the associated error
norm_match = (d13C_org = (val = [], e = []), HC = (val = [], e = []))

# Find index of appropriate age bin (should be same for isotopes and kerogen)
bi = findclosest(match.Age, org_rs.c)
bh = findclosest(match.Age, hc_rs.c)

# Subtract period d13C and H/C values, push value and error into the new tuple
append!(norm_match.d13C_org.val, (match.d13C_org .- org_rs.m[bi]))   # Values
append!(norm_match.HC.val, (match.HC .- hc_rs.m[bh]))

append!(norm_match.d13C_org.e, org_rs.e[bi])                             # Errors
append!(norm_match.HC.e, hc_rs.e[bh])


## -- Plot data
# Plot δ13C (org) data
normed = plot(xlabel="H/C Ratio", ylabel="δC13", title="normalized with resampled period average",
    framestyle=:box, size=(750, 750), legend=:bottomright
)

plot!(normed, norm_match.HC.val, norm_match.d13C_org.val,
    xerror=norm_match.HC.e, yerror=norm_match.d13C_org.e,
    color=:"#038cfc", label="", seriestype=:scatter, msc=:auto
)

# Add linear regression
# Find indices where both a and b are not NaN
t = all(!isnan, hcat(norm_match.HC.val, norm_match.d13C_org.val); dims=2) |>vec |> findall

# Initialize matrices a (x values) and b (y values)
a = [norm_match.HC.val[t] ones(length(norm_match.HC.val[t]))]
b = norm_match.d13C_org.val[t]

# Compute linear regression
aTa = *(transpose(a), a)
aTb = *(transpose(a), b)

# Convert data type
aTa = Array{Float64}(aTa)
aTb = Array{Float64}(aTb)

# Calculate slope, intercept and from those get vectors to plot
mb = aTa\aTb
lim = xlims(normed)
x = [lim[1], lim[2]]
y = mb[1] .*x .+ mb[2]

#Calculate R²
ssr_val = 0         # Initialize sum of residuals
sst_val = 0         # Initialize sum of squares
ybar = mean(norm_match.d13C_org.val[t])     # Mean of y values

# For each data point, ssr and sst
for i in t
    xi = norm_match.HC.val[i]
    yi = norm_match.d13C_org.val[i]

    # Find expected y value
    yhat = mb[1] * xi + mb[2]

    # Add to sums
    global ssr_val = ssr_val + (yi-yhat)^2
    global sst_val = sst_val + (yi-ybar)^2
end

r2_val = 1 - (ssr_val/sst_val)
r2_val = round(r2_val, digits=3)

# Add linear regression to plot
plot!(normed, x, y, line=:dash, label="linear regression \n R²=$r2_val")
display(normed)


## -- Find minimum using metropolis and subtract from each data point
# base metropolis code is from calculate_forg.jl, added HC calculation here
nsteps = 10000
dist = ones(10)

binedges = 0:100:3800
bincenters = cntr(binedges)
d13C_min = similar(bincenters)          #δ13C
d13C_min_sigma = similar(bincenters)
HC_min = similar(bincenters)            # H/C
HC_min_sigma = similar(bincenters)

# Iterate through each bin center
for i = 1:length(bincenters)
    # Get indexes for data that falls in the bin and is not NaN
    local t = @. (binedges[i] < isotopes.Age < binedges[i+1]) & !isnan(isotopes.d13C_org)
    local s = @. (binedges[i] < kerogen.Age < binedges[i+1]) & !isnan(kerogen.HC_Ratio)

    # δ13C: only do metropolis calculation if there are more than 2 data points: δ13C
    if count(t) > 2
        d13C_t = isotopes.d13C_org[t]
        # Error is generally 0.2‰. Create a vector the same length as d13C_t, fill with 0.2
        sigma_t = similar(d13C_t)
        for i in eachindex(sigma_t)
            sigma_t[i] = d13C_error
        end

        mindist = metropolis_min(nsteps, dist, d13C_t, sigma_t, burnin=5000)
        d13C_min[i] = nanmean(mindist)
        d13C_min_sigma[i] = nanstd(mindist)
    else
        d13C_min[i] = d13C_min_sigma[i] = NaN
    end

    # HC: only do metropolis calculation if there are more than 2 data points
    if count(s) > 2
        HC_s = kerogen.HC_Ratio[s]
        sigma_s = abs.(HC_s .* 0.01)  # Error at 1% of value

        mindist = metropolis_min(nsteps, dist, HC_s, sigma_s, burnin=5000)
        HC_min[i] = nanmean(mindist)
        HC_min_sigma[i] = nanstd(mindist)
    else
        HC_min[i] = HC_min_sigma[i] = NaN
    end

end

# Subtract metropolis calculated min from all values
# Create new tuple to hold normed values
norm_metro = (d13C = (val = [], err = []), HC = (val = [], err = []))

# Find index of closest bin
bi = findclosest(match.Age, bincenters)

# Subtract minimum δ13C and H/C values
append!(norm_metro.d13C.val, (match.d13C_org .- d13C_min[bi]))  # Values
append!(norm_metro.HC.val, (match.HC .- HC_min[bi]))

append!(norm_metro.d13C.err, d13C_min_sigma[bi])                # Error
append!(norm_metro.HC.err, HC_min_sigma[bi])


## -- Plot data
# Normalize both axes
plot_metro = plot(xlabel="H/C Ratio", ylabel="δC13", title="both axes normalized with metropolis_min",
    framestyle=:box, size=(750, 750), legend=:bottomright
)

plot!(plot_metro, norm_metro.HC.val, norm_metro.d13C.val,
    xerror=norm_metro.HC.err, yerror=norm_metro.d13C.err,
    color=:"#038cfc", label="", seriestype=:scatter, msc=:auto
)


# Normalize only δC13
metro_d13C = plot(xlabel="H/C Ratio", ylabel="δC13", title="δ13C normalized with metropolis_min",
framestyle=:box, size=(750, 750), legend=:bottomright
)

plot!(metro_d13C, match.HC, norm_metro.d13C.val,
    yerror=norm_metro.d13C.err,
    color=:"#038cfc", label="", seriestype=:scatter, msc=:auto
)


# Normalize only HC
metro_HC = plot(xlabel="H/C Ratio", ylabel="δC13", title="HC normalized with metropolis_min",
framestyle=:box, size=(750, 750), legend=:bottomright
)

plot!(metro_HC, norm_metro.HC.val, match.d13C_org,
    xerror=norm_metro.HC.err,
    color=:"#038cfc", label="", seriestype=:scatter, msc=:auto
)

# Display plots
display(plot_metro)
display(metro_d13C)
display(metro_HC)


## -- Plot all values (partially replicate plot from calculate_forg.jl)
# Plot all data (not resampled) and distinguish between data from different sources
plot_array = []

filters = ["All Data", "PMCID 1.1a", "Strauss 1992", "Krissansen-Totton 2015"]
# Add subplots
for i in filters
    # Filter duplicates for all plots, by external reference if applicable
    if i == "All Data"
        local t = @. (isotopes.Duplicate == "FALSE") & (isotopes.Circular_Citation == "FALSE")
    else
        local t = @. (isotopes.Duplicate == "FALSE") & (isotopes.Circular_Citation == "FALSE") & (isotopes.External_Ref == i)
    end

    active = plot(xlabel="Age (Ma)", ylabel="δC13", title="$i",
        ylims=(-60,25), xlims=(0, 3800), size=(600,1200), msc=:auto, framestyle=:box
    )
    plot!(active, isotopes.Age[t], isotopes.d13C_carb[t], label="δ13C carb", color=:red,
        msc=:auto, alpha=0.3, seriestype=:scatter, markersize=2
    )
    plot!(active, isotopes.Age[t], isotopes.d13C_org[t], label="δ13C carb", color=:blue,
        msc=:auto, alpha=0.3, seriestype=:scatter, markersize=2
    )

    push!(plot_array, active)
end
all_iso = plot(plot_array..., layout = (length(filters),1))
display(all_iso)
=#