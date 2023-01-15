## Matches H/C ratios and d13C-org values from the same formations
# Imports
using Plots
using StatGeochem
using Statistics
using Chron

# Global variables
global d13C_error = 0.2       # Strauss et al. (1992): uncertainty generally < 0.2‰


## -- Import data
isotopes = importdataset("data/Compiled_CarbonIsotopes_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)
kerogen = importdataset("data/Compiled_hcRatio_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)
schopf = importdataset("data/schopf_5-5_data.csv", ',', importas=:Tuple)


## -- Identify formations that have both H/C and δ13C_org data
# Set of formations that appear in both datasets
uniques = intersect(Set(isotopes.Std_Fm_Name), Set(kerogen.Std_Fm_Name))

# Create a new NamedTuple to store values
match = (
    Std_Fm_Name = String[],
    Age = Float64[],
    d13C_org = Float64[],
    d13C_carb = Float64[],
    HC = Float64[],
    TOC = Float64[]
)

# For each element in match, get the average H/C ratio, δC13_org, and δC13_carb value
for i in uniques
    i_iso = findall(x -> x == i, isotopes.Std_Fm_Name)  # Returns an array of indices
    i_kgn = findall(x -> x == i, kerogen.Std_Fm_Name)

    # Average δC13_org, δC13_carb, and H/C values and add to NamedTuple
    push!(match.Std_Fm_Name, i)
    push!(match.Age, nanmean(isotopes.Age[i_iso]))
    push!(match.d13C_org, nanmean(isotopes.d13C_org[i_iso]))
    push!(match.d13C_carb, nanmean(isotopes.d13C_carb[i_iso]))
    push!(match.HC, nanmean(kerogen.HC_Ratio[i_kgn]))
    push!(match.TOC, nanmean(isotopes.TOC[i_iso]))
end

# Only values younger than 1.6Ga
t = (match.Age .< 1600)


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
plot!(rel, match.HC[.!t], match.d13C_org[.!t], # Greater than 1.6 Ga
    color=:"#038cfc", seriestype=:scatter, label="δ13C (org) > 1.6Ga", msc=:auto
)

plot!(rel, match.HC[t], match.d13C_org[t],     # Less than 1.6 Ga
    color=:"#003763", seriestype=:scatter, label="δ13C (org) < 1.6Ga", msc=:auto
)

# Plot δC13 (carb) data for comparison
plot!(rel, match.HC, match.d13C_carb,
    color=:"#3dbd46", seriestype=:scatter, label="δ13C (carb)", msc=:auto
)


## -- Fit Rayliegh-fractionation style equation to δ13C-vs-HC
using LsqFit

s = (schopf.Age .< 1600)

@. r₀(HC, p) = p[1]/(HC/p[2] + p[3])^(p[4]-1) + p[5]
p₀ = Float64[1, 2, 0, 1.3, -30]
lb = Float64[0, 0, 0, 1, -50]
ub = Float64[Inf, 10, 1, 10, 0]

fit = curve_fit(r₀, schopf.HC[s], schopf.d13C[s], p₀, lower=lb, upper=ub)
display(fit.param)

x = 0:0.01:1.5
plot!(rel, x, r₀(x, fit.param), label="Rayleigh model")
display(rel)
savefig("HC ratio and δ13C.pdf")


## -- Repeat above using key system
# Compares individual samples rather than averages for formations


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