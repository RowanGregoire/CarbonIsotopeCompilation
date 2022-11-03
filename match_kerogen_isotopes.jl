## Matches H/C ratios and d13C-org values from the same formations
# Imports 
using Plots
using StatGeochem
using Statistics


## -- Import data
isotopes = importdataset("data/Compiled_CarbonIsotopes_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)
kerogen = importdataset("data/Compiled_hcRatio_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)
schopf = importdataset("data/schopf_5-5_data.csv", ',', importas=:Tuple)


## -- Identify formations that have both H/C and δ13C_org data
# Set of formations that appear in both datasets
uniques = intersect(Set(isotopes.Std_Fm_Name), Set(kerogen.Std_Fm_Name))

# Create a new NamedTuple to store values
match = (Std_Fm_Name = [], Age = [], d13C_org = [], d13C_carb = [], HC = [])

# For each element in match, get the average H/C ratio, δC13_org, and δC13_carb value
for i in uniques
    i_iso = findall(x -> x == i, isotopes.Std_Fm_Name)  # Returns an array of indices
    i_kgn = findall(x -> x == i, kerogen.Std_Fm_Name)

    # Average δC13_org, δC13_carb, and H/C values and add to NamedTuple
    push!(match.Std_Fm_Name, i)
    push!(match.Age, mean(filter(!isnan, isotopes.Age[i_iso])))
    push!(match.d13C_org, mean(filter(!isnan, isotopes.d13C_org[i_iso])))
    push!(match.d13C_carb, mean(filter(!isnan, isotopes.d13C_carb[i_iso])))
    push!(match.HC, mean(filter(!isnan, kerogen.HC_Ratio[i_kgn])))
    
end

# Only values younger than 1.6Ga
t = (match.Age .< 1600)


## -- Plot results
rel = plot(xlabel="H/C Ratio", ylabel="δC13",
    xlims=(0, 1.4), ylims=(-55, 15), framestyle=:box, size=(500,500), 
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


## -- Normalize data
# Resample d13C-org, HC values
org_rs = NamedTuple{(:c, :m, :e)}(bin_bsr_means(isotopes.Age, isotopes.d13C_org, 0, 3700, 37,
    x_sigma=isotopes.Uncertainty, y_sigma=isotopes.d13C_org*0.01, sem=:sigma
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

# bi -> for each element in match, what is the index of the nearest bin

# Subtract period d13C and H/C values, push value and error into the new tuple
append!(norm_match.d13C_org.val, (match.d13C_org[bi] .- org_rs.m[bi]))   # Values
append!(norm_match.HC.val, (match.HC[bh] .- hc_rs.m[bh]))
                             
append!(norm_match.d13C_org.e, org_rs.e[bi])                             # Errors
append!(norm_match.HC.e, hc_rs.e[bh])    


## -- Plot data
# Plot δ13C (org) data
normed = plot(xlabel="H/C Ratio", ylabel="δC13",
    framestyle=:box, size=(500,500), legend=:bottomright
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

# Calculate slope, intercept and get vectors to plot
mb = aTa\aTb
lim = xlims(normed)
x = [lim[1], lim[2]]
y = mb[1] .*x .+ mb[2]

#Calculate R²
ssr = 0         # Initialize sum of residuals
sst = 0         # Initialize sum of squares
ybar = mean(norm_match.d13C_org.val[t])     # Mean of y values

# For each data point, ssr and sst
for i in t
    xi = norm_match.HC.val[i]
    yi = norm_match.d13C_org.val[i]
    
    # Find expected y value
    yhat = mb[1] * xi + mb[2]

    # Add to sums
    ssr = ssr + (yi-yhat)^2
    sst = sst + (yi-ybar)^2
end

r2 = 1 - (ssr/sst)
r2 = trunc(r2, digits=3)

# Add linear regression to plot
plot!(normed, x, y, line=:dash, label="linear regression \n R²=$r2")
