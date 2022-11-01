## Matches H/C ratios and d13C-org values from the same formations
# Imports 
using Plots
using StatGeochem
using Statistics

## -- Import data
isotopes = importdataset("data/Compiled_CarbonIsotopes_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)
kerogen = importdataset("data/Compiled_hcRatio_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)

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
    xlims=(0, 1.4), ylims=(-55, 15), framestyle=:box, legend=:bottomright
)

plot!(rel, match.HC[.!t], match.d13C_org[.!t], # Greater than 1.6 Ga
    color=:"#038cfc", seriestype=:scatter, label="δ13C (org) older than 1.6Ga", msc=:auto
)

plot!(rel, match.HC[t], match.d13C_org[t],     # Less than 1.6 Ga
    color=:"#003763", seriestype=:scatter, label="δ13C (org) younger than 1.6Ga", msc=:auto
)

# Plot δC13 (carb) data for comparison
plot!(rel, match.HC, match.d13C_carb,
    alpha=:0.5, color=:green,
    seriestype=:scatter, label="δ13C (carb)", msc=:auto
)
