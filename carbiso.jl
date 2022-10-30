#####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     
# Takes code from BootstrapResamplingDemo.ipynb to resample carbon isotope and H/C data

# Load packages
try
  using StatGeochem
catch
  using Pkg
  Pkg.add("StatGeochem")
  using StatGeochem
end

using Statistics, StatsBase, DelimitedFiles
using Plots; gr();


#####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     
# Import data and prepare for resampling

# Import data as dictionary
println("importing data")

d13C = importdataset(raw"C:\\Users\\katie\\OneDrive\\Desktop\\Rocks\\Data\\Compiled_CarbonIsotopes_Data_UncertaintyUpdate.csv", ',', importas=:Dict);
HC = importdataset(raw"C:\\Users\\katie\\OneDrive\\Desktop\\Rocks\\Data\\Compiled_hcRatio_Data_UncertaintyUpdate.csv", ',', importas=:Dict);

# Create a list of all elements (must match column names in .csv file) we want to resample
# Elements that already have uncertanties must be listed first
d13Celements = ["Age", "Lat", "Long", "d13C-org", "d13C-carb"];
HCelements = ["Age", "Lat", "Long", "HC_Ratio"];

println("adding uncertainties")
# Add age uncertainties
d13C["Age_sigma"] = d13C["Uncertainty"];
HC["Age_sigma"] = HC["Uncertainty"];

# Add uncertainties - 1% for all elements with unspecified uncertainty (i.e. all except "Age")
for i=2:length(d13Celements)  # d13C data
  d13C[d13Celements[i]*"_sigma"] = d13C[d13Celements[i]] * 0.01
end

for i=2:length(HCelements)    # HC data
  HC[HCelements[i]*"_sigma"] = HC[HCelements[i]] * 0.01
end


#####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     
# Resample data

# Compute proximity coefficients (inverse weights)
println("resampling data")
k = invweight(d13C["Lat"], d13C["Long"], d13C["Age"])
k = invweight(d13C["Lat"], d13C["Long"], d13C["Age"])

# Probability of keeping a given data point when sampling
pd13C = 1.0./((k.*median(5.0./k)) .+ 1.0) # Keep rougly one-fith of the data in each resampling
pHC = 1.0./((k.*median(5.0./k)) .+ 1.0)

# Resample a few hundred times (all elements!)
nresamplings = 1000
mcd13C = bsresample(d13C, nresamplings*length(d13C["Lat"]), d13Celements, pd13C)
mcHC = bsresample(HC, nresamplings*length(HC["Lat"]), HCelements, pHC)

(c,m,el,eu) = bin_bsr_means(d13C["Age"],d13C["d13C-org"],0,3700,37, p=pd13C, x_sigma=d13C["Age_sigma"], nresamplings=10000)
display(plot(c,m,yerror=(el,eu),label="",xlabel="Age", ylabel="new resample d13C-org", framestyle=:box))


#####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####   
## --- This is a code cell  

println("this is a code cell!")
c = 1+1
println("Here is the results of your calculation: $c")

## ---

# Adjust d13C-org values down based on H/C ratios
println("Adjusting d13C-org values")

# Calculate mean d13C-org and H/C for 0-3700Ma, 37 bins, from resampled dataset
# (c = bin centers, m = means, e = 1-sigma S.E.M)
(c_d13C,m_d13C,e_d13C) = binmeans(mcd13C["Age"],mcd13C["d13C-org"],0,3700,37; resamplingratio=nresamplings)
(c_HC,m_HC,e_HC) = binmeans(mcHC["Age"],mcHC["HC_Ratio"],0,3700,37; resamplingratio=nresamplings)

# For each binned mean in d13C, adjust for the corresponding H/C ratio (from Des Marais 1992)
# How do I make sure that the ratios match? i.e. that bin 34 for d13C is the same time period as bin 34 for HC
adj_org = []  # Create a new array for the adjusted means

for i=1:length(m_d13C)
  r = m_HC[i];
  delta_org = 4.05 - 3.05*r + (0.785/r) + (0.0165/(r^2)) - (8.79e-4/(r^3));

  # Adjust
  push!(adj_org, m_d13C[i] - delta_org)
end

#=
I want some method to adjust the values down where I get a lot of data points so I can get standard deviation etc...
Also this does not generate a plot that looks anything at all similar to Des Marais so maybe think about that
=#


#####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     #####     
# Plot results using the resampled dataset
println("plotting data")

# d13C-org
display(plot(c_d13C,m_d13C,yerror=2*e_d13C,label="",xlabel="Age (Ma)", ylabel="d13C_org",framestyle=:box))
println("here - d13C")

# H/C ratio
display(plot(c_HC,m_HC,yerror=2*e_HC,label="",xlabel="Age (Ma)", ylabel="H/C ratio",framestyle=:box))
println("here - H/C")

# Adjusted d13C-org
# Does not at all match the Des Marais paper lmao
# Error bars? Probably from using a different method
display(plot(c_d13C, adj_org,label="",xlabel="Age (Ma)", ylabel="Adjusted d13C-org ratio", framestyle=:box))
println("here - adjusted d13C")
