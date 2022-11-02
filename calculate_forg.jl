using StatGeochem, Plots
cd(@__DIR__)


isotopes = importdataset("data/Compiled_CarbonIsotopes_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)
t = @. (isotopes.Flagged != "X") & (isotopes.Duplicate == "FALSE")


kerogen = importdataset("data/Compiled_hcRatio_Data_UncertaintyUpdate.csv", ',', importas=:Tuple)

Δδ(r) = 4.05 - 3.05r + 0.785/r + 0.0165/r^2 - (8.79E-4)/r^3


## --- Resample and plot average H/C ratio


hk = plot(xlabel="Age [Ma]", ylabel="H/C ratio")
HC_ratio = NamedTuple{(:c,:m,:el,:eu)}(bin_bsr_means(kerogen.Age, kerogen.HC_Ratio, 0, 3800, 38, x_sigma=kerogen.Uncertainty, y_sigma=kerogen.HC_Ratio.*0.02))
plot!(hk, HC_ratio.c,HC_ratio.m, yerror=(HC_ratio.el,HC_ratio.eu), label="d13C carb", seriestype=:scatter, msc=:auto)


## -- Resample d13C_org
h = plot(xlabel="Age [Ma]", ylabel="d13 C", framestyle=:box)


nresamplings = 10000
d13C_org = NamedTuple{(:c,:m,:el,:eu)}(bin_bsr_means(isotopes.Age[t], isotopes.d13C_org[t], 0, 3800, 38, x_sigma=isotopes.Uncertainty[t], y_sigma=isotopes.d13C_org[t].*0.01))
plot!(h, d13C_org.c,d13C_org.m, yerror=(d13C_org.el,d13C_org.eu), label="d13C org", seriestype=:scatter, msc=:auto)


d13C_carb = NamedTuple{(:c,:m,:el,:eu)}(bin_bsr_means(isotopes.Age[t], isotopes.d13C_carb[t], 0, 3800, 38, x_sigma=isotopes.Uncertainty[t], y_sigma=isotopes.d13C_org[t].*0.01))
plot!(h, d13C_carb.c,d13C_carb.m, yerror=(d13C_carb.el,d13C_carb.eu), label="d13C carb", seriestype=:scatter, msc=:auto)


d13C_org_intial = d13C_org.m .- Δδ.(HC_ratio.m)
plot!(h, d13C_org.c, d13C_org_intial, label="Corrected d13c org")

plot!(isotopes.Age[t], isotopes.d13C_org[t], color=:blue, alpha=0.1, markersize=1, seriestype=:scatter, label="")
plot!(isotopes.Age[t], isotopes.d13C_carb[t], color=:red, alpha=0.1, markersize=1, seriestype=:scatter, label="")


## ---

d13C_mantle = -5.5

forg = @. (d13C_mantle - 0) / (d13C_org_intial - 0)

plot(d13C_carb.c, forg, xflip=true, xlabel="Age [Ma]", ylabel="f_org", label="", framestyle=:box, legend=:bottomleft)

forg = @. (d13C_mantle - d13C_carb.m) / (d13C_org_intial - d13C_carb.m)
plot!(d13C_carb.c, forg, xflip=true, xlabel="Age [Ma]", ylabel="f_org", label="", framestyle=:box, legend=:bottomleft)



des_marais = (;
    age=[2650.0, 2495.5516180173463, 2047.2862072181294, 1949.6316329385436, 1853.1940688240234, 1747.3141844633037, 1646.8618856663252, 1553.2220460691974, 1451.8744754266536, 1350.582859274457, 1251.9162470079896, 1051.5664666811736, 957.5075512657118, 850.7471608277119, 756.0325287284863, 656.6550224336061],
    forg=[0.08958593677142582, 0.10020889676396522, 0.1494628368926606, 0.19049706238925668, 0.17004010071808257, 0.11977338431409115, 0.12975286766763028, 0.14035064813951315, 0.14039261400727404, 0.14105567471789607, 0.16009238967121547, 0.1497514947102282, 0.19076697804304346, 0.2210044219590288, 0.21133867232428727, 0.14991501911778413]
)

plot!(des_marais.age, des_marais.forg, label="des marais", color=:black)



## ----

strauss_kerogen = importdataset("data/SourceFiles/strauss_1992_17.10_orgHC_Summary.csv", ',', importas=:Tuple)
plot(strauss_kerogen.H_C, strauss_kerogen.delC13_PDB, seriestype=:scatter, xlabel="H/C", ylabel="d13C")


# ## --- Resample age and organic carbon
#
# t = @. !(isnan(isotopes.Age) | isnan(isotopes.d13C_org))
# data = unelementify(isotopes, (:Age, :d13C_org), floatout=true, rows=t)
# sigma = abs.(data.*0.01)
# sigma[:,1] .= isotopes.Uncertainty[t]
#
# resampled = elementify(bsresample(data, sigma, 10^6), ("Age", "d13C_org"))

## --- Calculate minimum real carbon isotope value for age bins

using Chron
nsteps = 10000
dist = ones(10)

binedges = 0:100:3800
bincenters = cntr(binedges)
d13C_min = similar(bincenters)
d13C_min_sigma = similar(bincenters)

for i = 1:length(bincenters)
    t = (binedges[i] .< isotopes.Age .< binedges[i+1]) .& !isnan.(isotopes.d13C_org)

    if count(t) > 2
        d13Cₜ = isotopes.d13C_org[t]
        sigmaₜ = abs.(d13Cₜ .* 0.01)
        # sigma = fill(2., size(d13C)) # ten per mil uncertainty?

        mindist = metropolis_min(nsteps, dist, d13Cₜ, sigmaₜ, burnin=5000)
        d13C_min[i] = nanmean(mindist)
        d13C_min_sigma[i] = nanstd(mindist)
    else
        d13C_min[i] = d13C_min_sigma[i] = NaN
    end
end

plot!(bincenters, d13C_min, yerror=d13C_min_sigma)

## -- Calculate r
