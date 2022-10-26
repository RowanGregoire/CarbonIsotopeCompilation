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
plot(d13C_carb.c, forg, xflip=true, xlabel="Age [Ma]", ylabel="f_org")
