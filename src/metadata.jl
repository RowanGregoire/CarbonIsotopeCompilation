#= 
Metadata - for all data, as well for some as arbitrary bins, print the number of organic and carbonate
isotopic data by country or region. Calculate the percentage of the total that the number of samples from
that region represents. Print the output to a .csv

Location uncertainties - given a list of latitude and longitudes which represent the extreme points of a country, calculate
the center of that country and the arc-degree uncertainty for that center

=#

using StatGeochem
using Test
isotopes = importdataset("data/compilation.csv", ',', importas=:Tuple)


"""
find_locats(samples::Vector, places::Vector, orgs::Vector, carbs::Vector)

For each unique value in `places`, count the number of non `NaN` values in `orgs` and `carbs`. 

Also takes `samples`, which is intended as `isotopes.key`, but can be any vector that is the same 
length or longer as the largest number of non `NaN` values in `orgs` and `carbs`.


# Example:
    (country_list, count_org, count_carb) = find_locats(isotopes.key, isotopes.country, isotopes.d13c_org, isotopes.d13c_carb)
"""
function find_locats(samples::Vector, places::Vector, orgs::Vector, carbs::Vector)
    placeunique = unique(places)
    n = length(placeunique)

    count_org = zeros(n)
    count_carb = zeros(n)

    # Iterate through every sample
    for i in 1:lastindex(samples)
        # For each sample, loop through the list of countries 
        for j in 1:lastindex(placeunique)
            # Stop current iteration and move to next if the country does not match
            if places[i] !== placeunique[j]
                continue

            # If the country matches, check if the sample is not NaN
            else
                if !isnan(orgs[i])
                    count_org[j] += 1
                end

                if !isnan(carbs[i])
                    count_carb[j] += 1
                end

                # Move to the next iteration - we know there will be no more matches
                # Fun fact! Adding this makes the function ~6x faster
                break
            end
        end
    end

    return placeunique, count_org, count_carb
end


"""
    getmetadata(compilation, maxage, minage)

Get sample location metadata from `compilation` between an inclusive lower bound `maxage` and exclusive upper 
bound `maxage`.

"""
function getmetadata(compilation, minage, maxage)
    # Get a boolean array to index into compilation at the right ages
    t = @. minage <= compilation.age < maxage

    # Get the right vectors at the right ranges and index into the compilation
    key = compilation.key[t]
    country = compilation.country[t]
    orgs = compilation.d13c_org[t]
    carbs = compilation.d13c_carb[t]

    # Find and return location metadata
    return find_locats(key, country, orgs, carbs)
end


## Get metadata for the Phanerozoic, Proterozoic, and Archean and print to a .csv
#=
Age divisions and bounds
    Phanerozoic         (0, 541]
    Neoproterozoic      (541, 1000]
    Mesoproterozoic     (1000, 1600]
    Paleoproterozoic    (1600, 2500]
    Archean             (2500, 4000]
=#

periods = ["All Data", "Phanerozoic", "Neoproterozoic", "Mesoproterozoic", "Paleoproterozoic", "Archean"]
ageconstraints = ((0,4000), (0,541), (541,1000), (1000, 1600), (1600, 2500), (2500, 4000))

file = open("output/metadata.csv", "w")

for i in 1:lastindex(ageconstraints)
    (country, org, carb) = getmetadata(isotopes, ageconstraints[i][1], ageconstraints[i][2])

    # Percentage of total data for organic and carbonate sample counts
    percentorg = round.((org ./ sum(org) * 100), digits=2)
    percentcarb = round.((carb ./ sum(carb) * 100), digits=2)

    # Total samples for each country
    total = org .+ carb
    percenttotal = round.((total ./ sum(total) * 100), digits=2)

    # Total samples in each time periods
    totalperiod = sum.([org, carb])

    # Headings
    write(file, periods[i] * "\n")
    write(file, "total samples in time period [organic, carbonate] $totalperiod \n")
    write(file, "country, count organic, % total org, count carbonate, % total carb, all samples, % total\n")

    # Data
    for j in 1:lastindex(country)
        countryj = country[j]
        orgj = string(org[j])
        percentorgj = string(percentorg[j])
        carbj = string(carb[j])
        percentcarbj = string(percentcarb[j])
        totalj = string(total[j])
        percenttotalj = string(percenttotal[j])

        write(file, "$countryj, $orgj, $percentorgj, $carbj, $percentcarbj, $totalj, $percenttotalj")
        write(file, "\n")
    end
    
    write(file, "\n")
end

close(file)


## -- Calculate location centers and uncertainties
latlons = importdataset("data/locations.csv", ',', importas=:Tuple)
file = open("output/updated_locations.csv", "w")
write(file, "Country, Latitude (deg), Longitude (deg), Uncertainty (arc-deg)\n")

# Calculate center and uncertainty and put into a .csv
for i in eachindex(latlons.Location)
    lats = [latlons.lat1[i], latlons.lat2[i], latlons.lat3[i], latlons.lat4[i]]
    lons = [latlons.lon1[i], latlons.lon2[i], latlons.lon4[i], latlons.lon4[i]]
    (latctr, lonctr, uncert) = dist_uncert(lats, lons)

    write(file, "$(latlons.Location[i]), $latctr, $lonctr, $uncert\n")
end

close(file)