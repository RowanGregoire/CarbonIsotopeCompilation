## Get spatiotemporal metadata about the samples in isotopes.csv

#= 
   Print a general data report for unbinned data, as well as a report for each 100Ma bin

   For each report, print the number of organic and carbonate samples for each country / region, 
   and the percentage of total samples in bin that the number represents

   Example (all samples):

    country, total, % total, count organic, % total, count carbonate, % total
    Australia, 4507.0, 17.46, 2053.0, 22.41, 2454.0, 14.74
    Russia, 2103.0, 8.15, 272.0, 2.97, 1831.0, 11.0
    ...
    Pakistan, 30.0, 0.12, 30.0, 0.33, 0.0, 0.0

=#

using StatGeochem
using Test
isotopes = importdataset("isotopes.csv", ',', importas=:Tuple)

## Count by country
"""
    find_locats1(isotopes)

Count the number of samples in `isotopes` by country.

# Benchmark:
    BenchmarkTools.Trial: 541 samples with 1 evaluation.
    Range (min … max):  7.499 ms …  13.012 ms  ┊ GC (min … max): 0.00% … 0.00%
    Time  (median):     9.221 ms               ┊ GC (median):    0.00%
    Time  (mean ± σ):   9.243 ms ± 364.904 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

    Memory estimate: 6.58 KiB, allocs estimate: 16.

# Example:
    (country_list, count_org, count_carb) = find_locats1(isotopes)
"""
function find_locats1(isotopes)
    country = unique(isotopes.country)
    n = length(country)

    count_org = zeros(n)
    count_carb = zeros(n)

    # Loop through all countries
    for i in 1:lastindex(country)
        # For each sample, increase counter if it is not NaN and in the country
        for j in 1:lastindex(isotopes.key)
            if (isotopes.country[j] == country[i]) && !isnan(isotopes.d13c_org[j])
                count_org[i] += 1
            end

            if (isotopes.country[j] == country[i]) && !isnan(isotopes.d13c_carb[j])
                count_carb[i] += 1
            end
        end
    end

    return country, count_org, count_carb
end

## second version
"""
    find_locats2(isotopes)

Count the number of samples in `isotopes` by country.

# Benchmark:

    BenchmarkTools.Trial: 970 samples with 1 evaluation.
    Range (min … max):  4.042 ms … 38.214 ms  ┊ GC (min … max): 0.00% … 0.00%
    Time  (median):     4.775 ms              ┊ GC (median):    0.00%        
    Time  (mean ± σ):   5.140 ms ±  1.691 ms  ┊ GC (mean ± σ):  0.99% ± 4.99%

    Memory estimate: 827.81 KiB, allocs estimate: 378.

# Example:
    (country_list, count_org, count_carb) = find_locats2(isotopes)
"""
function find_locats2(isotopes)
    country = unique(isotopes.country)
    n = length(country)

    count_org = zeros(n)
    count_carb = zeros(n)

    for i in 1:lastindex(country)
        t = @. isotopes.country == country[i]
        count_org[i] = count(!isnan, isotopes.d13c_org[t])
        count_carb[i] = count(!isnan, isotopes.d13c_carb[t])
    end 
    
    return country, count_org, count_carb
end


## third version. very similar to version 1
"""
find_locats3(samples::Vector, places::Vector, orgs::Vector, carbs::Vector)

For each unique value in `places`, count the number of non `NaN` values in `orgs` and `carbs`. 

Also takes `samples`, which is intended as `isotopes.key`, but can be any vector that is the same 
length or longer as the largest number of non `NaN` values in `orgs` and `carbs`.

Some modifications have been made since running the first benchmark test (without `break`).

# Benchmark 
Without `break`:
    BenchmarkTools.Trial: 796 samples with 1 evaluation.
    Range (min … max):  5.619 ms …   9.812 ms  ┊ GC (min … max): 0.00% … 0.00%
    Time  (median):     6.276 ms               ┊ GC (median):    0.00%
    Time  (mean ± σ):   6.280 ms ± 322.188 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

    Memory estimate: 6.58 KiB, allocs estimate: 16.

With `break`
    BenchmarkTools.Trial: 2572 samples with 1 evaluation.
    Range (min … max):  1.639 ms …   5.072 ms  ┊ GC (min … max): 0.00% … 0.00%
    Time  (median):     1.949 ms               ┊ GC (median):    0.00%        
    Time  (mean ± σ):   1.941 ms ± 155.353 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

    Memory estimate: 6.61 KiB, allocs estimate: 17.

# Example:
    (country_list, count_org, count_carb) = find_locats3(isotopes.key, isotopes.country, isotopes.d13c_org, isotopes.d13c_carb)
"""
function find_locats3(samples::Vector, places::Vector, orgs::Vector, carbs::Vector)
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


## Make sure the functions do the same thing :)
(c1, co1, cc1) = find_locats1(isotopes)
(c2, co2, cc2) = find_locats2(isotopes)

(c3, co3, cc3) = find_locats3(isotopes.key, isotopes.country, isotopes.d13c_org, isotopes.d13c_carb)

@test c1 == c2 == c3       # Country
@test co1 == co2 == co3    # Organic
@test cc1 == cc2 == cc3    # Carbonate


## Print this data to the file
# I've changed how the functions return things so this won't actually work anymore

#=
io = open("metadata.txt", "w")
write(io, "country, total, % total, count organic, % total, count carbonate, % total\n")
all = locats.count_org .+ locats.count_carb
per_all = all ./ sum(all) .* 100
per_org = locats.count_org ./ sum(locats.count_org) .* 100
per_carb = locats.count_carb ./ sum(locats.count_carb) .* 100

for i in 1:lastindex(locats2.country)
    write(io, "$(locats2.country[i]), $(all[i]), $(round(per_all[i], digits = 2)), ")
    write(io, "$(locats2.count_org[i]), $(round(per_org[i], digits = 2)), ")
    write(io, "$(locats2.count_carb[i]), $(round(per_carb[i], digits = 2))\n")
end
close(io)
=#

## Write a function that calls find_locat but only passes a subsection of the compilation
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
    return find_locats3(key, country, orgs, carbs)
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

# I need to make a data structure that can hold this data and also the data from percentages
# Then I can iterate through the structure as I'm printing
# Try iterating and printing to the file at the same time...
periods = ["All Data", "Phanerozoic", "Neoproterozoic", "Mesoproterozoic", "Paleoproterozoic", "Archean"]
ageconstraints = ((0,4000), (0,541), (541,1000), (1000, 1600), (1600, 2500), (2500, 4000))

file = open("metadata.csv", "w")

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

