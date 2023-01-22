#=
    Functions used for carbon isotope compilation intake and initial analysis!
=#

## -- Load the packages we'll use
using StatGeochem
using Statistics
using Chron
using LsqFit

## -- Fuctions!
"""
    convert_vpdb!(d13c_vpdb)

Converts δ13C measurement `d13c_vpdb` relative to the VPDB standard with the measurement 
relative to the PDB standard

# Example

intake one value that's relative to vpdb and return relative to pdb. can do this in the repl
"""
function convert_vpdb!()
    #=
        Do something
    =#
end


"""
    build_match(kerogen, isotopes)

Find common samples between data compilations `kerogen` and `isotopes`. Return a `Tuple` of 
samples with kerogen and δ13C organic values where both kerogen and organic carbon data are 
not `NaN`

This algorithm assumes that there are no duplicate keys. This is not true; distinct analyses
done on the same sample (e.g. isotopic analysis of bulk TOC vs. kerogen vs. aromatic bitumen) 
have the same key to indicate that they are from the same sample. In theory, one could write
an algorithm to ID duplicate keys and prefer one isotope source over another. Currently, this one 
finds whichever analysis is listed first in the dataset.

# Example

This function requires `kerogen` and `isotopes` to have elements `key`, `hc` and `d13c_org`.

I could probably change the arguments of the function to make it more general...

    julia> kerogen
    NamedTuple with 2 elements:
      key                     = Vector{Float64}(455,)       [3321.0 ... 16593.0]
      hc                      = Vector{Float64}(455,)       [1.3 ... 0.04]

    julia> isotopes
    NamedTuple with 29 elements:
      key                         = Vector{Float64}(24239,) [1.0 ... 16593.0]
      d13c_org                    = Vector{Float64}(24239,) [-32.1 ... -18.2]

    julia> build_match(kerogen, isotopes)
    NamedTuple with 3 elements:
      key       = Vector{Int64}(287,)       [3321 ... 16593]
      d13c_org  = Vector{Float64}(287,)     [-36.2 ... -18.2]
      hc        = Vector{Float64}(287,)     [1.3 ... 0.04]
"""
function build_match(kerogen, isotopes)
    @assert typeof(isotopes) <: NamedTuple
    @assert typeof(kerogen) <: NamedTuple

    # Find common sample keys between the kerogen and isotope datasets
    keymatch = intersect(kerogen.key, isotopes.key)
    matches = (key = Int64[], d13c_org = Float64[], hc = Float64[])

    # Excluding NaNs, find the H/C and d13Corg values for each key
    for i in 1:lastindex(keymatch)
        keyi = keymatch[i]
        orgi = isotopes.d13c_org[findmatches(keyi, isotopes.key)]
        hci = kerogen.hc[findmatches(keyi, kerogen.key)]

        if !isnan(orgi) && !isnan(hci)
            push!(matches.key, keyi)
            push!(matches.d13c_org, orgi)
            push!(matches.hc, hci)
        else
            continue
        end
        
    end

    return matches
end


"""
    fit_rayleigh(d13c, hc)

Fit a Rayleigh-style fractionation curve to δ13C org `d13c` and H/C ratio `hc` data.

Currently returns a weird domain error that I can't track down :(

# Example

    (x, y) = fit_rayleigh(d13c, hc)
"""
function fit_rayleigh(d13c, hc)
    @. r₀(HC, p) = p[1]/(HC/p[2] + p[3])^(p[4]-1) + p[5]
    p₀ = Float64[1, 2, 0, 1.3, -30]
    lb = Float64[0, 0, 0, 1, -50]
    ub = Float64[Inf, 10, 1, 10, 0]

    fitted = curve_fit(r₀, hc, d13c, p₀, lower=lb, upper=ub)

    return 0:0.01:maximum(hc), r₀(x, fitted.param)
end


"""
    sync_vectors(data1, data2)

Replace values in `data1` that do not have a corresponding value in `data2` with `NaN`. 
Returns the new versions of `data1` and `data2`.

# Example

Removing data values that do not have an uncertainty:

    (new_data, new_uncertainty) = sync_vectors(data, uncertainty)    
"""
function sync_vectors(data1::Vector{Float64}, data2)
    @assert length(data1) == length(data2)

    # We want to make a copy of the data, not modify the original data
    localdata1 = deepcopy(data1)
    localdata2 = deepcopy(data2)

    # Remove values in data1 where data2 is NaN
    for i in 1:lastindex(localdata2)
        if isnan(localdata1[i])
            continue
        elseif isnan(localdata2[i])
            localdata1[i] = NaN
        end
    end

    return localdata1, localdata2
end

function sync_vectors(data1::Vector, data2)
    @assert length(data1) == length(data2)

    # We want to make a copy of the data, not modify the original data
    localdata1 = deepcopy(data1)
    localdata2 = deepcopy(data2)

    # Must be of type Float64 to add NaNs
    localdata1 = convert(Vector{Float64}, localdata1)

    # Remove values in data1 where data2 is NaN
    for i in 1:lastindex(localdata2)
        if isnan(localdata1[i])
            continue
        elseif isnan(localdata2[i])
            localdata1[i] = NaN
        end
    end

    return localdata1, localdata2
end


"""
    count_data(age::Vector, data::Vector, edges::AbstractRange)

Count the number of datapoints in `data` between each bin as defined by `edges`. Return
a vector of the counts and bin centers.
"""
function count_data(age::Vector, data::Vector, edges::AbstractRange)
    nbins = length(edges) - 1
    bin_width = step(edges)
    counts = zeros(nbins)

    # Count data in each bin
    for i in 1:lastindex(data)
        if !isnan(data[i]) && !isnan(age[i])
            bin_index = convert(Int, age[i] ÷ bin_width) + 1
            counts[bin_index] += 1
        end
    end

    return counts, cntr(edges)
end


"""
    knit_vectors(source)


    knit_vectors(source)


    fmt_data(source)

Given a vector of vectors `source`, return a matrix where each column is one of the vectors
from `source`

# Example

    julia> data1 = [1,2,3,4,5];

    julia> data2 = [6,7,8,9,10];

    julia> data3 = [data1, data2]
    2-element Vector{Vector{Int64}}:
    [1, 2, 3, 4, 5]
    [6, 7, 8, 9, 10]

    julia> knit_vectors(data3)
    5x2 Matrix{Float64}:
    1.0   6.0
    2.0   7.0
    3.0   8.0
    4.0   9.0
    5.0  10.0
"""
function knit_vectors(source)
    # Make sure all vectors in source are the same length
    #=
    @assert ... ?
    =#

    # Initialize target array
    ncols = length(source)
    nrows = length(source[1])
    target = zeros(nrows, ncols)

    # For each vector in source, add it to a column in target
    for j in 1:ncols
        for i in 1:nrows
            target[i, j] = source[j][i]
        end
    end

    return target
end


"""
    Δδ(r, params)

Use the parameters `param` from the Rayleigh-style fractionation curve to determine the
change in δ13C org values based on the H/C ratios `r`

Currently uses the Δδ equation defined by Des Marais et al. 1992
"""
function Δδ(r, param)
    #=
        Do something with the parameters from the Rayleigh fractionation...

        At least try to figure this one out given the parameters of the one that works
        Should I do this with the resampled H/C dataset? If we're using the non-resampled
        dataset to create the curve, it seems like we shouldn't need to resample the H/C 
        data... or do we need it for something else?
    =#

    # Equation from Des Marais 1992
    return 4.05 - 3.05r + 0.785/r + 0.0165/r^2 - (8.79E-4)/r^3
end