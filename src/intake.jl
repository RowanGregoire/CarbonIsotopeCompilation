## --- Functions used for carbon isotope compilation intake and initial analysis!


## -- Load the packages we'll use
using StatGeochem
using Statistics
using Chron
using LsqFit
using ProgressMeter

## -- Fuctions!

@. r₀(HC, p) = p[1]/(HC/p[2] + p[3])^(p[4]-1) + p[5]

"""
    fit_rayleigh(d13c, hc)

Fit a Rayleigh-style fractionation curve to δ13C org `d13c` and H/C ratio `hc` data.

# Example

    (x, y) = fit_rayleigh(d13c, hc)
"""
function fit_rayleigh(d13c, hc)
  p₀ = Float64[1, 2, 1e-3, 1.3, -30]
  lb = Float64[0, 1e-3, 1e-3, 1, -50]
  ub = Float64[Inf, 10, 1, 10, 0]

  fitted = curve_fit(r₀, hc, d13c, p₀, lower=lb, upper=ub)
  return fitted.param
end

function rayleigh_curve(p, hc)
  x = 0:0.01:maximum(hc)
  y = r₀(x, p)
  return x, y
end


"""
    sync_vectors(data1, data2)

Replace values in `data1` that do not have a corresponding value in `data2` with `NaN`. 
Returns the new versions of `data1` and `data2`.

kgn_age = (age = copy(kerogen.age), uncert = copy(kerogen.age_uncertainty))

# Example

Removing data values that do not have an uncertainty:

    (new_data, new_uncertainty) = sync_vectors(data, uncertainty)    
"""
function sync_vectors(data1::Vector{Float64}, data2)
    @assert length(data1) == length(data2)

    # We want to make a copy of the data, not modify the original data
    localdata1 = copy(data1)
    localdata2 = copy(data2)

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
    localdata1 = copy(data1)
    localdata2 = copy(data2)

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
```julia
add_uncert!(data, data_u;
    err = 0.05,
    method = :percent,
    addnote = false, 
    flag::Vector{String} = fill("", length(data)), 
    comment::Vector{String} = fill("", length(data)),
    note::String = "added est. uncertainty"
  )
```

Add estimated uncertainty to samples in `data` where the corresponding uncertainty in
`data_u` is `NaN`.

Optional keyword arguments and defaults:

  err = 0.05

The uncertainty to add or multiply to the data.

  method = :percent

How the uncertainty is calculated. If `:percent` is chosen, error will the value in `data`
multiplied by `err`. If `:plusminus` is chosen, error will be `err` with no other operations.

  addnote = false

Option to add a flag and comment to samples with estimated uncertainty.

  flag::Vector{String} = fill("", length(data))

If `addnote` is `true`, add "X" to the `flag` vector. If no `flag` vector specified, create a blank one.

  comment::Vector{String} = fill("", length(data))

If `addnote` is `true`, add `note` to the comment vector. If no `comment` vector specified, create a blank one.

  note::String = "added est. uncertainty"

If `addnote` is `true`, specify a comment to add to the `comment` vector.

### Example

add_uncert!(x, x_uncertainty, err=0.05, method=:percent, addnote=true, 
    flag=newdata.flag, comment=newdata.comments, note="est. 5% age uncertainty"
)
"""
function add_uncert!(data, data_u;
    err = 0.05,
    method = :percent,
    addnote = false, 
    flag::Vector{String} = fill("", length(data)), 
    comment::Vector{String} = fill("", length(data)),
    note::String = "added est. uncertainty"
  )
  for i in 1:lastindex(data_u)
      if isnan(data_u[i])
        if method == :percent
          # Percentage of the data value
          data_u[i] = data[i] * err
        elseif method == :plusminus
          # ± set number
          data_u[i] = err
        else
          error("method must be :percent or :plusminus")
        end

        # Add notes to samples
        if addnote
          # Add flag, if necessary
          if flag[i] != "X"
            flag[i] = "X"
          end

          # Add or append to existing comment
          if comment[i] == ""
            comment[i] = comment[i] * note
          else
            comment[i] = comment[i] * "; " * note
          end
        end
      end
  end
end


"""
    count_data(age::Vector, data::Vector, edges::AbstractRange)

Count the number of datapoints in `data` between each bin as defined by `edges`. Return
a vector of the counts and bin centers.
"""
function count_data(age::Vector, data::Vector, edges::AbstractRange)
    nbins = length(edges) - 1
    bin_width = step(edges)
    agemin = first(edges)
    counts = zeros(nbins)

    # Count data in each bin
    for i in 1:lastindex(data)
        if !isnan(data[i]) && !isnan(age[i])
            bin_index = convert(Int, (age[i] - agemin) ÷ bin_width) + 1
            counts[bin_index] += 1
        end
    end

    return counts, cntr(edges)
end


"""
    Δδ(r, params)

Use the parameters `param` from the Rayleigh-style fractionation curve to determine the
change in δ13C org values based on the H/C ratios `r`

Currently uses the Δδ equation defined by Des Marais et al. 1992
"""
function Δδ(hc, param::Symbol)
    #=
       Des Marais 1992 correction based on H/C ratio
    =#
    @assert param===:DesMarais

    # Equation from Des Marais 1992
    return @. 4.05 - 3.05hc + 0.785/hc + 0.0165/hc^2 - (8.79E-4)/hc^3
end


function Δδ(hc, param::Vector, hc₀::Number=1.5)
  #=
      Given an H/C ratio, what is the correction that we should make?
      We know from the Rayleigh fractionation curve we built earlier that lower H/C ratios
      mean more correction      
  =#

  @assert length(param) == 5
  return r₀(hc, param) .-  r₀(hc₀, param)
end


"""
```julia
(x::Vector, y::Vector, edges::LinRange)
```

Returns the means `m` and 1-σ uncertainties `e` for a variable `y` binned by independent variable `x`
into bins defined by `edges`.

### Example
```
(c, m, e) = bin_means(x, y, edges)
```
"""
function bin_means(x::Vector, y::Vector, edges)
    nbins = length(edges) - 1
    bin_width = step(edges)
    xmin = first(edges)

    # Preallocate
    sums = zeros(Float64, nbins)
    counts = zeros(Int, nbins)
    ordered = Array{Float64}(undef, length(y), nbins)
    sigma = zeros(Float64, nbins)
    
    # Fill sums and counts vectors to calculate means. Add data to stdevs array
    # Good canidate for @inbounds but only AFTER it's been bug tested so it doesn't segfault?
    for i in eachindex(x)
        if !isnan(x[i]) && !isnan(y[i])
            bin_index = Int((x[i] - xmin) ÷ bin_width) + 1

            # Make sure bin_index is valid
            if 1 <= bin_index <= nbins
                sums[bin_index] += y[i]
                counts[bin_index] += 1


                # Keeping the data points in this array means we can index into a specific
                # bin later to compute the standard deviation
                ordered[counts[bin_index], bin_index] = y[i]
            end
        end
    end

    # Calculate standard deviation for each bin. If there is no data, set error to NaN
    for i in 1:nbins
        sigma[i] = ifelse(counts[i] == 0, NaN, std(ordered[1:counts[i], i]))   # I don't think this is grabbing the things I want it to
    end

    return (edges[1:end-1] + edges[2:end])/2, sums ./ counts, sigma
end


"""
```julia
estimate_hc(age, hc, bound::Number=10)
```

Estimates H/C ratios for elements in `hc` that are `NaN`. For each `NaN` element
in `hc`, randomly select a non-`NaN` H/C value from the points with an age weighting
`bound` from the age of the sample.
"""
function estimate_hc!(age, hc, bound::Number=10)
  p = Progress(length(age) ÷ 10 , desc = "Finding H/C Ratios: ")

  @inbounds for i in eachindex(hc)
      if isnan(hc[i])
          # Get all samples in the age range
          target_age = age[i]
          t = @. (target_age - bound) < age < (target_age + bound)

          # Make sure the bound has data
          s = @. !isnan(hc[t])
          newbound = bound
          while isempty(hc[t][s])
              # L + ratio + no data in your subset
              # Is there a way to make this recursive? seems to throw me into a infinite loop
              # Is that even useful?
              newbound += 5
              t = @. (target_age - newbound) < age < (target_age + newbound)
              s = @. !isnan(hc[t])
          end

          # Randomly pick one and assign it to be the H/C value for that sample
          hc[i] = rand(hc[t][s])
      end
      (i % 10 == 0) && next!(p)
  end
end

