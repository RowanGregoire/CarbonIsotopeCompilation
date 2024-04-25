## --- Definitions 

    colors = (;
        org_light = :navajowhite,
        org_dark = :darkorange,
        carb_light = :powderblue,
        carb_dark = :mediumblue,
    )


## --- Basic utilities for carbon isotope analysis 

    """
    ```julia
    f_org([in], carb, org)
    ```

    Calculate the fraction of carbon buried as organic given δ¹³C compositions for mantle 
    carbon (defaults to -5.5 ‰), carbonate / inorganic carbon, and organic carbon.
    """
    f_org(in::Number, carb::Number, org::Number) = (in - carb)/(org - carb)
    f_org(carb::Number, org::Number) = (-5.5 - carb)/(org - carb)


## --- Rayleigh fractionation functions 

    # Function which describes Rayleigh curve
    @. r₀(HC, p) = p[1]/(HC/p[2] + p[3])^(p[4]-1) + p[5]

    # Find parameters to fit data to curve
    function fit_rayleigh(d13c, hc)
        p₀ = Float64[1, 2, 1e-3, 1.3, -30]
        lb = Float64[0, 1e-3, 1e-3, 1, -50]
        ub = Float64[Inf, 10, 1, 10, 0]

        fitted = curve_fit(r₀, hc, d13c, p₀, lower=lb, upper=ub)
        return fitted.param
    end

    # Calculate Rayleigh curve
    function rayleigh_curve(p, hc)
        x = 0:0.01:maximum(hc)
        y = r₀(x, p)
        return x, y
    end


## --- Likelihood and sample matching computations 

    # Calculate likelihood of match between two samples. 
    # Return the index of the matched sample.
    function likelihood(bulkage::AbstractArray, sampleage::Number,
            bulklat::AbstractArray, bulklon::AbstractArray, 
            samplelat::Number, samplelon::Number,
            index::AbstractArray,
        )

        # Preallocate
        npoints = length(bulkage)
        ll_total = zeros(npoints, 1)

        @turbo for i in 1:npoints
            ll_total[i] -= ((bulkage[i] - sampleage)^2)/(38^2)
            ll_total[i] -= ((haversine(samplelat, samplelon, bulklat[i], bulklon[i]))^2)/(18.0^2)
        end

        matched_sample = rand_prop_liklihood(ll_total)
        return index[matched_sample]
    end

    # Weighted-random selection of an index based on log-likelihoods `ll`.
    function rand_prop_liklihood(ll)
        log_sum_ll = logsumexp(ll)
        r = rand()*exp(log_sum_ll)
        s = zero(typeof(log_sum_ll))
        @inbounds for i in eachindex(ll)
            s += exp(ll[i])
            if s > r
                return i
            end
        end
        return lastindex(ll)
    end


## --- End of file 