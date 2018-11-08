## -- Load required packages

    using Chron
    using Plots; gr(); default(fmt=:svg);

    if VERSION>=v"0.7"
        using Statistics, DelimitedFiles, SpecialFunctions
    else
        using Compat
    end

## --- Prepare inversion
    # Bayesian intercalibration

    d = [1.9190E-03 1.00E-06  6.7890E-05  5.5000E-07  "79CE"
        17.42000    0.03000   0.60556  0.00249  "IP"
        125.65000   0.17000   4.54940  0.00610  "P1T-2"
        133.06000   0.26000   4.82050  0.01530  "PR94-7"
        201.27000   0.13000   7.46360  0.00720  "NMB"
        242.14000   0.45000   9.06190  0.01220  "MSG-09"
        252.17000   0.36000   9.49830  0.01020  "SH-10"
        252.40000   0.33000   9.49180  0.00410  "D3T"
        252.60000   0.33000   9.51020  0.00590  "SH-09"
        253.02000   0.50000   9.49580  0.01220  "PD97-2"
        253.24000   0.33000   9.52950  0.00680  "SH-27"
        253.68000   0.30000   9.53910  0.00850  "SH-16"
        257.30000   0.30000   9.69870  0.00530  "SH-08"
        327.99000   0.45000  12.60390  0.01480  "PaV"
        454.59000   0.56000  18.09000  0.05230  "DK-LQ"
        1094.20000  1.15000  52.90110  0.11620  "F-239"
        2067.50000  1.40000 135.63890  0.18210  "EGB-032"
    ]

    data = Dict()
    data["Sample"] = d[:,5]
    data["t"] = Array{Float64}(1e6*d[:,1])
    data["t_sigma"] = Array{Float64}(1e6*d[:,2])
    data["R"] = Array{Float64}(d[:,3])
    data["R_sigma"] = Array{Float64}(d[:,4])

    # Other constraints -- order: K λe λb λtot
    data["constants"]       = [1.6408E-03, 5.8000E-11, 4.8840E-10, 5.5545E-10]
    data["constants_sigma"] = [4.7000E-06, 7.0000E-13, 4.9000E-12, 1.0900E-12]

    # Paul and Greg's original residence time correction
    # data["t"][2:end] = data["t"][2:end] .- 90000
    # data["t_sigma"][2:end] = sqrt.( (data["t_sigma"][2:end]).^2 .+ 77000^2)

    # For comparison
    results=Dict()
    results["Description"] = ["Literature\nconstraints", "Renne et al.\n2010", "Balco\nno lambda t"]
    results["K"] =           [1.6408E-03, 1.6418E-03, 1.6417E-03]
    results["K_sigma"] =     [0.0047E-03, 0.0045E-03, 0.0046E-03]
    results["λe"] =          [0.5800E-10, 0.5755E-10, 0.5756E-10]
    results["λe_sigma"] =    [0.0070E-10, 0.0016E-10, 0.0017E-10]
    results["λb"] =          [4.8840E-10, 4.9737E-10, 4.9551E-10]
    results["λb_sigma"] =    [0.0490E-10, 0.0093E-10, 0.0135E-10]
    results["λ40"] =         [5.5545E-10, 5.5492E-10, 5.5308E-10]
    results["λ40_sigma"] =   [0.0109E-10, 0.0093E-10, 0.0136E-10]


    function LL(K, λe, λb, Ti, data)
        # Objective function to determine Ar decay parameters
        # data is a dict containing the t's and r's

        # Calculate combined decay constant
        λ = λb + λe

        # Compute R's
        Ri = (λe./(K.*λ)) .* (exp.(Ti.*λ) .- 1)
        # Or compute t's
        # Ri = Ti
        # Ti = (1./(λe+λb)).*log(((λe+λb)./λe).*K.*Ri + 1)

        # Match to t's and r's
        SSr = 0.5 .* ( (Ri .- data["R"]) ./ data["R_sigma"] ).^2
        SSt = 0.5 .* ( (Ti .- data["t"]) ./ data["t_sigma"] ).^2
        # SSr = ( (Ri .- data["R"]) ./ data["R_sigma"] ).^2
        # SSt = ( (Ti .- data["t"]) ./ data["t_sigma"] ).^2

        # Match to other constraints
        #SSconstants = (([K λe λb λ]' - data.a)./data.dela).^2

        # Ignore lambda total
        SSconstants = (([K, λe, λb] .- data["constants"][1:3]) ./ data["constants_sigma"][1:3]).^2

        return - sum(SSr) - sum(SSt) - sum(SSconstants)
    end

    function runMC(nsteps, sieve, data)
        # Use physical constraints as first proposal
        (K, λe, λb) = data["constants"][1:3]
        (K_sigma, λe_sigma, λb_sigma) = data["constants_sigma"][1:3]
        t = data["t"]
        t_sigma = data["t_sigma"]

        # Calculate log likelihood of first proposal
        ll = LL(K, λe, λb, t, data)

        n = 0
        ll_dist = Array{Float64}(undef, nsteps)
        K_dist = Array{Float64}(undef, nsteps)
        λe_dist = Array{Float64}(undef, nsteps)
        λb_dist = Array{Float64}(undef, nsteps)
        t_dist = Array{Float64}(undef, length(t), nsteps)
        for i=1:(nsteps*sieve)
            # Reset proposal
            K_prop = copy(K)
            λe_prop = copy(λe)
            λb_prop = copy(λb)
            t_prop = copy(t)

            # Choose which element to adjust
            r = rand(Float64)
            if r<0.25
                # Adjust K
                K_prop += randn(Float64) .* K_sigma
            elseif r<0.5
                # Adjust λe
                λe_prop += randn(Float64) .* λe_sigma
            elseif r<0.75
                # Adjust λb
                λb_prop += randn(Float64) .* λb_sigma
            else
                # Adjust t
                r_index = ceil(Int, rand(Float64)*length(t))
                t_prop[r_index] += randn(Float64) .* t_sigma[r_index]
            end

            # Calculate log likelihood of propsal
            ll_prop = LL(K_prop, λe_prop, λb_prop, t_prop, data)

            # Accept with probability p = exp(ll_prop-ll)
            if log(rand(Float64)) < (ll_prop - ll)
                # Accept proposal
                ll = copy(ll_prop)
                K = copy(K_prop)
                λe = copy(λe_prop)
                λb = copy(λb_prop)
                t = copy(t_prop)
            end

            # record current accepted proposal
            if mod(i,sieve)==0
                n += 1
                ll_dist[n] = ll
                K_dist[n] = K
                λe_dist[n] = λe
                λb_dist[n] = λb
                t_dist[:,n] = t
            end
        end
        return (ll_dist, K_dist, λe_dist, λb_dist, t_dist)
    end

    # Run Monte Carlo
    nsteps = 10^6
    burnin = 10000
    sieve = 32
    (ll_dist, K_dist, λe_dist, λb_dist, t_dist) = runMC(nsteps, sieve, data)
    λ40_dist = λe_dist + λb_dist

    push!(results["Description"], "Bayesian\nimplementation")
    push!(results["K"], mean(K_dist[burnin:end]))
    push!(results["K_sigma"], std(K_dist[burnin:end]))
    push!(results["λe"], mean(λe_dist[burnin:end]))
    push!(results["λe_sigma"], std(λe_dist[burnin:end]))
    push!(results["λb"], mean(λb_dist[burnin:end]))
    push!(results["λb_sigma"], std(λb_dist[burnin:end]))
    push!(results["λ40"], mean(λ40_dist[burnin:end]))
    push!(results["λ40_sigma"], std(λ40_dist[burnin:end]))

    h = plot(ll_dist, xlabel="step number", ylabel="log likelihood", label="")
    savefig(h,"Burnin1.png")

    t_mean_b = mean(t_dist, dims=2)
    t_sigma_b = std(t_dist, dims=2)
    hresidU = plot(data["t"] - t_mean_b, yerror=data["t_sigma"], seriestype=:scatter, label="Bayesian")
    plot!(hresidU, xticks=(1:length(data["Sample"]), data["Sample"]), xrotation=90, ylabel="U-Pb age - t_opt", legend=:bottomleft, fg_color_legend=:white)
    savefig(hresidU,"residualsU.pdf")

    λe = results["λe"][4]
    λb = results["λb"][4]
    K  = results["K"][4]
    tAr = (1 ./ (λe+λb)).*log.(((λe+λb)./λe).*K.*data["R"] .+ 1)
    tAr_p = (1 ./ (λe+λb)).*log.(((λe+λb)./λe).*K.*(data["R"]+data["R_sigma"]) .+ 1)
    hresidAr = plot(tAr - t_mean_b, yerror=abs.(tAr_p-tAr), seriestype=:scatter, label="Bayesian")
    plot!(hresidAr, xticks=(1:length(data["Sample"]), data["Sample"]), xrotation=90, ylabel="Ar-Ar age - t_opt", legend=:topleft, fg_color_legend=:white)
    savefig(hresidAr,"residualsAr.pdf")



## ---
    h = plot(results["K"]*1E3, yerror=results["K_sigma"]*1E3, seriestype=:scatter, ylabel="K_FCS *1E3", xticks=(1:length(results["Description"]), results["Description"]), label="")
    savefig(h, "Kfcs_comparison1.pdf")

    h = plot(results["λe"]*1E10, yerror=results["λe_sigma"]*1E10, seriestype=:scatter, ylabel="lambda e *1E10", xticks=(1:length(results["Description"]), results["Description"]), label="")
    savefig(h, "lambda_e_comparison1.pdf")

    h = plot(results["λb"]*1E10, yerror=results["λb_sigma"]*1E10, seriestype=:scatter, ylabel="lambda b *1E10", xticks=(1:length(results["Description"]), results["Description"]), label="")
    savefig(h, "lambda_b_comparison1.pdf")
    display(h);

    h = plot(results["λ40"]*1E10, yerror=results["λ40_sigma"]*1E10, seriestype=:scatter, ylabel="lambda K-40 *1E10", xticks=(1:length(results["Description"]), results["Description"]), label="")
    savefig(h, "lambda_40_comparison1.pdf")


## --- Uranium decay constraints / reset results

    results["Description"] = results["Description"][1:4]
    results["K"] = results["K"][1:4]
    results["K_sigma"] = results["K_sigma"][1:4]
    results["λe"] = results["λe"][1:4]
    results["λe_sigma"] = results["λe_sigma"][1:4]
    results["λb"] = results["λb"][1:4]
    results["λb_sigma"] = results["λb_sigma"][1:4]
    results["λ40"] = results["λ40"][1:4]
    results["λ40_sigma"] = results["λ40_sigma"][1:4]

    λU238Jaffey     = 746.19 * 525960 * 1000 * 238 / 6.022E23
    λU238Jaffey_sigma = 0.41 * 525960 * 1000 * 238 / 6.022E23 * 2 / 2 # One-sigma, but including "systematic factor" of 2 from Jaffey abstract

    λU238Tashi = 4.92923E-18 * 31557600
    λU238Tashi_sigma = λU238Tashi * 0.462/100 / 2 # One-sigma

    U_results=Dict()
    U_results["Description"] = ["Jaffey\n1971", "Parsons-Davis\n2018"]
    U_results["λ238"] =        [λU238Jaffey, λU238Tashi]
    U_results["λ238_sigma"] =  [λU238Jaffey_sigma, λU238Tashi_sigma]


## --- Full bayesian, propagating Jaffey 1971 λ238 uncertainty

    # Independent constraints: λ238, K, λe, λb, λtot
    data["constants"]       = [      λU238Jaffey, 1.6408E-03, 5.8000E-11, 4.8840E-10, 5.5545E-10]
    data["constants_sigma"] = [λU238Jaffey_sigma, 4.7000E-06, 7.0000E-13, 4.9000E-12, 1.0900E-12]


    function LL_full(constants, t, data)
        # Extract relevant constants
        λ238 = constants[1]
        λr = data["constants"][1] / λ238
        K  = constants[2]
        λe = constants[3]
        λb = constants[4]
        λtot = constants[5]

        # Compute R's
        R = (λe./(K.*λtot)) .* (exp.(t.*λtot) .- 1)

        # Sum squared relative error = log likelihood for Gaussian distributions
        SSr = ( (R .- data["R"]) ./ data["R_sigma"] ).^2
        SSt = ( (t .- data["t"]*λr) ./ data["t_sigma"] ).^2
        SSconstants = ((constants .- data["constants"]) ./ data["constants_sigma"]).^2

        return - sum(SSr) - sum(SSt) - sum(SSconstants)
    end

    function runMC_full(nsteps, sieve, data)
        # Use physical constraints as first proposal
        t = data["t"]
        t_sigma = data["t_sigma"]
        constants = data["constants"]
        constants_sigma = data["constants_sigma"]

        # Calculate log likelihood of first proposal
        ll = LL_full(constants, t, data)

        n = 0
        accepted = 0
        ll_dist = Array{Float64}(undef, nsteps)
        t_dist = Array{Float64}(undef, length(t), nsteps)
        constants_dist = Array{Float64}(undef, length(constants), nsteps)
        for i=1:(nsteps*sieve)

            # Reset proposal
            t_prop = copy(t)
            constants_prop = copy(constants)

            # Choose which element to adjust
            r = rand(Float64)
            if r<0.5
                # Adjust t
                rand_index = ceil(Int, rand(Float64)*length(t))
                t_prop[rand_index] += randn(Float64) .* t_sigma[rand_index] .* 2
            else
                # Adjust constants
                rand_index = ceil(Int, rand(Float64)*(length(constants)-1))
                constants_prop[rand_index] += randn(Float64) .* constants_sigma[rand_index] .* 2
                # Make sure lambda total = lambda e plus lambda b
                constants_prop[5] = constants_prop[3] + constants_prop[4]
            end

            # Calculate log likelihood of propsal
            ll_prop = LL_full(constants_prop, t_prop, data)

            # Accept with probability p = exp(ll_prop-ll)
            if log(rand(Float64)) < (ll_prop - ll)
                # Accept proposal
                accepted += 1
                ll = copy(ll_prop)
                t = copy(t_prop)
                constants = copy(constants_prop)
            end

            # record current accepted proposal
            if mod(i,sieve)==0
                n += 1
                ll_dist[n] = ll
                constants_dist[:,n] = constants
                t_dist[:,n] = t
            end
        end
        return (ll_dist, constants_dist, t_dist, accepted)
    end


    (ll_dist, constants_dist, t_dist, accepted) = runMC_full(nsteps, sieve, data)

    constants_mean = mean(constants_dist, dims=2)
    constants_sigma = std(constants_dist, dims=2)

    push!(results["Description"], "Full bayesian\n Jaffey")
    push!(results["K"], constants_mean[2])
    push!(results["K_sigma"], constants_sigma[2])
    push!(results["λe"], constants_mean[3])
    push!(results["λe_sigma"], constants_sigma[3])
    push!(results["λb"], constants_mean[4])
    push!(results["λb_sigma"], constants_sigma[4])
    push!(results["λ40"], constants_mean[5])
    push!(results["λ40_sigma"], constants_sigma[5])

    push!(U_results["Description"], "Intercal\nw/Jaffey")
    push!(U_results["λ238"], constants_mean[1])
    push!(U_results["λ238_sigma"], constants_sigma[1])


    h = plot(ll_dist, xlabel="step number", ylabel="log likelihood", label="")
    savefig(h,"Burnin2.png")

    t_mean_J = mean(t_dist, dims=2)
    t_sigma_J = std(t_dist, dims=2)
    plot!(hresidU, data["t"]*λU238Jaffey/U_results["λ238"][3] - t_mean_J, yerror=data["t_sigma"], seriestype=:scatter, label="Bayesian + Jaffey")
    savefig(hresidU,"residualsU.pdf")

    λe = results["λe"][5]
    λb = results["λb"][5]
    K  = results["K"][5]
    tAr = (1 ./ (λe+λb)).*log.(((λe+λb)./λe).*K.*data["R"] .+ 1)
    tAr_p = (1 ./ (λe+λb)).*log.(((λe+λb)./λe).*K.*(data["R"]+data["R_sigma"]) .+ 1)
    plot!(hresidAr, tAr - t_mean_J, yerror=abs.(tAr_p-tAr), seriestype=:scatter, label="Bayesian + Jaffey")
    savefig(hresidAr,"residualsAr.pdf")

## --- Full bayesian, propagating Parsons-Davis 2018 λ238 uncertainty

    # Independent constraints: λ238, K, λe, λb, λtot
    data["constants"]       = [      λU238Tashi, 1.6408E-03, 5.8000E-11, 4.8840E-10, 5.5545E-10]
    data["constants_sigma"] = [λU238Tashi_sigma, 4.7000E-06, 7.0000E-13, 4.9000E-12, 1.0900E-12]

    (ll_dist, constants_dist, t_dist, accepted) = runMC_full(nsteps, sieve, data)

    constants_mean = mean(constants_dist, dims=2)
    constants_sigma = std(constants_dist, dims=2)

    push!(results["Description"], "Full bayesian\n Parsons-Davis")
    push!(results["K"], constants_mean[2])
    push!(results["K_sigma"], constants_sigma[2])
    push!(results["λe"], constants_mean[3])
    push!(results["λe_sigma"], constants_sigma[3])
    push!(results["λb"], constants_mean[4])
    push!(results["λb_sigma"], constants_sigma[4])
    push!(results["λ40"], constants_mean[5])
    push!(results["λ40_sigma"], constants_sigma[5])

    push!(U_results["Description"], "Intercal\nw/Parsons-Davis")
    push!(U_results["λ238"], constants_mean[1])
    push!(U_results["λ238_sigma"], constants_sigma[1])

    h = plot(ll_dist, xlabel="step number", ylabel="log likelihood", label="")
    savefig(h,"Burnin3.png")

    t_mean_PD = mean(t_dist, dims=2)
    t_sigma_PD = std(t_dist, dims=2)
    plot!(hresidU, data["t"]*λU238Jaffey/U_results["λ238"][4] - t_mean_PD, yerror=data["t_sigma"], seriestype=:scatter, label="Bayesian + Parsons-Davis")
    savefig(hresidU,"residualsU.pdf")

    λe = results["λe"][6]
    λb = results["λb"][6]
    K  = results["K"][6]
    tAr = (1 ./ (λe+λb)).*log.(((λe+λb)./λe).*K.*data["R"] .+ 1)
    tAr_p = (1 ./ (λe+λb)).*log.(((λe+λb)./λe).*K.*(data["R"]+data["R_sigma"]) .+ 1)
    plot!(hresidAr, tAr - t_mean_PD, yerror=abs.(tAr_p - tAr), seriestype=:scatter, label="Bayesian + Parsons-Davis")
    savefig(hresidAr,"residualsAr.pdf")

## ---

    h = plot(results["K"]*1E3, yerror=results["K_sigma"]*1E3, seriestype=:scatter, ylabel="K_FCS *1E3", xticks=(1:length(results["Description"]), results["Description"]), label="")
    savefig(h, "Kfcs_comparison2.pdf")

    h = plot(results["λe"]*1E10, yerror=results["λe_sigma"]*1E10, seriestype=:scatter, ylabel="lambda e *1E10", xticks=(1:length(results["Description"]), results["Description"]), label="")
    savefig(h, "lambda_e_comparison2.pdf")

    h = plot(results["λb"]*1E10, yerror=results["λb_sigma"]*1E10, seriestype=:scatter, ylabel="lambda b *1E10", xticks=(1:length(results["Description"]), results["Description"]), label="")
    savefig(h, "lambda_b_comparison2.pdf")

    h = plot(results["λ40"]*1E10, yerror=results["λ40_sigma"]*1E10, seriestype=:scatter, ylabel="lambda K-40 *1E10", xticks=(1:length(results["Description"]), results["Description"]), label="")
    savefig(h, "lambda_40_comparison2.pdf")

    h = plot(U_results["λ238"]*1E10, yerror=U_results["λ238_sigma"]*1E10, seriestype=:scatter, ylabel="lambda U-238 *1E10", xticks=(1:length(U_results["Description"]), U_results["Description"]), label="")
    savefig(h, "lambda_238_comparison.pdf")
    display(h);


## --- Chron it:
# # Sample information
#
#     # Basic information about the data
#     Path = "U:Pb 2-sigma/" # Where are the data files?
#     inputSigmaLevel = 2; # i.e., are the data files 1-sigma or 2-sigma. Integer.
#     AgeUnit = "Ma" # Unit of measurement for ages and errors in the data files
#
#     # Get filenames
#     system("ls '$Path' | grep .csv | sed -e 's/.csv//' > filenames")
#     Name = tuple(readdlm("filenames")...)
#
#     # Make an instance of a Chron StratAgeData object with the data we've just enterd
#     smpl = NewStratAgeData(Name,Path,inputSigmaLevel)
#     smpl.Age_Sidedness[:] = zeros(length(Name)) # Sidedness (zeros by default, geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
#
# # (Optional) Calculate bootstrapped distribution
#
#     # Bootstrap a KDE of the pre-eruptive (or pre-deposition) zircon distribution
#     # shape from individual sample datafiles using a KDE of stacked sample data
#     BootstrappedDistribution = BootstrapCrystDistributionKDEfromStrat(smpl)
#     h = plot(BootstrappedDistribution, xlabel="Time (arbitrary units)", ylabel="Probability Density", label="Bootstrapped distribution", fg_color_legend=:white)
#     savefig(h, joinpath(Path,"BootstrappedDistribution.pdf"))
#
# # Configure Distribution Model
#     # Number of steps to run in distribution MCMC
#     distSteps = 10^6;
#     distBurnin = floor(Int,distSteps/100);
#
#     # Choose the form of the prior distribution to use
#     # A variety of potentially useful distributions are provided in DistMetropolis.jl
#     # Options include UniformDisribution, TriangularDistribution,
#     # BootstrappedDistribution, and MeltsVolcanicZirconDistribution
#     # or you can define your own.
#     dist = TriangularDistribution;
#
#     # Run MCMC to estimate saturation and eruption/deposition age distributions
#     smpl = tMinDistMetropolis(smpl,distSteps,distBurnin,dist);
#
#     # # (Optional) Save the sample struct for later use
#     # using JLD: @save, @load
#     # @save "smpl.jld" smpl
#
#     # Print results to file
#     results = vcat(["Sample" "Age" "2.5% CI" "97.5% CI" "sigma"], hcat(collect(smpl.Name),smpl.Age,smpl.Age_025CI,smpl.Age_975CI,smpl.Age_Sigma))
#     writedlm("ChronResults.csv", results, ',')
