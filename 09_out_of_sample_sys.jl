include("main.jl")

REWRITE_DATA = true

BUSES, LINES, GENERATORS, rPTDF = load_case_data()
ERROR_VARIANCES_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()

# Test dependency on variance estimation confidence
casename = "true_out_of_sample_systematic"
case_path = "results/$(casename)"

if REWRITE_DATA && isdir(case_path)
    rm(case_path, recursive=true)
end

if !(isdir(case_path)) || REWRITE_DATA
    # Only run if there are no saved results
    mkdir(case_path)

    testvalues = [0.05, 0.03]
    var_conf = [0.1, 0.01, 0.005, 0.0001]
    oos_size = [0.2, 0.5, 0.8]
    
    cctest_dict = DataFrame()
    ddtest_dict = DataFrame()
    ldtest_dict = DataFrame()

    for oos in oos_size
        sample_size = 750
        ξ = .001
        ERROR_VARIANCES_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()
        oos_vars = ERROR_VARIANCES_TRUE + ((ERROR_VARIANCES_WC - ERROR_VARIANCES_TRUE) * oos)
        global OBSERVATIONS = create_data(oos_vars, sample_size)

        objective, bus_results, line_results = Solve_LinDist(BUSES, LINES, GENERATORS)
        det_alpha = vcat([1], zeros(length(BUSES)-1)) 
        # In Deterministic Case all imbalance is controlled by the substation
        bus_results = hcat(bus_results, DataFrame(α=det_alpha))
        ldtest = test_results(bus_results, line_results, objective; samples = sample_size, message="Standard CC", print_on=false)

        for v in testvalues
            # Set Confidence Intervals
            η_v = v # 1 - Confidence for voltage limit cc
            η_g = v # 1 - Confidence for generation limit cc
            # SND Quantiles for global use
            snd = Normal(0,1)
            z_v = quantile(snd, 1-η_v)
            z_g = quantile(snd, 1-η_g)
            z_g2 = quantile(snd, 1-(η_g/2))

            buses_ad = change_load_same_pf(3.5)
            objective, bus_results, line_results, solvetimecc = Solve_DD_DR_CCLinDist(buses_ad, LINES, GENERATORS, ERROR_VARIANCES_SAMPLE)
            cctest = test_results(bus_results, line_results, objective; samples = sample_size, message="Standard CC", print_on=false)
            

            for vc in var_conf
                info("Runnig dr for η = $v and ξ=$vc and oos=$oos")

                ξ = vc
                ERROR_VARIANCES_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity(create_new_data=false)

                ddobjective, ddbus_results, ddline_results, solvetimedr = Solve_DD_DR_CCLinDist(buses_ad, LINES, GENERATORS, ERROR_VARIANCES_WC)
                ddtest = test_results(ddbus_results, ddline_results, ddobjective; samples = sample_size, message="DD CC", print_on=false)

                param = DataFrame(oos = oos, xi=vc, eta=v)
                cctest_dict = vcat(cctest_dict, hcat(param, cctest)) # Redundant Data Saving for easier plotting
                ddtest_dict = vcat(ddtest_dict, hcat(param, ddtest)) 

            end
        end
        ldtest_dict = vcat(ldtest_dict, hcat(DataFrame(oos=oos), ldtest))
    end

    writetable("$(case_path)/cctest_dict.csv", cctest_dict)
    writetable("$(case_path)/ddtest_dict.csv", ddtest_dict)
    writetable("$(case_path)/ldtest_dict.csv", ldtest_dict)

else
    # If there is saved data load it
    info("Data already exists at $case_path")
    cctest_dict = readtable("$(case_path)/cctest_dict.csv")
    ddtest_dict = readtable("$(case_path)/ddtest_dict.csv")
end
