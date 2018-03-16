include("src/main.jl")

REWRITE_DATA = true

BUSES, LINES, GENERATORS, rPTDF = load_case_data()
ERROR_VARIANCES_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()

# Test dependency on variance estimation confidence
casename = "risk_level_and_varconf_variation"
case_path = "results/$(casename)"

if REWRITE_DATA && isdir(case_path)
    rm(case_path, recursive=true)
end

if !(isdir(case_path)) || REWRITE_DATA
    # Only run if there are no saved results
    mkdir(case_path)

    # testvalues = [0.15, 0.12, 0.1, 0.08, 0.05, 0.03]
    testvalues = [0.1, 0.08, 0.05, 0.03, 0.01, 0.005]
    # var_conf = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
    var_conf = [0.1, 0.01, 0.005, 0.0001]

    sample_size = 750
    # Create global set of true observations to make tests comparable
    OBSERVATIONS = create_data(ERROR_VARIANCES_TRUE, sample_size)

    cctest_dict = DataFrame()
    ddtest_dict = DataFrame()
    ldtest_dict = DataFrame()
    
    buses_ad = change_load_same_pf(1)

    objective, bus_results, line_results = Solve_LinDist(BUSES, LINES, GENERATORS)
    det_alpha = vcat([1], zeros(length(BUSES)-1)) 
    # In Deterministic Case all imbalance is controlled by the substation
    bus_results = hcat(bus_results, DataFrame(α=det_alpha))
    ldtest = test_results(bus_results, line_results, objective; samples = sample_size, message="Standard CC", print_on=false)
    ldtest_dict = vcat(ldtest_dict, ldtest)


    for v in testvalues
        # Set Confidence Intervals
        η_v = v # 1 - Confidence for voltage limit cc
        η_g = v # 1 - Confidence for generation limit cc
        # SND Quantiles for global use
        snd = Normal(0,1)
        z_v = quantile(snd, 1-η_v)
        z_g = quantile(snd, 1-η_g)
        z_g2 = quantile(snd, 1-(η_g/2))

        objective, bus_results, line_results, solvetimecc = Solve_DD_DR_CCLinDist(buses_ad, LINES, GENERATORS, ERROR_VARIANCES_SAMPLE)
        cctest = test_results(bus_results, line_results, objective; samples = sample_size, message="Standard CC", print_on=false)
        

        for vc in var_conf
            info("Runnig dr for η = $v and ξ=$vc")

            ξ = vc
            ERROR_VARIANCES_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity(create_new_data=false)

            ddobjective, ddbus_results, ddline_results, solvetimedr = Solve_DD_DR_CCLinDist(buses_ad, LINES, GENERATORS, ERROR_VARIANCES_WC)
            ddtest = test_results(ddbus_results, ddline_results, ddobjective; samples = sample_size, message="DD CC", print_on=false)

            param = DataFrame(xi=vc, eta=v)
            cctest_dict = vcat(cctest_dict, hcat(param, cctest)) # Redundant Data Saving for easier plotting
            ddtest_dict = vcat(ddtest_dict, hcat(param, ddtest)) 

        end
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
