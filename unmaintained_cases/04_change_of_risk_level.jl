include("main.jl")

REWRITE_DATA = true

BUSES, LINES, GENERATORS, rPTDF = load_case_data()
ERROR_VARIANCES_TRUE, ERROR_DIST_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()

# Test dependency on variance estimation confidence
casename = "risk_level_variation"
case_path = "results/$(casename)"

if REWRITE_DATA && isdir(case_path)
    rm(case_path, recursive=true)
end

if !(isdir(case_path)) || REWRITE_DATA
    # Only run if there are no saved results
    mkdir(case_path)

    testvalues = [0.15, 0.12, 0.1, 0.08, 0.05]

    sample_size = 750
    # Create global set of true observations to make tests comparable
    OBSERVATIONS = rand(ERROR_DIST_TRUE, sample_size)

    cctest_dict = DataFrame(eta=Any[], v_violation=Any[], f_violation=Any[], median_costs=Any[], mean_costs=Any[], min_costs=Any[], max_costs=Any[], quantile10_cost=Any[], quantile90_cost=Any[], median_delta=Any[], mean_delta=Any[])
    ddtest_dict = DataFrame(eta=Any[], v_violation=Any[], f_violation=Any[], median_costs=Any[], mean_costs=Any[], min_costs=Any[], max_costs=Any[], quantile10_cost=Any[], quantile90_cost=Any[], median_delta=Any[], mean_delta=Any[])

    for v in testvalues
        # Set Confidence Intervals
        η_v = v # 1 - Confidence for voltage limit cc
        η_g = v # 1 - Confidence for generation limit cc
        # SND Quantiles for global use
        snd = Normal(0,1)
        z_v = quantile(snd, 1-η_v)
        z_g = quantile(snd, 1-η_g)
        z_g2 = quantile(snd, 1-(η_g/2))

        info("Runnig for η = $v")
        ξ = 0.005
        ERROR_VARIANCES_TRUE, ERROR_DIST_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()

        objective, bus_results, line_results, solvetimecc = Solve_DD_DR_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_SAMPLE)
        cctest = test_results(bus_results, line_results, objective; samples = sample_size, message="Standard CC", print_on=false)
        push!(cctest_dict, vcat(v, cctest))

        ddobjective, ddbus_results, ddline_results, solvetimedr = Solve_DD_DR_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_WC)
        ddtest = test_results(ddbus_results, ddline_results, ddobjective; samples = sample_size, message="DD CC", print_on=false)

        push!(ddtest_dict, vcat(v, ddtest))
    end

    writetable("$(case_path)/cctest_dict.csv", cctest_dict)
    writetable("$(case_path)/ddtest_dict.csv", ddtest_dict)

else
    # If there is saved data load it
    info("Data already exists at $case_path")
    cctest_dict = readtable("$(case_path)/cctest_dict.csv")
    ddtest_dict = readtable("$(case_path)/ddtest_dict.csv")
end
