include("main.jl")

REWRITE_DATA = false

BUSES, LINES, GENERATORS, rPTDF = load_case_data()
ERROR_VARIANCES_TRUE, ERROR_DIST_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()

# Test dependency on variance estimation confidence
casename = "var_estimation_sensitivity"
case_path = "results/$(casename)"

if REWRITE_DATA
    rm(case_path, recursive=true)
end

if !(isdir(case_path)) || REWRITE_DATA
    # Only run if there are no saved results
    mkdir(case_path)

    testvalues = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005]

    sample_size = 750
    # Create global set of true observations to make tests comparable
    OBSERVATIONS = rand(ERROR_DIST_TRUE, sample_size)

    cctest_dict = DataFrame(v_violation=Any[], f_violation=Any[], median_costs=Any[], min_costs=Any[], max_costs=Any[], quantile10_cost=Any[], quantile90_cost=Any[], median_delta=Any[])
    ddtest_dict = DataFrame(xi=Any[], v_violation=Any[], f_violation=Any[], median_costs=Any[], min_costs=Any[], max_costs=Any[], quantile10_cost=Any[], quantile90_cost=Any[], median_delta=Any[])

    # No parameter change for normal cc
    objective, bus_results, line_results, solvetimecc = Solve_DD_DR_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_SAMPLE)
    cctest = test_results(bus_results, line_results, objective; samples = sample_size, message="Standard CC", print_on=false)
    push!(cctest_dict, cctest)

    solvetimes_dr = []
    for v in testvalues
        info("Runnig for ξ = $v")
        ξ = v
        ERROR_VARIANCES_TRUE, ERROR_DIST_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()

        ddobjective, ddbus_results, ddline_results, solvetimedr = Solve_DD_DR_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_WC)
        ddtest = test_results(ddbus_results, ddline_results, ddobjective; samples = sample_size, message="DD CC", print_on=false)

        push!(solvetimes_dr, solvetimedr)
        push!(ddtest_dict, vcat(v, ddtest))
    end

    case_arr = vcat("cc", testvalues)
    time_arr = vcat(solvetimecc, solvetimes_dr)
    solvetime_dict = DataFrame(case=case_arr, time=time_arr)

    writetable("$(case_path)/cctest_dict.csv", cctest_dict)
    writetable("$(case_path)/ddtest_dict.csv", ddtest_dict)

else
    # If there is saved data load it
    info("Data already exists at $case_path")
    cctest_dict = readtable("$(case_path)/cctest_dict.csv")
    ddtest_dict = readtable("$(case_path)/ddtest_dict.csv")
end
