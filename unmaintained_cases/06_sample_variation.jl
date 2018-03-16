include("main.jl")

REWRITE_DATA = true

BUSES, LINES, GENERATORS, rPTDF = load_case_data()
ERROR_VARIANCES_TRUE, ERROR_DIST_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()

# Test dependency on variance estimation confidence
casename = "sample_number_variation"
case_path = "results/$(casename)"

if REWRITE_DATA && isdir(case_path)
    rm(case_path, recursive=true)
end

if !(isdir(case_path)) || REWRITE_DATA
    # Only run if there are no saved results
    mkdir(case_path)

    testvalues = [20, 50, 100, 500, 1000, 5000]
    # var_conf = [0.1, 0.05, 0.01, 0.005, 0.001]
    var_conf = [0.05]

    sample_size = 750
    # Create global set of true observations to make tests comparable
    OBSERVATIONS = rand(ERROR_DIST_TRUE, sample_size)

    cctest_dict = DataFrame(xi= Any[], eta=Any[], v_violation=Any[], f_violation=Any[], median_costs=Any[], mean_costs=Any[], min_costs=Any[], max_costs=Any[], quantile10_cost=Any[], quantile90_cost=Any[], median_delta=Any[], mean_delta=Any[], v_violation_sample=Any[])
    ddtest_dict = DataFrame(xi = Any[], eta=Any[], v_violation=Any[], f_violation=Any[], median_costs=Any[], mean_costs=Any[], min_costs=Any[], max_costs=Any[], quantile10_cost=Any[], quantile90_cost=Any[], median_delta=Any[], mean_delta=Any[], v_violation_sample=Any[])

    solvetimes_dr = []
    
    for v in testvalues
    for vc in var_conf
        info("Runnig dr for N = $v and ξ=$vc")

        ξ = vc
        ERROR_VARIANCES_TRUE, ERROR_DIST_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity(sample_N=v)

        objective, bus_results, line_results, solvetimecc = Solve_DD_DR_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_SAMPLE)
        cctest = test_results(bus_results, line_results, objective; samples = sample_size, message="Standard CC", print_on=false)

        ddobjective, ddbus_results, ddline_results, solvetimedr = Solve_DD_DR_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_WC)
        ddtest = test_results(ddbus_results, ddline_results, ddobjective; samples = sample_size, message="DD CC", print_on=false)

        push!(cctest_dict, vcat([vc, v], cctest)) # Redundant Data Saving for easier plotting
        push!(ddtest_dict, vcat([vc, v], ddtest))
    end
    end

    writetable("$(case_path)/cctest_dict.csv", cctest_dict)
    writetable("$(case_path)/ddtest_dict.csv", ddtest_dict)

else
    # If there is saved data load it
    info("Data already exists at $case_path")
    cctest_dict = readtable("$(case_path)/cctest_dict.csv")
    ddtest_dict = readtable("$(case_path)/ddtest_dict.csv")
end
