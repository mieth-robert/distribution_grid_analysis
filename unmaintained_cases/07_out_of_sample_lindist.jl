include("main.jl")

REWRITE_DATA = true

BUSES, LINES, GENERATORS, rPTDF = load_case_data()
ERROR_VARIANCES_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()

# Test dependency on variance estimation confidence
casename = "out_of_sample_lindist"
case_path = "results/$(casename)"

if REWRITE_DATA && isdir(case_path)
    rm(case_path, recursive=true)
end

if !(isdir(case_path)) || REWRITE_DATA
    # Only run if there are no saved results
    mkdir(case_path)

    sample_size = 750
    # Create global set of true observations to make tests comparable
    OBSERVATIONS = create_data(ERROR_VARIANCES_TRUE, sample_size)
    
    objective, bus_results, line_results = Solve_LinDist(BUSES, LINES, GENERATORS)
    det_alpha = vcat([1], zeros(length(BUSES)-1)) 
    # In Deterministic Case all imbalance is controlled by the substation
    bus_results = hcat(bus_results, DataFrame(α=det_alpha))
    ldtest = test_results(bus_results, line_results, objective; samples = sample_size, message="Standard CC", print_on=false)

    ldtest_dict = DataFrame()
    
    ldtest_dict = vcat(ldtest_dict, ldtest)  

    writetable("$(case_path)/ldtest_dict.csv", ldtest_dict)

else
    # If there is saved data load it
    info("Data already exists at $case_path")
    cctest_dict = readtable("$(case_path)/cctest_dict.csv")
    ddtest_dict = readtable("$(case_path)/ddtest_dict.csv")
end
