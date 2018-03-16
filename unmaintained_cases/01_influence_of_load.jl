include("main.jl")





BUSES, LINES, GENERATORS, rPTDF = load_case_data()
ERROR_VARIANCES_TRUE, ERROR_DIST_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()



# objective, bus_results, line_results = Solve_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_SAMPLE)
# ddobjective, bus_results, line_results = Solve_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_WC)

buses_lc = change_load_same_pf(1)


objective, bus_results, line_results = Solve_DD_DR_CCLinDist(buses_lc, LINES, GENERATORS, ERROR_VARIANCES_SAMPLE)
ddobjective, ddbus_results, ddline_results = Solve_DD_DR_CCLinDist(buses_lc, LINES, GENERATORS, ERROR_VARIANCES_WC)

info("Results of standard cc. Objective: $objective")
println(bus_results)
println(line_results)

info("Results of dd cc. Objective: $ddobjective")
println(ddbus_results)
println(ddline_results)

sample_size = 1000 # max = 10000
# Create global set of true observations to make tests comparable
OBSERVATIONS = rand(ERROR_DIST_TRUE, sample_size)
test_results(bus_results, line_results, objective; samples = sample_size, message="Standard CC")
test_results(ddbus_results, ddline_results, ddobjective; samples = sample_size, message="DD CC")