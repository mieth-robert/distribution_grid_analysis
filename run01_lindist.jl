include("src/main.jl")

BUSES, LINES, GENERATORS, rPTDF =  load_case_data(datafile="basecase")
# basecas = 15 Bus example from A. Papavasiliou, 2016, "Analysis of Distribution Locational Marginal Prices"

ld_objective, ld_bus_results, ld_line_results = Solve_LinDist(BUSES, LINES, GENERATORS)
info("Results of standard lindist. Objective: $ld_objective")
println(ld_bus_results)
println(ld_line_results)