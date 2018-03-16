# DATA_DIR = "simplecase"
DATA_DIR = "basecase"
# DATA_DIR = "basecase_lv" # same as basecase but with lower voltage limits
# DATA_DIR = "33buscase_pu"

include("data_handler.jl")

# Provides functions to read and manipulate case data
# Creates BUSES, LINES, GENERATORS as Dict
include("stochasticity.jl")
# Settings and inital calculations for load stochasticity
include("lindist.jl")
# Provides function Solve_LinDist(buses, lines, generators) for feasability testing
include("standard_cc.jl")
# Common CC OPF Formulation in JuMPChance
# provides function Solve_CCLinDist(buses, lines, generators, error_variances)
include("dd_dist_robust_cc.jl")
# provides function Solve_CCLinDist(buses, lines, generators, error_variances)
include("result_tester.jl")
# Provides function test_results(bus_results, line_results, samples)
include("soc_flow.jl")
# Provides function Solve_Dist_SOC(buses, lines, generators) with a SOC dist fow implementation
include("lin_soc_flow.jl")
# Provides the function function Solve_Linearized_Dist_SOC(buses, lines, generators, fp0, fq0, a0)
include("result_plotters.jl")
# Provedes some functiuons to vizualize results
