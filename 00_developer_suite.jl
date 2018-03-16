# Use to provide proof that everything works
include("src/main.jl")

BUSES, LINES, GENERATORS, rPTDF = load_case_data(datafile="basecase")
# BUSES, LINES, GENERATORS, rPTDF = load_case_data(datafile="33buscase_pu")

η_v = .05 # 1 - Confidence for voltage limit cc
η_g = .05 # 1 - Confidence for generation limit cc
ξ = .0005# 1 - Confidence for variance estimator

# SND Quantiles for global use
snd = Normal(0,1)
z_v = quantile(snd, 1-η_v)
z_g = quantile(snd, 1-η_g)
z_g2 = quantile(snd, 1-(η_g/2))

# Set True Load Standard Deviation
ld = 0.2


ERROR_VARIANCES_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE =
    init_stochasticity(variance_opt="implicit")


To test feasability
ld_objective, ld_bus_results, ld_line_results = Solve_LinDist(BUSES, LINES, GENERATORS)
info("Results of standard lindist. Objective: $ld_objective")
println(ld_bus_results)
println(ld_line_results)

# To test SOC
# socobjective, socbus_results, socline_results = Solve_Dist_SOC(BUSES, LINES, GENERATORS)
# info("Results of soc dist. Objective: $fobjective")
# println(fbus_results)
# println(fline_results)

# linobjective, linbus_results, linline_results = Solve_Linearized_Dist_SOC(BUSES, LINES, GENERATORS,
#                         socline_results[:f_P], socline_results[:f_Q], socline_results[:a_squared])
# info("Results of soc dist. Objective: $fobjective")
# println(fbus_results)
# println(fline_results)

# println("SOC Objective = $(socobjective)")
# println("LIN Objective = $(linobjective)")
# println("LinDist Objective = $(ld_objective)")
# compare_voltages(socbus_results, linbus_results, ld_bus_results, label1="soc", label2="linearized", label3="lindist")
# compare_voltages(socbus_results, linbus_results, label1="soc", label2="linearized")
# compare_voltages(ld_bus_results, linbus_results, label1="lindist", label2="linearized")
# compare_voltages(fbus_results, b_ref, label1="soc_dist", label2="ref" )
# compare_voltages(fbus_results, ld_bus_results, label1="soc_dist", label2="lin_dist" )


ccobjective, ccbus_results, ccline_results = Solve_DD_DR_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_SAMPLE)
info("Results of standard cc. Objective: $ccobjective")
println(ccbus_results)
println(ccline_results)

objective, bus_results, line_results = Solve_DD_DR_CCLinDist(BUSES, LINES, GENERATORS, ERROR_VARIANCES_WC)
info("Results of dd cc. Objective: $objective")
println(bus_results)
println(line_results)


sample_size = 750
# Create global set of true observations to make tests comparable
OBSERVATIONS = create_data(ERROR_VARIANCES_TRUE, sample_size)
test_results(ccbus_results, ccline_results, ccobjective; samples = sample_size, message="Standard CC")
test_results(bus_results, line_results, objective; samples = sample_size, message="DD CC")
