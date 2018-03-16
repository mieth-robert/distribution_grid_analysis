include("main.jl")

REWRITE_DATA = false

BUSES, LINES, GENERATORS, rPTDF = load_case_data()
ERROR_VARIANCES_TRUE, ERROR_DIST_TRUE, ERROR_VARIANCES_WC, ERROR_VARIANCES_SAMPLE = init_stochasticity()


# CASE: Numerical Experiments on the difference between SOC and LinDist
# 1) Voltage Deviations
# 2) Sum of losses and respective change in objective
# => With different demands

casename = "lindist_vs_soc"
case_path = "results/$(casename)"

if REWRITE_DATA && isdir(case_path)
    rm(case_path, recursive=true)
end

if !(isdir(case_path)) || REWRITE_DATA
    # Only run if there are no saved results
    mkdir(case_path)

    demand_cases = collect(0.5:0.1:1.3)

    voltage_result_soc = DataFrame()
    voltage_result_ld = DataFrame()
    voltage_deltas = DataFrame()
    objective_result = DataFrame(case=["SOC", "LinDist", "delta"])

    for α in demand_cases 
        buses_part_load = change_load_same_pf(α)

        soc_objective, soc_bus_results, soc_line_results = Solve_Dist_SOC(buses_part_load, LINES, GENERATORS)
        ld_objective, ld_bus_results, ld_line_results = Solve_LinDist(buses_part_load, LINES, GENERATORS)

        voltage_result_soc[Symbol(α)] = sqrt.(soc_bus_results[:v_squared])
        voltage_result_ld[Symbol(α)] = sqrt.(ld_bus_results[:v_squared])
        voltage_deltas[Symbol(α)] = (voltage_result_ld[Symbol(α)] -  voltage_result_soc[Symbol(α)]) ./ voltage_result_soc[Symbol(α)]

        objective_result[Symbol(α)] = [soc_objective, ld_objective, (ld_objective - soc_objective) / soc_objective]
    end

    bus_summary = DataFrame(bus=Any[], quant_10=Any[], quant_90=Any[], median=Any[])
    for i in 1:nrow(voltage_deltas)
        arr = [voltage_deltas[i,a] for a in names(voltage_deltas)]
        b_sum = [i, exp_arr_quantile(arr, 0.1), exp_arr_quantile(arr, 0.9), median(arr)]
        push!(bus_summary, b_sum)
    end

    writetable("$(case_path)/voltage_result_soc.csv", voltage_result_soc)
    writetable("$(case_path)/voltage_result_ld.csv", voltage_result_ld)
    writetable("$(case_path)/voltage_deltas.csv", voltage_deltas)
    writetable("$(case_path)/bus_summary.csv", bus_summary)
    writetable("$(case_path)/objective_result.csv", objective_result)

end