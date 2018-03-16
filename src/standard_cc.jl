using JuMP, JuMPChance
using Ipopt
using Gurobi
using Distributions


function Solve_CCLinDist(buses, lines, generators, error_variances)
info("Starting Chance Constraint LinDist Model")

bus_set = collect(keys(buses))
line_set = collect(keys(lines))
gen_set = collect(keys(generators))
gen_bus = []
for g in keys(generators)
    push!(gen_bus, generators[g].bus_idx)
end
root_bus = 0
for b in keys(buses)
    if buses[b].is_root
        root_bus = b
    end
end

lines_to = Dict()
for l in keys(lines)
    lines_to[lines[l].to_node] = lines[l]
end

v_root = 1

error_var_sum = sum(error_variances)

# m = ChanceModel(solver=IpoptSolver())
m = ChanceModel(solver=GurobiSolver())

@variable(m, v[bus_set] >= 0) #variable for voltage square
@variable(m, fp[bus_set]) #variable for active power flow
@variable(m, fq[bus_set]) # variable for reactive power flow
@variable(m, gp[bus_set]) #variable for active power generation
@variable(m, gq[bus_set]) #variable for reactive power generation
@variable(m, α[bus_set] >= 0) # variable to distribute load deviations (affine frequency control)

@indepnormal(m, e[1:length(bus_set)], mean=zeros(length(bus_set)), var=error_variances)
# @indepnormal(m, e2, mean=0, var=1)

@objective(m, Min, sum((gp[g]^2 + α[g]^2*error_var_sum) * buses[g].generator.cost for g in gen_bus))
# Expected Costs

# Control variable
@constraint(m, sum(α) == 1)

#Voltage and flow constraints for root node
@constraint(m, v[root_bus] == v_root)
@constraint(m, fp[root_bus] == 0)
@constraint(m, fq[root_bus] == 0)

for b in setdiff(bus_set, gen_bus)
# All buses without generation
    @constraint(m, α[b] == 0)
    @constraint(m, gp[b] == 0)
    @constraint(m, gq[b] == 0)
end

for b in bus_set
# All buses
    @constraint(m, buses[b].d_P - gp[b] + sum(fp[k] for k in buses[b].children) == fp[b])
    @constraint(m, buses[b].d_Q - gq[b] + sum(fq[k] for k in buses[b].children) == fq[b])
end

for b in setdiff(bus_set, [root_bus])
# All buses without root
    b_ancestor = buses[b].ancestor[1]
    @constraint(m, v[b] == v[b_ancestor] - 2*(lines_to[b].r * fp[b] + lines_to[b].x * fq[b]))
    @constraint(m, fp[b]^2 + fq[b]^2 <= lines_to[b].s_max^2)
end

# CHANCE CONSTRAINTS:
for b in gen_bus
# All buses with generation
    @constraint(m, gp[b] + α[b]*sum(e[bb] for bb in bus_set) <= buses[b].generator.g_P_max, with_probability = (1 - η_g))
    @constraint(m, gp[b] + α[b]*sum(e[bb] for bb in bus_set) >= 0, with_probability = (1 - η_g))

    @constraint(m, gq[b] + α[b]*sum(e[bb] * buses[bb].tanphi for bb in bus_set) <=  buses[b].generator.g_Q_max, with_probability = (1 - η_g))
    @constraint(m, gq[b] + α[b]*sum(e[bb] * buses[bb].tanphi for bb in bus_set) >= -buses[b].generator.g_Q_max, with_probability = (1 - η_g))
end

for b in setdiff(bus_set, [root_bus])
# All buses without root
    l = lines_to[b]
    line_arr = [lines[l] for l in line_set]
    @constraint(m, v[b] - 2 * (sum(rPTDF[l.index, b] * l.r * sum(rPTDF[l.index, j] * (e[j] - α[j]*sum(e[z] for z in bus_set)) for j in bus_set) for l in line_arr) +
                    sum(rPTDF[l.index, b] * l.x * sum(rPTDF[l.index, j] * (e[j]*buses[j].tanphi - α[j]*sum(e[z]*buses[j].tanphi for z in bus_set)) for j in bus_set) for l in line_arr))
                    <= buses[b].v_max, with_probability = (1 - η_v))
    @constraint(m, v[b] - 2 * (sum(rPTDF[l.index, b] * l.r * sum(rPTDF[l.index, j] * (e[j] - α[j]*sum(e[z] for z in bus_set)) for j in bus_set) for l in line_arr) +
                    sum(rPTDF[l.index, b] * l.x * sum(rPTDF[l.index, j] * (e[j]*buses[j].tanphi - α[j]*sum(e[z]*buses[j].tanphi for z in bus_set)) for j in bus_set) for l in line_arr))
                    >= buses[b].v_min, with_probability = (1 - η_v))
end

method = :Cuts
# method = :Reformulate
status = solve(m, method=method)
objective = getobjectivevalue(m)
if status == :Optimal
    info("CC solved to optimality with method $(method)")
    info("Objective: $(objective)")
else
    warn("CC not optimal, termination status $(status)")
end

# Prepare Results
bus_results = DataFrame(bus = Any[], d_P = Any[], d_Q = Any[], cosphi=Any[], tanphi = [],  g_P = Any[], g_Q = Any[], g_cost=Any[], v_squared = Any[], α=Any[])
line_results = DataFrame(line = Any[], from = Any[], to = Any[], f_P = Any[], f_Q = Any[], a = Any[])
for b in bus_set
    cost = 0
    if in(b, gen_bus)
        cost = buses[b].generator.cost
    end

    # Numerical smoothing
    α_res = getvalue(α[b])
    α_res < 1e-10 ? α_res = 0 : nothing

    res = [b, buses[b].d_P, buses[b].d_Q, buses[b].cosphi, buses[b].tanphi,  getvalue(gp[b]), getvalue(gq[b]), cost, getvalue(v[b]), α_res]
    push!(bus_results, res)

    # Calc square current on line
    a_res = 0
    fp_res = getvalue(fp[b])
    fq_res = getvalue(fq[b])
    if b != root_bus
        v_res = getvalue(v[buses[b].ancestor[1]])
        a_res = (fp_res^2 + fq_res^2)/v_res^2
        lres = [lines_to[b].index, buses[b].ancestor[1], b, fp_res, fq_res, a_res]
        push!(line_results, lres)
    end
end

sort!(bus_results, cols=:bus)
sort!(line_results, cols=:to)

return objective, bus_results, line_results
end #Solve_CCLinDist
