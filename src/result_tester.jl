

function progress(t, max_t)
    max_bars = 50
    r = t/max_t
    n = Int(floor(r*max_bars))
    p = Int(floor(r*100))
    bar = repeat("=", n)
    res = repeat(" ", max_bars - n)
    out = "[" * bar * res * "] " * string(p) * "%"
    (t == max_t ? print("\r" * out * "\n") : print("\r" * out))
    flush(STDOUT)
end

function exp_arr_quantile(a, q)
    sort!(a)
    n = length(a)
    i = max(1,floor(n*q))
    i = convert(Int, i)
    return(a[i])
end


function test_results(bus_results, line_results, exp_costs;
                        samples=100, message="", print_on=true)
# Tests the calculated results via Monte Carlo Method
if message != ""
    message = "($message)"
end
info("Testing Results $message")

n_buses = nrow(bus_results)
n_lines = nrow(line_results)

# Prepare Base Data
dP = bus_results[:d_P]
dQ = bus_results[:d_Q]
gP = bus_results[:g_P]
gQ = bus_results[:g_Q]
α = bus_results[:α]
v_base = bus_results[:v_squared]

fP_base = line_results[:f_P]
fQ_base = line_results[:f_Q]

v_violation = 0
g_violation = 0
f_violation = 0
v_violateion_in_sample = 0

# Some helping Data
gen_bus = []
for g in keys(GENERATORS)
    push!(gen_bus, GENERATORS[g].bus_idx)
end
lines_to = Dict()
for l in keys(LINES)
    lines_to[LINES[l].to_node] = LINES[l]
end
root_bus = 0
for b in keys(BUSES)
    if BUSES[b].is_root
        root_bus = b
    end
end
line_arr = [LINES[l] for l in 1:n_lines]

real_costs = zeros(samples)
for s in 1:samples
    # if print_on
    progress(s, samples)
    # end
    # Create Load Deviation from true distribution
    if s > length(OBSERVATIONS[1,:])
        warn("Requested sample size $(samples) larger than available sample set.")
        break
    end

    v_viol_flag = false

    load_dev_P = OBSERVATIONS[:,s]
    tanphi = [BUSES[b].tanphi for b in keys(BUSES)]
    load_dev_Q = load_dev_P .* tanphi
    load_dev_P_sum = sum(load_dev_P)
    load_dev_Q_sum = sum(load_dev_Q)

    # Adapt load and generation
    dP_real = dP + load_dev_P
    dQ_real = dQ + load_dev_Q

    gP_real = gP + (α .* load_dev_P_sum)
    gQ_real = gQ + (α .* load_dev_Q_sum)

    # Calculate Flow and Voltage Based on LinDistFlow
    net_load_P = dP_real - gP_real
    net_load_Q = dQ_real - gQ_real

    fP_real = rPTDF*net_load_P
    fQ_real = rPTDF*net_load_Q
    R = [l.r for l in line_arr]
    X = [l.x for l in line_arr]

    v_real = ones(n_buses) - 2 * rPTDF' * (R.*fP_real + X.*fQ_real)

    # Double Check Energy Balance
    if abs(sum(dP_real) - sum(gP_real)) > 1e-8
        warn("Active Power unbalanced at test $s")
        println(abs(sum(dP_real) - sum(gP_real)))
    end
    if abs(sum(dQ_real) - sum(gQ_real)) > 1e-8
        warn("Reactive Power unbalanced at test $s")
    end
    eb_P = zeros(n_buses)
    for b in 2:n_buses
        eb_P[b] = net_load_P[b] - fP_real[b-1]
        if length(BUSES[b].children) > 0
            eb_P[b] = eb_P[b] + sum(fP_real[c-1] for c in BUSES[b].children)
        end
    end
    eb_Q = zeros(n_buses)
    for b in 2:n_buses
        eb_Q[b] = net_load_Q[b] - fQ_real[b-1]
        if length(BUSES[b].children) > 0
            eb_Q[b] = eb_Q[b] + sum(fQ_real[c-1] for c in BUSES[b].children)
        end
    end
    for b in 1:n_buses
        if abs(eb_P[b]) > 1e-8
            warn("Active Power unbalanced at node $i in sample $s")
        end
        if abs(eb_Q[b]) > 1e-8
            warn("Reactive Power unbalanced at node $i in sample $s")
        end
    end

    # Calculate cost
    real_costs[s] = sum(BUSES[b].generator.cost * gP_real[b]^2 for b in gen_bus)

    # Check for constraint violations
    for i in gen_bus
        gP_real[i] > BUSES[i].generator.g_P_max ? g_violation += 1 : nothing
        gQ_real[i] > BUSES[i].generator.g_Q_max ? g_violation += 1 : nothing
        gP_real[i] < 0 ? g_violation += 1 : nothing
        gQ_real[i] < - BUSES[i].generator.g_Q_max ? g_violation += 1 : nothing
    end
    for i in 1:n_buses
        if v_real[i] > BUSES[i].v_max
            v_viol_flag = true
            v_violation += 1
        elseif v_real[i] < BUSES[i].v_min
            v_viol_flag = true
            v_violation += 1
        end
    end
    for i in length(LINES)
        (fP_real[i]^2 + fQ_real[i]^2) > (LINES[i].s_max)^2 ? f_violation += 1 : nothing
    end

    if v_viol_flag
        v_violateion_in_sample +=1
    end
end

#Check the costs
median_costs = median(real_costs)
mean_costs = mean(real_costs)
min_costs = minimum(real_costs)
max_costs = maximum(real_costs)
std_dev_cost = std(real_costs)
quantile10_cost = exp_arr_quantile(real_costs, 0.1)
quantile90_cost = exp_arr_quantile(real_costs, 0.9)

delta_costs = real_costs .- exp_costs
median_delta = median(delta_costs)
mean_delta = mean(delta_costs)


if print_on
println(" ++ Test Results with $samples samples ++ ")
@printf("%.2f%% of voltage constraints hold (%d violations).  \n", ((1-(v_violation/(2*samples*n_buses)))*100), v_violation)
@printf("%.2f%% of generation constraints hold (%d violations) . \n", ((1-(g_violation/(4*samples*length(gen_bus))))*100), g_violation)
@printf("%.2f%% of flow constraints hold (%d violations). \n", ((1-(f_violation/(samples*length(LINES))))*100), f_violation)
println()
@printf("Expected Cost: %.4f, Median Costs: %.4f (Min: %.4f, Max: %.4f). \n", exp_costs, median_costs, min_costs, max_costs)
@printf("Median deviation from expectation: %.4f. \n", median_delta)
println()
end

test_results = DataFrame(
                v_violation = v_violation,
                f_violation = f_violation,
                g_violation = g_violation,
                median_costs = median_costs,
                mean_costs = mean_costs,
                min_costs = min_costs,
                max_costs = max_costs,
                std_dev_cost = std_dev_cost,
                quantile10_cost = quantile10_cost,
                quantile90_cost = quantile90_cost,
                median_delta = median_delta,
                mean_delta = mean_delta,
                v_violation_sample =  v_violateion_in_sample)

return test_results
end
