using Distributions

# Set Confidence Intervals
η_v = .05 # 1 - Confidence for voltage limit cc
η_g = .05 # 1 - Confidence for generation limit cc
ξ = .005# 1 - Confidence for variance estimator

# SND Quantiles for global use
snd = Normal(0,1)
z_v = quantile(snd, 1-η_v)
z_g = quantile(snd, 1-η_g)
z_g2 = quantile(snd, 1-(η_g/2))

# Set True Load Standard Deviation
ld = 0.2


# Save Historic Data Globally 
global ERROR_HIST = 0

# Calculate variance estimation confidence
function chisq_interval(s, n, alpha)
    # println("s = $s, n = $n, alpha = $alpha")
    # s... sample variance
    # n... degrees of freedom
    # alpha... confidence
    # Returns lower and upper bound of the 1-alpha-confidence interval
    # of the estmated variance
    ch = Chisq(n)
    lower = (n * s)/quantile(ch, 1- alpha)
    upper = (n * s)/quantile(ch, alpha)
    return [lower, upper]
end

function create_data(var_true, N)
    obs = zeros(length(var_true), N)
    for (i,v) in enumerate(var_true)
        if v > 0
            dist = Normal(0,sqrt(v))
            obs[i,:] = rand(dist, N)
        end
    end
    return obs
end

function init_stochasticity(;create_new_data=true, sample_N=100, variance_opt="implicit")

# Set the ture error distribution
# Same Variance for each bus
n_buses = length(BUSES)
if variance_opt == "implicit"
# Create individual node variance from load 
    loads = [BUSES[b].d_P for b in 1:n_buses]
    var_true_vector = (loads.*ld).^2
    # Variance = (ld% of load)^2
elseif variance_opt == "explicit"
# Same variance for all nodes
    if DATA_DIR == "basecase"
        f = 0.01
    elseif DATA_DIR == "simplecase"
        f = 0.1
    elseif DATA_DIR == "33buscase_pu"
        f = 4
    else
        f = 1
        warn("Data Dir. $(DATA_DIR) unkown for uncertainty setting")
    end
    var_true_vector = ones(n_buses)*f
elseif variance_opt == "manual"
# Manually create an array of variances
    var_true_vector = []
    if length(var_true_vector) < n_buses
        warn("Variance Vector to short")
        var_true_vector = ones(n_buses).*1e-8
    end
else
    warn("Unknown Variance Setting")
    var_true_vector = ones(n_buses).*1e-8
end

# Create or load data for samples
if create_new_data || ERROR_HIST == 0 
    info("Creating new Historical Observations for DR Implementation")
    N = sample_N
    error_hist = create_data(var_true_vector, N)
    global ERROR_HIST = error_hist
else
    N = length(ERROR_HIST[1,:])
    error_hist = ERROR_HIST
end


# Work with the data
# Calculate sample variances
var_sample_vector =[]
for i in 1:n_buses
    sn = 1/N * sum(error_hist[i,:].^2)
    push!(var_sample_vector, sn)
end

# var_sample_covmat = 1/N * error_hist * error_hist'
# var_sample_vector =  diag(var_sample_covmat) # Drop covariance entries

intervals = [chisq_interval(s_n, N, ξ) for s_n in var_sample_vector]
var_sample_lower = [intervals[i][1] for i in 1:n_buses]
var_sample_upper = [intervals[i][2] for i in 1:n_buses]

# return var_true_vector, error_dist_true, var_sample_upper, var_sample_vector
return var_true_vector, var_sample_upper, var_sample_vector
end