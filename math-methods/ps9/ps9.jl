using StatsBase
using Plots
gr()

type MarkovChain
    P::Array{Float64, 2}
    pi::Array{Float64}
end

MarkovChain(P::Array{Float64, 2}) = MarkovChain(P, Array{Float64, 1}())
find_initial_state(mc::MarkovChain) = indmax(mc.pi) #this should be generalized

function find_stationary_distribution(mc::MarkovChain)
    val, vec = eig(transpose(mc.P))
    b = [(val[i] - 1) * (val[i] - 1) for i=1:size(val)[1]];
    (min_val, ind) = findmin(real(b));
    stationary_dist = [Real(vec[i, ind]) for  i=1:size(vec)[1]]
    stationary_dist = stationary_dist / sum(stationary_dist)
    return stationary_dist
end

function compute_bernoulli_lapace_matrix(r::Int64)
    a = zeros(r+1, r+1);
    a[1, 2] = 1.0;
    a[r+1, r] = 1.0;
    for i = 2:r
        j = i - 1; # because our indexing starts at 1, not 0
        @inbounds a[i, i-1] = (j / r)^2;
        @inbounds a[i, i+1] = ((r - j) / r)^2;
        @inbounds a[i, i] = 2 * j * (r - j) / r^2;    
    end
    return a
end

get_states(a::Array{Float64, 1}) = 1:size(a)[1]
function pick_next_state(probs::Array{Float64, 1})
    wn = WeightVec(probs)
    return sample(get_states(probs), wn)
end

is_high_return_state(current_state::Int64) = current_state == 1 || current_state == 3;
is_low_return_state(current_state::Int64) = current_state == 2 || current_state == 4;
is_high_liquidity_state(current_state::Int64) = current_state == 1 || current_state == 2;
is_low_liquidity_state(current_state::Int64) = current_state == 3 || current_state == 4;
high_return_result(cur::Float64) = cur * 0.06 + cur;
low_return_result(cur::Float64) = cur * -0.20 + cur;

function iterate_markov_chain(mc::MarkovChain, num_iter::Int64)
    current_state = find_initial_state(mc);
    state_history = [current_state];
    for i = 1:num_iter
        current_state = pick_next_state(mc.P[current_state,:])
        append!(state_history, copy(current_state))
    end
    return state_history
end

m = MarkovChain([0 1 0 0; 0 0 0.5 0.5; 1 0 0 0; 0 1 0 0]);
one_a_c = find_stationary_distribution(m);
a = MarkovChain(compute_bernoulli_lapace_matrix(6));
b = find_stationary_distribution(a);
c = MarkovChain(compute_bernoulli_lapace_matrix(7));
d = find_stationary_distribution(c);

sick_no_coffee = [0.99 0.01; 0.3 0.7];
sick_with_coffee = [0.995 0.005; 0.14 0.82];
two_weeks_transition_no_coffee = sick_no_coffee^14;
two_weeks_transition_coffee = sick_with_coffee^14;
long_run_no_coffee = find_stationary_distribution(MarkovChain(sick_no_coffee));
long_run_coffee = find_stationary_distribution(MarkovChain(sick_with_coffee));


function run_analysis_for_four(mc::MarkovChain, num_simulations::Int64, T::Int64)
    results = simulate(mc, num_simulations, T);
    return analyze_results(results, num_simulations, T);
end

function simulate(mc::MarkovChain, num_simulations::Int64, T::Int64)
    results = []
    for i=1:num_simulations
        result = iterate_markov_chain(mc, T);
        append!(results, [result])
    end
    return results
end

function analyze_results(results, num_simulations::Int64, T::Int64)
    amount_histories = [];
    liquidity_histories = [];
    num_time_in_low_liquidity = 0
    num_time_in_low_return = 0
    total_num_periods = T * num_simulations
    for result in results
        amount_history = [];
        liquidity_history = [];
        amount = 100.0;
        for state in result
            append!(amount_history, copy(amount));
            if is_high_return_state(state)
                amount = high_return_result(amount);
            else is_low_return_state(state)
                num_time_in_low_return += 1;
                amount = low_return_result(amount);
            end
            if is_high_liquidity_state(state)
                append!(liquidity_history, 1);
            else is_low_liquidity_state(state)
                num_time_in_low_liquidity += 1;
                append!(liquidity_history, 0);
            end
        end
        append!(liquidity_histories, [liquidity_history]);
        append!(amount_histories, [amount_history]);
    end
    fraction_of_time_in_low_liquidity = num_time_in_low_liquidity / total_num_periods
    fraction_of_time_in_low_return = num_time_in_low_return / total_num_periods
    return (amount_histories, liquidity_histories, fraction_of_time_in_low_return, fraction_of_time_in_low_liquidity)
end

P_a = [.9603 .0297 .0097 .0003; .891 .099 .009 .001; .125 .125 .375 .375; .025 .225 .075 .675];
pi = [1 0 0 0];
st_b = MarkovChain(P_a, pi);
res_b = run_analysis_for_four(st_b, 4, 100);
long_run_return_b = find_stationary_distribution(st_b);
long_run_low_return_b = long_run_return_b[2] + long_run_return_b[4];
long_run_low_liquidity_b = long_run_return_b[3] + long_run_return_b[4];

P_c = [.9603 .0297 .0097 .0003; .9603 .0297 .0097 .0003; .005 .005 .495 .495; .005 .005 .495 .495];
st_c = MarkovChain(P_c, pi);
res_c = run_analysis_for_four(st_c, 4, 100);
long_run_return_c = find_stationary_distribution(st_c);
long_run_low_return_c = long_run_return_c[2] + long_run_return_c[4];
long_run_low_liquidity_c = long_run_return_c[3] + long_run_return_c[4];
