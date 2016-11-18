data = csvread('/Users/garidor/Desktop/first-year/macro/data.csv', 1);

realgdp = data(:, 3);
realwage = data(:, 2);
hours_worked = data(:, 5);

% 2.2
empirical_labor = figure();
plot(data(:, 1), hpfilter(log(realgdp), 1600))
hold on;
plot(data(:, 1), hpfilter(log(hours_worked), 1600))
xlabel('Quarters')
legend('y', 'n')
title('Empirical Hours vs GDP')

model_labor = figure();
plot(1:200, hpfilter(log(sim_y_series), 1600))
hold on;
plot(1:200, hpfilter(log(sim_n_series), 1600))
title('Simulated Hours vs GDP')
legend('sim_y', 'sim_n')
xlabel('Quarters')

corrcoef(realgdp, hours_worked)
corrcoef(sim_y_series, sim_n_series)

%2.3
corrcoef(realwage, realgdp)
corrcoef(sim_y_series, sim_wage)

%2.4
in_model = compute_growth(sim_y_series);
in_model_corr = autocorr(in_model, 1)
in_data = compute_growth(realgdp);
in_data_cor = autocorr(in_data, 1)

%2.5
solow_residual_empirical = log(realgdp) + (1 - alfa) * log(data(:, 4)) - alfa * log(data(:, 5));
coeff_of_var_empirical = std(solow_residual_empirical) / mean(solow_residual_empirical)
coeff_of_var_model = std(solow_residual) / mean(solow_residual)
coeff_of_var_empirical_y = std(realgdp) / mean(realgdp)
coeff_of_var_model_y = std(sim_y_series) / mean(sim_y_series)
solow_empirical = figure();
plot(1:length(solow_residual_empirical), solow_residual_empirical)
title('Solow Residual Empirical')
xlabel('Quarters')
solow_model = figure();
plot(1:length(solow_residual), solow_residual)
title('Solow Residual Model')
xlabel('Quarters')

function growth = compute_growth(d) 
growth = zeros(1, length(d)-1);
for i = 2:(length(d))
    growth(i-1) = d(i) / d(i-1);
end
end