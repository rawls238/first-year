data = csvread('/Users/garidor/Desktop/first-year/macro/data.csv', 1);

realgdp = data(:, 3)
corrcoef(data(:,2), realgdp)
corrcoef(sim_y_series, sim_wage_2)

in_model = compute_growth(sim_y_series);
in_model_corr = autocorr(in_model, 1)
in_data = compute_growth(realgdp);
in_data_cor = autocorr(in_data, 1)

function growth = compute_growth(d) 
growth = zeros(1, length(d)-1);
for i = 2:(length(d))
    growth(i-1) = d(i) / d(i-1);
end
end