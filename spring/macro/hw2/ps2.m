load data_hwk2.mat;

estimates = zeros(length(data_hwk2), 3);
r = 0.08;
for i=1:length(data_hwk2)
   dat = data_hwk2{i,2};
   estimate = estimate_dat(log(dat(2,:)));
   estimates(i,:) = [estimate(2); estimate(3); is_counter_cyclical(estimate, r)];
end

x = linspace(0, 2);
c = figure();
scatter(estimates(:,1), estimates(:,2));
hold on;
plot(x, -1*ones(size(x)));
hold on;
plot(x, 1-x);
hold on;
plot(x, 1+x);
hold on;
plot(x, (1+r)*(1-x));

str_country = strrep({data_hwk2{:,1}}, ' ', '_');
printmat(estimates, 'AR(2) estimates', strjoin(str_country), 'rho_1 rho_2 countercyclical');

function [res] = is_counter_cyclical(vals, r)
    res = vals(2) > 1 && vals(3) < 0 && vals(3) > (1+r)*(1-vals(2));
end

function [estimates] = estimate_dat(data)  
    s = size(data);
    t = s(2);
    
    % Detrend, 1.b
    detrended = detrend_dat(t, data);
    dat = exp(detrended); %1.c
    
    %1.d
    VARn = 1;
    VARp = 2;
    Z = [ones(t, 1) lagmatrix(dat',1) lagmatrix(dat',2)];
    Z = Z(VARp+1:length(Z),:)';
    
    dat = dat(:,VARp+1:length(dat));

    % GLS estimator
    estimates  = kron(inv(Z*Z')*Z,eye(VARn))*dat(:);
end

function [detrended] = detrend_dat(t, dat)
    t1 = 1:t;
    t2 = t1.^2;
    y = dat';
    X = [ones(1, t); t1; t2];
    beta = inv(X * X') * (X * y);
    proj = beta(1) + beta(2) * t1 + beta(3) * t2;
    detrended = dat - proj;
end