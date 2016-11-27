clear all;
 close all;
 
  % Parameter values
 mu    = 0.34;
 beta = 0.99;
 alpha = 1.0;
 gamma = -1.0;
 theta = 0.36;
 v = 3.0;
 sigma = 0.01;
 delta = 0.025;
 J = 4;
 A = [.906 .088; .088 .906];
 
 % Steady state
 yss = 1;    % Output is normalized to 1
 zss = yss * ((1 - beta) / (beta * sigma))^(1/(-v - 1));
 kss_denom = 1 + (delta - 1)*beta;
 kss = (theta * beta / kss_denom) * yss^(v+1) * ((yss)^(-v) - sigma * (zss)^(-v))^(-1);
 nss_first_term = ((1 - theta) * (yss)^(v+1)*((yss)^(-v) - sigma*(zss)^(-v)));
 nss_second_term = (1/mu - 1) * (yss - delta * kss);
 nss = (nss_first_term * nss_second_term + 1)^-1;
 css = yss - delta * kss;
 
 tmp = sigma * yss^v * zss^-v;
 xi_c_c = gamma*mu - 1;
 xi_c_l = gamma - gamma * mu;
 xi_l_c = gamma * mu;
 xi_l_l = -1 + gamma - gamma * mu;
 xi_n_n = xi_l_l * (nss / (nss - 1));
 xi_c_n = gamma * (1 - mu) * nss / (nss - 1);
 xi_n_c = gamma * mu;
 zita_k = theta * (1 - tmp);
 zita_n = (1 - theta) * (1 - tmp);
 zita_z = tmp;
 zita_lambda = 1 - tmp;
 tau_k_k = (v + 1) * theta * (1 - tmp) - v * theta - 1;
 tau_k_n = (v + 1) * (1 - theta) * (1 - tmp) + v * (theta - 1);
 tau_k_z = (v + 1) * tmp;
 tau_n_z = tau_k_z;
 tau_n_k = (v + 1) * theta * (1 - tmp) - v * theta;
 tau_n_n = (v + 1) * (1 - theta) * (1 - tmp) + v * (theta - 1) - 1;
 tau_z_k = (v + 1) * theta * (1 - tmp);
 tau_z_n = (v + 1) * (1 - theta) * (1 - tmp);
 tau_z_z = (v + 1) * tmp - v - 1;
 tau_z_lambda = (v + 1) * (1 - tmp);
 tau_n_lambda = (v+1)*(1 - (sigma)*(yss)^v*(zss)^(-v)) - v;
 tau_k_lambda = (v+1)*(1 - (sigma)*(yss)^v*(zss)^(-v)) - v; 

num_countries = 2;
country_offset = 4;
state_keys = {'Kt', 'Zt', 'lambdat', 'Kt1', 'Kt2', 'Kt3', 'Zt1', 'Zt2', 'Zt3', 'lambdat1'};
state_vals = zeros(1, length(state_keys));
for i=1:length(state_keys)
    state_vals(i) = i;
end
state_offset = length(state_vals);
control_keys = {'Ct', 'Ct1', 'Ct2', 'Ct3', 'Yt', 'Yt1', 'Yt2', 'Yt3', 'Nt1', 'Nt2', 'Nt3', 'lambdat2', 'lambdat3'};
control_vals = zeros(1, length(control_keys));
for i=1:length(control_keys)
    control_vals(i) = i;
end
control_offset = length(control_vals);
val_offset = [state_vals length(state_vals) + control_vals];
vars = containers.Map([state_keys control_keys], val_offset);
total_vals = [state_vals control_vals];
total_eqs = num_countries * country_offset + 4;
total_controls = num_countries * control_offset;
total_states = num_countries * state_offset;

var_offset = zeros(num_countries, length(total_vals));
for i=1:num_countries
    for j=1:length(state_vals)
        var_offset(i, j) = state_offset * (i - 1) + state_vals(j);
    end
    start = length(state_vals) + 1;
    done = length(state_vals) + length(control_vals);
    for j=start:done
        var_offset(i, j) = control_offset * (i - 1) + control_vals(j-length(state_vals));
    end
end


Fy=zeros(total_eqs,total_controls);  Fx=zeros(total_eqs,total_states);
Fyp=zeros(total_eqs,total_controls); Fxp=zeros(total_eqs,total_states);

RHS_flag = -1;

for country_num=1:num_countries
    
Ct = var_offset(country_num, vars('Ct'));
Kt = var_offset(country_num, vars('Kt'));
lambdat = var_offset(country_num, vars('lambdat'));
Zt = var_offset(country_num, vars('Zt'));
Nt = var_offset(country_num, vars('Nt'));
Yt = var_offset(country_num, vars('Yt'));

Ct1 = var_offset(country_num, vars('Ct1'));
Kt1 = var_offset(country_num, vars('Kt1'));
lambdat1 = var_offset(country_num, vars('lambdat1'));
Zt1 = var_offset(country_num, vars('Zt1'));
Nt1 = var_offset(country_num, vars('Nt1'));
Yt1 = var_offset(country_num, vars('Yt1'));

Ct2 = var_offset(country_num, vars('Ct2'));
Kt2 = var_offset(country_num, vars('Kt2'));
lambdat2 = var_offset(country_num, vars('lambdat2'));
Zt2 = var_offset(country_num, vars('Zt2'));
Nt2 = var_offset(country_num, vars('Nt2'));
Yt2 = var_offset(country_num, vars('Yt2'));

Ct3 = var_offset(country_num, vars('Ct3'));
Kt3 = var_offset(country_num, vars('Kt3'));
lambdat3 = var_offset(country_num, vars('lambdat3'));
Zt3 = var_offset(country_num, vars('Zt3'));
Nt3 = var_offset(country_num, vars('Nt3'));
Yt3 = var_offset(country_num, vars('Yt3'));


%2. Inventory FOC
eqn = 1 + (country_num - 1) * country_offset;
Fy(eqn,C) = xi_c_c;
Fy(eqn,N) = xi_c_n;
Fyp(eqn,C) = xi_c_c;
Fyp(eqn,N) = xi_c_n;
Fxp(eqn,lambda) = RHS_flag * (1 - beta) * tau_z_lambda;
Fxp(eqn,K) = RHS_flag * (1 - beta) * tau_z_k;
Fxp(eqn,Z) = RHS_flag * (1 - beta) * tau_z_z;
Fyp(eqn,N) = RHS_flag * (1 - beta) * tau_z_n;

%3. Euler
mid = 1 / (1 + (delta - 1) * beta);
eqn = 2 + (country_num - 1) * country_offset;
Fyp(eqn, C) = RHS_flag * (xi_c_c - mid * (delta - 1) * beta * xi_c_c);
Fyp(eqn, N) = RHS_flag * (xi_c_n - mid * (delta - 1) * beta * xi_c_n);
Fxp(eqn, lambda) = RHS_flag * tau_k_lambda;
Fxp(eqn, K) = RHS_flag * tau_k_k;
Fxp(eqn, Z) = RHS_flag * tau_k_z;
Fy(eqn, C) = mid * xi_c_c;
Fy(eqn, N) = mid * xi_c_n;

%Production Function
eqn = 3 + (country_num - 1) * country_offset;
Fy(eqn,Y) = RHS_flag;
Fx(eqn,lambda) = zita_lambda;
Fy(eqn,N) = zita_n;
Fx(eqn,K) = zita_k;
Fx(eqn,Z) = zita_z;

% Labor
eqn = 4 + (country_num - 1) * country_offset;
Fy(eqn, C) = xi_n_c - xi_c_c;
Fy(eqn, N) = xi_n_n - xi_c_n;
Fx(eqn, lambda) = RHS_flag * tau_n_lambda;
Fx(eqn, K) = RHS_flag * tau_n_k;
Fy(eqn, N) = tau_n_n;
Fx(eqn, Z) = tau_n_z;
end

Ct_h = var_offset(1, vars('Ct'));
Ct_f = var_offset(2, vars('Ct'));
Kt_h = var_offset(1, vars('Kt'));
Kt_f = var_offset(2, vars('Kt'));
lambdat_h = var_offset(1, vars('lambdat'));
lambdat_f = var_offset(2, vars('lambdat'));
Zt_h = var_offset(1, vars('Zt'));
Zt_f = var_offset(2, vars('Zt'));
Yt_h = var_offset(1, vars('Yt'));
Yt_f = var_offset(2, vars('Yt'));
Nt_h = var_offset(1, vars('Nt'));
Nt_f = var_offset(2, vars('Nt'));

Ct1_h = var_offset(1, vars('Ct1'));
Ct1_f = var_offset(2, vars('Ct1'));
Kt1_h = var_offset(1, vars('Kt1'));
Kt1_f = var_offset(2, vars('Kt1'));
lambdat1_h = var_offset(1, vars('lambdat1'));
lambdat1_f = var_offset(2, vars('lambdat1'));
Zt1_h = var_offset(1, vars('Zt1'));
Zt1_f = var_offset(2, vars('Zt1'));
Yt1_h = var_offset(1, vars('Yt1'));
Yt1_f = var_offset(2, vars('Yt1'));
Nt1_h = var_offset(1, vars('Nt1'));
Nt1_f = var_offset(2, vars('Nt1'));


Kt2_h = var_offset(1, vars('Kt2'));
Kt2_f = var_offset(2, vars('Kt2'));

Kt3_h = var_offset(1, vars('Kt3'));
Kt3_f = var_offset(2, vars('Kt3'));


eqn = 9;
Fxp(eqn,lambda_h)  = RHS_flag;
Fx(eqn,lambda_h)  = A(1, 1);
Fx(eqn,lambda_f) = A(1, 2);

eqn = 10;
Fxp(eqn,lambda_f) = RHS_flag;
Fx(eqn,lambda_h) = A(2, 1);
Fx(eqn,lambda_f) = A(2, 2);

%1. Marginal Utility eq
eqn = 11;
Fy(eqn, C_h) = xi_c_c;
Fy(eqn, N_h) = xi_c_n;
Fy(eqn, C_f) = RHS_flag * xi_c_c;
Fy(eqn, N_f) = RHS_flag * xi_c_n;

%Resource constraint
eqn = 12;
Fy(eqn,Y_h) = yss;
Fy(eqn,Y_f) = yss;
Fy(eqn,C_h) = css;
Fy(eqn,C_f) = css;
Fx(eqn,Kt_h) = 0.25 * (1 - delta) * kss;
Fx(eqn,Kt_f) = 0.25 * (1 - delta) * kss;
Fx(eqn,Kt1_h) = 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt1_f) = 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt2_h) = 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt2_f) = 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt3_h) = 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt3_f) = 0.25 * ((1 - delta) - 1) * kss;
Fxp(eqn,Kt3_h) = RHS_flag * 0.25 * -1 * kss;
Fxp(eqn,Kt3_f) = RHS_flag * 0.25 * -1 * kss;


At = [-Fxp -Fyp]; Bt = [Fx Fy];
[H,G]=solab(At,Bt,size(Fx,2));