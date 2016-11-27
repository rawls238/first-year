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
country_offset = 18;
state_keys = {'Kt', 'Kt1', 'Kt2', 'Kt3', 'Zt', 'Zt1', 'Zt2', 'Zt3', 'lambdat'};
state_vals = zeros(1, length(state_keys));
for i=1:length(state_keys)
    state_vals(i) = i;
end
state_offset = length(state_vals);
control_keys = {'Ct', 'Ct1', 'Ct2', 'Ct3', 'Nt', 'Nt1', 'Nt2', 'Nt3', 'lambdat1', 'lambdat2', 'lambdat3'};
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

m = 1 + delta * beta + delta * beta^2 + delta * beta^3 + (delta-1) * beta^4;
m_inv = 1 / m;

for country_num=1:num_countries
    Ct = var_offset(country_num, vars('Ct'));
    Kt = var_offset(country_num, vars('Kt'));
    lambdat = var_offset(country_num, vars('lambdat'));
    Zt = var_offset(country_num, vars('Zt'));
    Nt = var_offset(country_num, vars('Nt'));

    Ct1 = var_offset(country_num, vars('Ct1'));
    Kt1 = var_offset(country_num, vars('Kt1'));
    lambdat1 = var_offset(country_num, vars('lambdat1'));
    Zt1 = var_offset(country_num, vars('Zt1'));
    Nt1 = var_offset(country_num, vars('Nt1'));

    Ct2 = var_offset(country_num, vars('Ct2'));
    Kt2 = var_offset(country_num, vars('Kt2'));
    lambdat2 = var_offset(country_num, vars('lambdat2'));
    Zt2 = var_offset(country_num, vars('Zt2'));
    Nt2 = var_offset(country_num, vars('Nt2'));

    Ct3 = var_offset(country_num, vars('Ct3'));
    Kt3 = var_offset(country_num, vars('Kt3'));
    lambdat3 = var_offset(country_num, vars('lambdat3'));
    Zt3 = var_offset(country_num, vars('Zt3'));
    Nt3 = var_offset(country_num, vars('Nt3'));

    offset = (country_num - 1) * country_offset;
    %1. Inventory FOC
    eqn = 1 + offset;
    Fy(eqn,Ct3) = xi_c_c;
    Fy(eqn,Nt3) = xi_c_n;
    Fyp(eqn,Ct3) = RHS_flag * xi_c_c;
    Fyp(eqn,Nt3) = RHS_flag * ((1 - beta) * tau_z_n + xi_c_n);
    Fyp(eqn,lambdat3) = RHS_flag * (1 - beta) * tau_z_lambda;
    Fxp(eqn,Kt3) = RHS_flag * (1 - beta) * tau_z_k;
    Fxp(eqn,Zt3) = RHS_flag * (1 - beta) * tau_z_z;

    %2. Euler
    eqn = 2 + offset;
    Fyp(eqn, Ct3) = RHS_flag * (xi_c_c - (m_inv * (delta - 1) * beta^4 * xi_c_c));
    Fyp(eqn, Nt3) = RHS_flag * (xi_c_n - (m_inv * (delta - 1) * beta^4 * xi_c_n));
    Fyp(eqn, lambdat3) = RHS_flag * tau_k_lambda;
    Fxp(eqn, Kt3) = RHS_flag * tau_k_k;
    Fxp(eqn, Zt3) = RHS_flag * tau_k_z;
    Fy(eqn, Ct) = m_inv * xi_c_c;
    Fy(eqn, Nt) = m_inv * xi_c_n;
    Fy(eqn, Ct1) = m_inv * xi_c_c * delta * beta;
    Fy(eqn, Nt1) = m_inv * xi_c_n * delta * beta;
    Fy(eqn, Ct2) = m_inv * xi_c_c * delta * beta^2;
    Fy(eqn, Nt2) = m_inv * xi_c_n * delta * beta^2;
    Fy(eqn, Ct3) = m_inv * xi_c_c * delta * beta^3;
    Fy(eqn, Nt3) = m_inv * xi_c_n * delta * beta^3;

    %4. Labor
    eqn = 3 + offset;
    Fy(eqn, Ct) = xi_n_c - xi_c_c;
    Fy(eqn, Nt) = xi_n_n - xi_c_n - tau_n_n;
    Fx(eqn, lambdat) = RHS_flag * tau_n_lambda;
    Fx(eqn, Kt) = RHS_flag * tau_n_k;
    Fp(eqn, Zt) = RHS_flag * tau_n_z;
    
    eqn = 4 + offset;
    Fx(eqn,Kt1) = 1;
    Fxp(eqn,Kt) = RHS_flag;

    eqn = 5 + offset;
    Fx(eqn,Kt2) = 1;
    Fxp(eqn,Kt1) = RHS_flag;

    eqn = 6 + offset;
    Fx(eqn,Kt3) = 1;
    Fxp(eqn,Kt2) = RHS_flag;

    eqn = 7 + offset;
    Fy(eqn, lambdat1) = 1;
    Fxp(eqn, lambdat) = RHS_flag;

    eqn = 8 + offset;
    Fy(eqn, lambdat2) = 1;
    Fyp(eqn, lambdat1) = RHS_flag;

    eqn = 9 + offset;
    Fy(eqn, lambdat3) = 1;
    Fyp(eqn, lambdat2) = RHS_flag;

    eqn = 10 + offset;
    Fx(eqn, Zt1) = 1;
    Fxp(eqn,Zt) = RHS_flag;

    eqn = 11 + offset;
    Fx(eqn, Zt2) = 1;
    Fxp(eqn,Zt1) = RHS_flag;

    eqn = 12 + offset;
    Fx(eqn, Zt3) = 1;
    Fxp(eqn,Zt2) = RHS_flag;

    eqn = 13 + offset;
    Fy(eqn, Ct1) = 1;
    Fyp(eqn, Ct) = RHS_flag;

    eqn = 14 + offset;
    Fy(eqn, Ct2) = 1;
    Fyp(eqn, Ct1) = RHS_flag;

    eqn = 15 + offset;
    Fy(eqn, Ct3) = 1;
    Fyp(eqn, Ct2) = RHS_flag;
    
    eqn = 16 + offset;
    Fy(eqn, Nt1) = 1;
    Fyp(eqn, Nt) = RHS_flag;

    eqn = 17 + offset;
    Fy(eqn, Nt2) = 1;
    Fyp(eqn, Nt1) = RHS_flag;

    eqn = 18 + offset;
    Fy(eqn, Nt3) = 1;
    Fyp(eqn, Nt2) = RHS_flag;
end

Ct_h = var_offset(1, vars('Ct'));
Ct_f = var_offset(2, vars('Ct'));
Kt_h = var_offset(1, vars('Kt'));
Kt_f = var_offset(2, vars('Kt'));
lambdat_h = var_offset(1, vars('lambdat'));
lambdat_f = var_offset(2, vars('lambdat'));
Zt_h = var_offset(1, vars('Zt'));
Zt_f = var_offset(2, vars('Zt'));
Nt_h = var_offset(1, vars('Nt'));
Nt_f = var_offset(2, vars('Nt'));
Kt1_h = var_offset(1, vars('Kt1'));
Kt1_f = var_offset(2, vars('Kt1'));

Kt2_h = var_offset(1, vars('Kt2'));
Kt2_f = var_offset(2, vars('Kt2'));

Kt3_h = var_offset(1, vars('Kt3'));
Kt3_f = var_offset(2, vars('Kt3'));

eqn = 37;
Fxp(eqn,lambdat_h)  = RHS_flag;
Fx(eqn,lambdat_h)  = A(1, 1);
Fx(eqn,lambdat_f) = A(1, 2);

eqn = 38;
Fxp(eqn,lambdat_f) = RHS_flag;
Fx(eqn,lambdat_h) = A(2, 1);
Fx(eqn,lambdat_f) = A(2, 2);

%1. Marginal Utility eq
eqn = 39;
Fy(eqn, Ct_h) = xi_c_c;
Fy(eqn, Nt_h) = xi_c_n;
Fy(eqn, Ct_f) = RHS_flag * xi_c_c;
Fy(eqn, Nt_f) = RHS_flag * xi_c_n;

%Resource constraint
eqn = 40;
Fx(eqn,lambdat_h) = RHS_flag * yss * zita_lambda;
Fx(eqn,lambdat_f) = RHS_flag * yss * zita_lambda;
Fy(eqn,Nt_h) = RHS_flag * yss * zita_n;
Fy(eqn,Nt_f) = RHS_flag * yss * zita_n;
Fx(eqn,Zt_h) = RHS_flag * yss * zita_z;
Fx(eqn,Zt_f) = RHS_flag * yss * zita_z;
Fy(eqn,Ct_h) = css;
Fy(eqn,Ct_f) = css;
Fx(eqn,Kt_h) = RHS_flag * (0.25 * (1 - delta) * kss + yss * zita_k);
Fx(eqn,Kt_f) = RHS_flag * (0.25 * (1 - delta) * kss + yss * zita_k);
Fx(eqn,Kt1_h) = RHS_flag * 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt1_f) = RHS_flag * 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt2_h) = RHS_flag * 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt2_f) = RHS_flag * 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt3_h) = RHS_flag * 0.25 * ((1 - delta) - 1) * kss;
Fx(eqn,Kt3_f) = RHS_flag * 0.25 * ((1 - delta) - 1) * kss;
Fxp(eqn,Kt3_h) = RHS_flag * 0.25 * -1 * kss;
Fxp(eqn,Kt3_f) = RHS_flag * 0.25 * -1 * kss;


At = [-Fxp -Fyp]; Bt = [Fx Fy];
[H,G]=solab(At,Bt,size(Fx,2));