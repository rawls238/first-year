% Replication of  Backus, Kehoe and Kydland (1992): International Real
% Business Cycles, Journal of Political Economy

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
xss = delta * kss;

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
country_offset = 22;
state_keys = {'Kt', 'Kt1', 'Kt2', 'Kt3', 'Zt', 'Zt1', 'Zt2', 'Zt3', 'lambdat'};
state_vals = zeros(1, length(state_keys));
for i=1:length(state_keys)
    state_vals(i) = i;
end
state_offset = length(state_vals);
control_keys = {'Ct', 'Ct1', 'Ct2', 'Ct3', 'Yt', 'Yt1', 'Yt2', 'Yt3', 'Nt', 'Nt1', 'Nt2', 'Nt3', 'lambdat1', 'lambdat2', 'lambdat3'};
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
    Fyp(eqn, Nt3) = RHS_flag * (xi_c_n + tau_k_n - (m_inv * (delta - 1) * beta^4 * xi_c_n));
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
    Fx(eqn, Zt) = RHS_flag * tau_n_z;
    
    eqn = 4 + offset;
    Fy(eqn,Yt) = RHS_flag;
    Fx(eqn,lambdat) = zita_lambda;
    Fy(eqn,Nt) = zita_n;
    Fx(eqn,Kt) = zita_k;
    Fx(eqn,Zt) = zita_z;
    
    eqn = 5 + offset;
    Fx(eqn,Kt1) = 1;
    Fxp(eqn,Kt) = RHS_flag;

    eqn = 6 + offset;
    Fx(eqn,Kt2) = 1;
    Fxp(eqn,Kt1) = RHS_flag;

    eqn = 7 + offset;
    Fx(eqn,Kt3) = 1;
    Fxp(eqn,Kt2) = RHS_flag;

    eqn = 8 + offset;
    Fy(eqn, lambdat1) = 1;
    Fxp(eqn, lambdat) = RHS_flag;

    eqn = 9 + offset;
    Fy(eqn, lambdat2) = 1;
    Fyp(eqn, lambdat1) = RHS_flag;

    eqn = 10 + offset;
    Fy(eqn, lambdat3) = 1;
    Fyp(eqn, lambdat2) = RHS_flag;

    eqn = 11 + offset;
    Fx(eqn, Zt1) = 1;
    Fxp(eqn,Zt) = RHS_flag;

    eqn = 12 + offset;
    Fx(eqn, Zt2) = 1;
    Fxp(eqn,Zt1) = RHS_flag;

    eqn = 13 + offset;
    Fx(eqn, Zt3) = 1;
    Fxp(eqn,Zt2) = RHS_flag;

    eqn = 14 + offset;
    Fy(eqn, Ct1) = 1;
    Fyp(eqn, Ct) = RHS_flag;

    eqn = 15 + offset;
    Fy(eqn, Ct2) = 1;
    Fyp(eqn, Ct1) = RHS_flag;

    eqn = 16 + offset;
    Fy(eqn, Ct3) = 1;
    Fyp(eqn, Ct2) = RHS_flag;
    
    eqn = 17 + offset;
    Fy(eqn, Nt1) = 1;
    Fyp(eqn, Nt) = RHS_flag;

    eqn = 18 + offset;
    Fy(eqn, Nt2) = 1;
    Fyp(eqn, Nt1) = RHS_flag;

    eqn = 19 + offset;
    Fy(eqn, Nt3) = 1;
    Fyp(eqn, Nt2) = RHS_flag;
    
    eqn = 20 + offset;
    Fy(eqn, Yt1) = 1;
    Fyp(eqn, Yt) = RHS_flag;

    eqn = 21 + offset;
    Fy(eqn, Yt2) = 1;
    Fyp(eqn, Yt1) = RHS_flag;
 
    eqn = 22 + offset;
    Fy(eqn, Yt3) = 1;
    Fyp(eqn, Yt2) = RHS_flag;
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
Yt_h = var_offset(1, vars('Yt'));
Yt_f = var_offset(2, vars('Yt'));
Kt1_h = var_offset(1, vars('Kt1'));
Kt1_f = var_offset(2, vars('Kt1'));

Kt2_h = var_offset(1, vars('Kt2'));
Kt2_f = var_offset(2, vars('Kt2'));

Kt3_h = var_offset(1, vars('Kt3'));
Kt3_f = var_offset(2, vars('Kt3'));

%Tech progress
eqn = 45;
Fxp(eqn,lambdat_h)  = RHS_flag;
Fx(eqn,lambdat_h)  = A(1, 1);
Fx(eqn,lambdat_f) = A(1, 2);

eqn = 46;
Fxp(eqn,lambdat_f) = RHS_flag;
Fx(eqn,lambdat_h) = A(2, 1);
Fx(eqn,lambdat_f) = A(2, 2);

%1. Marginal Utility eq
eqn = 47;
Fy(eqn, Ct_h) = xi_c_c;
Fy(eqn, Nt_h) = xi_c_n;
Fy(eqn, Ct_f) = RHS_flag * xi_c_c;
Fy(eqn, Nt_f) = RHS_flag * xi_c_n;

%Resource constraint
eqn = 48;
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

%% Impulse responses
 shock(:,1) = [0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0];
 for i=1:20
    shock(:,i+1) = G*shock(:,i);
 end

z = (H*shock)';
shock = shock';

lambda_h = shock(:,lambdat_h);
lambda_f = shock(:,lambdat_f);
k_h = shock(:,Kt_h);
k1_h = shock(:,Kt1_h);
k2_h = shock(:,Kt2_h);
k3_h = shock(:,Kt3_h);
k_f = shock(:,Kt_f);
k1_f = shock(:,Kt1_f);
k2_f = shock(:,Kt2_f);
k3_f = shock(:,Kt3_f);
c_h = z(:,Ct_h);
c_f= z(:,Ct_f);
zi_h = shock(:,Zt_h);
zi_f = shock(:,Zt_f);

y_h = z(:,Yt_h);
y_f = z(:,Yt_f);
i_h = zeros(20, 1);
i_f = zeros(20, 1);
x_h = zeros(20, 1);
x_f = zeros(20, 1);
nx_h = zeros(20, 1);
nx_f = zeros(20, 1);
for i=1:20
    x_h(i) = (1/(delta*4))* ((delta-1)*k_h(i) + delta*k1_h(i) + delta*k2_h(i) + delta*k3_h(i) + k3_h(i+1));
    x_f(i) = (1/(delta*4))* ((delta-1)*k_f(i) + delta*k1_f(i) + delta*k2_f(i) + delta*k3_f(i) + k3_f(i+1));
    i_h(i) = x_h(i);
    i_f(i) = x_f(i);
    nx_h(i) = y_h(i) - c_h(i) - i_h(i);
    nx_f(i) = y_f(i) - c_f(i) - i_f(i);
end

 % Plot Impulse Responses
homeResponse = figure();
plot(0:21, [0 c_h'],'-d','MarkerSize',3,'Color',[0,0,0])
hold on 
plot(0:21, [0 y_h'],'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
hold on 
plot(0:20, [0 i_h'],'-.','MarkerSize',3,'Color',[0,.5,0])
hold on
plot(0:21, [0 lambda_h'], '-x','MarkerSize',3,'Color',[0,0,.5])
hold on
plot(0:20, [0 nx_h'], '-y', 'MarkerSize', 3,'Color',[0.5,0.5,.5])
hold off
ylabel('Percent Deviations')
xlabel('Quarters')
title('Home response')
legend('c', 'y','i','lambda', 'nx');

foreignResponse = figure();
plot(0:21, [0 c_f'],'-d','MarkerSize',3,'Color',[0,0,0])
hold on 
plot(0:21, [0 y_f'],'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
hold on 
plot(0:20, [0 i_f'],'-.','MarkerSize',3,'Color',[0,.5,0])
hold on
plot(0:21, [0 lambda_f'], '-x','MarkerSize',3,'Color',[0,0,.5])
hold on
plot(0:20, [0 nx_f'], '-y', 'MarkerSize', 3,'Color',[0.5,0.5,.5])
hold off
ylabel('Percent Deviations')
xlabel('Quarters')
title('Foreign response')
legend('c', 'y','i','lambda', 'nx');


period = 100;
sim_num = 50;
total_periods = period * sim_num;
mu = [0,0];
sigma = [.00852^2,.00852^2 * .258; .00852^2 * .258, .00852^2];
e = mvnrnd(mu,sigma,total_periods+1);
s = [0;0;0;0;0;0;0;0;e(1,1);0;0;0;0;0;0;0;0;e(1,2)];
for i=1:total_periods
    s(:,i+1) = G*s(:,i) + [0;0;0;0;0;0;0;0;e(i,1);0;0;0;0;0;0;0;0;e(i,2)];
end
z = (H*s)';
s = s';


k_h = s(:,Kt_h);
c_h = z(:,Ct_h);
n_h = z(:,Nt_h);
z_h = s(:,Zt_h);
k_f = s(:,Kt_f);
c_f = z(:,Ct_f);
n_f = z(:,Nt_f);
z_f = s(:,Zt_f);
k1_h = s(:,Kt1_h);
k2_h = s(:,Kt2_h);
k3_h = s(:,Kt3_h);
k1_f = s(:,Kt1_f);
k2_f = s(:,Kt2_f);
k3_f = s(:,Kt3_f);
y_h = z(:,Yt_h);
y_f = z(:,Yt_f);
i_h = zeros(total_periods, 1);
i_f = zeros(total_periods, 1);
x_h = zeros(total_periods, 1);
x_f = zeros(total_periods, 1);
nx_h = zeros(total_periods, 1);
nx_f = zeros(total_periods, 1);
for i=1:total_periods
    x_h(i) = (1/(delta*4))* ((delta-1)*k_h(i) + delta*k1_h(i) + delta*k2_h(i) + delta*k3_h(i) + k3_h(i+1));
    x_f(i) = (1/(delta*4))* ((delta-1)*k_f(i) + delta*k1_f(i) + delta*k2_f(i) + delta*k3_f(i) + k3_f(i+1));
    i_h(i) = x_h(i);
    i_f(i) = x_f(i);
    nx_h(i) = y_h(i) - c_h(i) - i_h(i);
    nx_f(i) = y_f(i) - c_f(i) - i_f(i);
end

nxss = yss - css - xss;
lag = 5;
total_cross_corr = lag * 2 + 1;
standard_deviations = zeros(total_periods/period, 7);
std_dev_relative = zeros(total_periods/period, 6);
cross_corr_y = zeros(total_periods/period, total_cross_corr);
cross_corr_n = zeros(total_periods/period, total_cross_corr);
cross_corr_k = zeros(total_periods/period, total_cross_corr);
cross_corr_c = zeros(total_periods/period, total_cross_corr);
cross_corr_x = zeros(total_periods/period, total_cross_corr);
cross_corr_nx = zeros(total_periods/period, total_cross_corr);
cross_corr_z = zeros(total_periods/period, total_cross_corr);
cross_corr_y_across_countries = zeros(total_periods/period, total_cross_corr);
cross_corr_c_across_countries = zeros(total_periods/period, total_cross_corr);
cross_corr_saving = zeros(total_periods/period, total_cross_corr);
for j = 1:total_periods/period
    yh_1 = log(yss * exp(y_h(1+(j-1)*period:j*period, 1)));
    lyh_1 = yh_1 - hpfilter(yh_1, 1600);
    yf_1 = log(yss * exp(y_f(1+(j-1)*period:j*period, 1)));
    lyf_1 = yf_1 - hpfilter(yf_1, 1600);
    ch_1 = log(css * exp(c_h(1+(j-1)*period:j*period, 1)));
    lch_1 = ch_1 - hpfilter(ch_1, 1600);
    cf_1 = log(css * exp(c_f(1+(j-1)*period:j*period, 1)));
    lcf_1 = cf_1 - hpfilter(cf_1, 1600);
    nh_1 = log(nss * exp(n_h(1+(j-1)*period:j*period, 1)));
    lnh_1 = nh_1 - hpfilter(nh_1, 1600);
    kh_1 = log(kss * exp(k_h(1+(j-1)*period:j*period, 1)));
    lkh_1 = kh_1 - hpfilter(kh_1, 1600);
    xh_1 = log(xss * exp(x_h(1+(j-1)*period:j*period, 1)));
    lxh_1 = xh_1 - hpfilter(xh_1, 1600);
    zh_1 = log(zss * exp(z_h(1+(j-1)*period:j*period, 1)));
    lzh_1 = zh_1 - hpfilter(zh_1, 1600);
    nx = exp(nx_h(1+(j-1)*period:j*period, 1));
    y1_h = exp(y_h(1+(j-1)*period:j*period, 1));
    nxh_1 = log(nxss/yss .* nx .* y1_h);
    lnxh_1 = nxh_1 - hpfilter(nxh_1, 1600);
    c1_h = exp(c_h(1+(j-1)*period:j*period, 1));
    saving_h = log(yss*(css)/yss .* (c1_h) .* y1_h);
    lsavingh_1 = saving_h - hpfilter(saving_h, 1600);
    c1_f = exp(c_f(1+(j-1)*period:j*period, 1));
    y1_f = exp(y_f(1+(j-1)*period:j*period, 1));
    saving_f = log(yss*(css)/yss .* (c1_f) .* y1_f);
    lsavingf_1 = saving_f - hpfilter(saving_f, 1600);
    V=cov([lyh_1,lch_1,lxh_1,lnh_1,lkh_1,lzh_1,lnxh_1]);
    
    
    sdzHP=sqrt(diag(V)); % standard deviations
    std_y = sdzHP(1);
    c_percent = sdzHP(2)/std_y;
    n_percent = sdzHP(3)/std_y;
    k_percent = sdzHP(4)/std_y;
    x_percent = sdzHP(5)/std_y;
    z_percent = sdzHP(6)/std_y;
    standard_deviations(j,:) = sdzHP;
    std_dev_relative(j,:) = [1.0, c_percent, n_percent, k_percent, x_percent, z_percent];

    corHP =  V./(sqrt(diag(V))*sqrt(diag(V))'); % cross-correlations
    cross_corr_y(j,:) = xcorr(lyh_1, lyh_1, lag, 'coeff');
    cross_corr_k(j,:) = xcorr(lkh_1, lyh_1, lag, 'coeff');
    cross_corr_c(j,:) = xcorr(lch_1, lyh_1, lag, 'coeff');
    cross_corr_n(j,:) = xcorr(lnh_1, lyh_1, lag, 'coeff');
    cross_corr_x(j,:) = xcorr(lxh_1, lyh_1, lag, 'coeff');
    cross_corr_z(j,:) = xcorr(lzh_1, lyh_1, lag, 'coeff');
    cross_corr_nx(j,:) = xcorr(lnxh_1, lyh_1, lag, 'coeff');

    cross_corr_y_across_countries(j,:) = xcorr(lyh_1, lyf_1, lag, 'coeff');
    cross_corr_c_across_countries(j,:) = xcorr(lch_1, lcf_1, lag, 'coeff');
    cross_corr_saving(j,:) = xcorr(lsavingh_1, lsavingf_1, lag, 'coeff');
end

standard_deviations_mean = mean(standard_deviations * 100);
standard_deviations_std = std(standard_deviations * 100);
std_dev_relative_mean = mean(std_dev_relative);
std_dev_relative_std = std(std_dev_relative);
cross_corr_y_mean = mean(cross_corr_y);
cross_corr_y_std = std(cross_corr_y);
cross_corr_c_mean = mean(cross_corr_c);
cross_corr_c_std = std(cross_corr_c);
cross_corr_x_mean = mean(cross_corr_x);
cross_corr_x_std = std(cross_corr_x);
cross_corr_n_mean = mean(cross_corr_n);
cross_corr_n_std = std(cross_corr_n);
cross_corr_k_mean = mean(cross_corr_k);
cross_corr_k_std = std(cross_corr_k);
cross_corr_z_mean = mean(cross_corr_z);
cross_corr_z_std = std(cross_corr_z);
cross_corr_nx_mean = mean(cross_corr_nx);
cross_corr_nx_std = std(cross_corr_nx);

cross_corr_y_across_countries_mean = mean(cross_corr_y_across_countries);
cross_corr_y_across_countries_std = std(cross_corr_y_across_countries);
cross_corr_c_across_countries_mean = mean(cross_corr_c_across_countries);
cross_corr_c_across_countries_std = std(cross_corr_c_across_countries);
cross_corr_saving_mean = mean(cross_corr_saving);
cross_corr_saving_std = std(cross_corr_saving);

