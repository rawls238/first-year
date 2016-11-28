% CH2 Solving the real business cycle model
 %
 % ECON 614, Karel Mertens, Cornell University
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 
 m = 1;
 s = xi_n_n - xi_c_n - tau_n_n;
 p_c = xi_c_c + xi_c_n * (xi_c_c - xi_n_c);
 q_c = xi_c_c * (xi_n_n - tau_n_n) - xi_c_n * xi_n_c;
 p_lambda = xi_c_n * tau_n_lambda;
 p_k = xi_c_n * tau_n_k;
 p_z = xi_c_n * tau_n_z;
 
 % Model Solution
 % The Linearized System is of the form
 % 
 % Fyp*Ey(t+1)+Fxp x(t+1)+Fy*y(t)+Fx*x(t)=0
 % 
 % where y are the controls and x are the states

num_countries = 2;
country_offset = 4;
state_keys = {'K', 'Z', 'lambda'};
state_vals = zeros(length(state_keys), 1);
for i=1:length(state_keys)
    state_vals(i) = i;
end
state_offset = length(state_vals);
control_keys = {'C', 'Y', 'N'};
control_vals = zeros(length(control_keys), 1);
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
    
C = var_offset(country_num, vars('C'));
K = var_offset(country_num, vars('K'));
lambda = var_offset(country_num, vars('lambda'));
Z = var_offset(country_num, vars('Z'));
N = var_offset(country_num, vars('N'));
Y = var_offset(country_num, vars('Y'));

%2. Inventory FOC
eqn = 1 + (country_num - 1) * country_offset;
Fy(eqn,C) = xi_c_c;
Fy(eqn,N) = xi_c_n;
Fyp(eqn,C) = RHS_flag * xi_c_c;
Fyp(eqn,N) = RHS_flag * (xi_c_n + (1 - beta) * tau_z_n);
Fxp(eqn,lambda) = RHS_flag * (1 - beta) * tau_z_lambda;
Fxp(eqn,K) = RHS_flag * (1 - beta) * tau_z_k;
Fxp(eqn,Z) = RHS_flag * (1 - beta) * tau_z_z;

%3. Euler
mid = 1 / (1 + (delta - 1) * beta);
eqn = 2 + (country_num - 1) * country_offset;
Fyp(eqn, C) = RHS_flag * (xi_c_c - mid * (delta - 1) * beta * xi_c_c);
Fyp(eqn, N) = RHS_flag * (xi_c_n + tau_k_n - mid * (delta - 1) * beta * xi_c_n);
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
Fy(eqn, N) = xi_n_n - tau_n_n - xi_c_n;
Fx(eqn, lambda) = RHS_flag * tau_n_lambda;
Fx(eqn, K) = RHS_flag * tau_n_k;
Fx(eqn, Z) = RHS_flag * tau_n_z;
end

C_h = var_offset(1, vars('C'));
C_f = var_offset(2, vars('C'));
K_h = var_offset(1, vars('K'));
K_f = var_offset(2, vars('K'));
lambda_h = var_offset(1, vars('lambda'));
lambda_f = var_offset(2, vars('lambda'));
Z_h = var_offset(1, vars('Z'));
Z_f = var_offset(2, vars('Z'));
Y_h = var_offset(1, vars('Y'));
Y_f = var_offset(2, vars('Y'));
N_h = var_offset(1, vars('N'));
N_f = var_offset(2, vars('N'));

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
Fy(eqn,Y_h) = RHS_flag * yss;
Fy(eqn,Y_f) = RHS_flag * yss;
Fy(eqn,C_h) = css;
Fy(eqn,C_f) = css;
Fx(eqn,K_h) = RHS_flag *(1 - delta) * kss;
Fx(eqn,K_f) = RHS_flag *(1 - delta) * kss;
Fxp(eqn,K_h) = RHS_flag * (-kss);
Fxp(eqn,K_f) = RHS_flag * (-kss);


At = [-Fxp -Fyp]; Bt = [Fx Fy];
[H,G]=solab(At,Bt,size(Fx,2));

%% Impulse responses
 shock(:,1) = [0;0;1;0;0;0];
 for i=1:20
    shock(:,i+1) = G*shock(:,i);
 end

z = (H*shock)';
shock = shock';

lambda_h = shock(:,3);
lambda_f = shock(:,6);
k_h = shock(:,1);
k_f = shock(:,4);
c_h = z(:,1);
c_f= z(:, 2);
n_h = z(:,5);
n_f = z(:,6);
y_h = z(:,3);
y_f = z(:,4);
z_h = shock(:,2);
z_f = shock(:,5);
i_h = zeros(19, 1);
i_f = zeros(19, 1);
x_h = zeros(19, 1);
x_f = zeros(19, 1);
nx_h = zeros(19, 1);
nx_f = zeros(19, 1);

a = xi_n_n - xi_c_n - tau_n_n;
b = 1/a;
n_h = b * ((xi_c_c - xi_n_c) * c_h + tau_n_lambda * lambda_h + tau_n_k * k_h + tau_n_z * z_h);
y_h = zita_lambda * lambda_h + zita_k * k_h + zita_n * n_h + zita_z * z_h;

for i=1:19
    x_h(i) = k_h(i+1) - (1 - delta) * k_h(i);
    x_f(i) = k_f(i+1) - (1 - delta) * k_f(i);
    i_h(i) = x_h(i) + z_h(i+1) - z_h(i);
    i_f(i) = x_f(i) + z_f(i+1) - z_f(i);
    nx_h(i) = y_h(i) - c_h(i) - i_h(i);
    nx_f(i) = y_f(i) - c_f(i) - i_f(i);
end

 % Plot Impulse Responses
plot(0:21, [0 c_h'],'-d','MarkerSize',3,'Color',[0,0,0])
hold on 
plot(0:21, [0 y_h'],'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
hold on 
plot(0:19, [0 x_h'],'-.','MarkerSize',3,'Color',[0,.5,0])
hold on
plot(0:21, [0 lambda_h'], '-x','MarkerSize',3,'Color',[0,0,.5])
hold on
plot(0:19, [0 nx_h'], '-y', 'MarkerSize', 3,'Color',[0.5,0.5,.5])
hold off
ylabel('Percent Deviations')
xlabel('Quarters')
title('Home response')
legend('c', 'y','i','lambda', 'nx');

plot(0:21, [0 c_f'],'-d','MarkerSize',3,'Color',[0,0,0])
hold on 
plot(0:21, [0 y_f'],'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
hold on 
plot(0:19, [0 i_f'],'-.','MarkerSize',3,'Color',[0,.5,0])
hold on
plot(0:21, [0 lambda_f'], '-x','MarkerSize',3,'Color',[0,0,.5])
hold on
plot(0:19, [0 nx_f'], '-y', 'MarkerSize', 3,'Color',[0.5,0.5,.5])
hold off
ylabel('Percent Deviations')
xlabel('Quarters')
title('Foreign response')
legend('c', 'y','i','lambda', 'nx');
% saveas(gcf,'foreign','psc2')

period = 100;
sim_num = 50;
total_periods = period * sim_num;
mu = [0,0];
sigma = [.00852^2,.00852^2 * .258; .00852^2 * .258, .00852^2];
e = mvnrnd(mu,sigma,total_periods+1);
s = [0;0;e(1,1);0;0;e(1,2)];
for i=1:total_periods+1
    s(:,i+1) = G*s(:,i) + [0;0;e(i,1);0;0;e(i,2)];
end
z = (H*s)';
s = s';

y_h = z(:,Y_h);
k_h = s(:,K_h);
c_h = z(:,C_h);
n_h = z(:,N_h);
z_h = s(:,Z_h);
y_f = z(:,Y_f);
k_f = s(:,K_f);
c_f = z(:,C_f);
n_f = z(:,N_f);
z_f = s(:,Z_f);

i_h = zeros(total_periods, 1);
i_f = zeros(total_periods, 1);
x_h = zeros(total_periods, 1);
x_f = zeros(total_periods, 1);
nx_h = zeros(total_periods, 1);
nx_f = zeros(total_periods, 1);
for i=1:total_periods
    x_h(i) = k_h(i+1) - (1 - delta) * k_h(i);
    x_f(i) = k_f(i+1) - (1 - delta) * k_f(i);
    i_h(i) = x_h(i) + z_h(i+1) - z_h(i);
    i_f(i) = x_f(i) + z_f(i+1) - z_f(i);
    nx_h(i) = y_h(i) - c_h(i) - i_h(i);
    nx_f(i) = y_f(i) - c_f(i) - i_f(i);
end

nxss = yss - css - xss;
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
    nxh_1 = log(nxss * exp(nx_h(1+(j-1)*period:j*period, 1)));
    lnxh_1 = nxh_1 - hpfilter(nxh_1, 1600);
    V=cov([lyh_1,lch_1,lnh_1,lkh_1,lxh_1,lzh_1,lnxh_1]);
end

sdzHP=sqrt(diag(V)); % standard deviations
y_percent = 1.0;
std_y = sdzHP(1);
c_percent = sdzHP(2)/std_y;
n_percent = sdzHP(3)/std_y;
k_percent = sdzHP(4)/std_y;
x_percent = sdzHP(5)/std_y;
z_percent = sdzHP(6)/std_y;
nx_percent = sdzHP(7)/std_y;

corHP =  V./(sqrt(diag(V))*sqrt(diag(V))'); % cross-correlations
cross_corr_y = xcorr(lyh_1, lyh_1, 5, 'coeff');
cross_corr_k = xcorr(lkh_1, lyh_1, 5, 'coeff');
cross_corr_c = xcorr(lch_1, lyh_1, 5, 'coeff');
cross_corr_n = xcorr(lnh_1, lyh_1, 5, 'coeff');
cross_corr_x = xcorr(lxh_1, lyh_1, 5, 'coeff');
cross_corr_z = xcorr(lzh_1, lyh_1, 5, 'coeff');
cross_corr_nx = xcorr(lnxh_1, lyh_1, 5, 'coeff');

cross_corr_y_across_countries = xcorr(lyh_1, lyf_1, 5, 'coeff');
cross_corr_c_across_countries = xcorr(lch_1, lcf_1, 5, 'coeff');