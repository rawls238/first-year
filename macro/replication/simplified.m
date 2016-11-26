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
 kss = (theta * beta / kss_denom) * yss^(v+1) * ((yss)^(-v) - sigma * (zss)^(-v));
 nss_first_term = ((1 - theta) * (yss)^(v+1)*((yss)^(-v) - sigma*(zss)^(-v)))^(-1);
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
country_offset = 2;
state_keys = {'K', 'Z', 'lambda'};
state_vals = [1, 2, 3];
state_offset = length(state_vals);
control_keys = {'C'};
control_vals = [1];
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

%2. Inventory FOC
eqn = 1 + (country_num - 1) * country_offset;
Fy(eqn,C) = q_c;
Fyp(eqn,C) = RHS_flag * (q_c + (1 - beta) * tau_z_n * (xi_c_c - xi_n_c));
Fx(eqn,lambda) = p_lambda;
Fxp(eqn,lambda) = RHS_flag * (xi_c_n * tau_n_lambda + (1 - beta) * (tau_z_n * tau_n_lambda + s * tau_z_lambda));
Fx(eqn,K) = p_k;
Fxp(eqn, K) = RHS_flag * (xi_c_n * tau_n_k + (1 - beta) * (tau_z_n * tau_n_k + s * tau_z_k));
Fx(eqn,Z) = p_z;
Fxp(eqn, Z) = RHS_flag * (xi_c_n * tau_n_z + (1 - beta) * (tau_z_n * tau_n_z + s * tau_z_z));
 
%3. Euler
eqn = 2 + (country_num - 1) * country_offset;
Fy(eqn, C) = p_c;
Fyp(eqn, C) = RHS_flag * (m * q_c + tau_k_n * (xi_c_c - xi_n_c) * (m + (delta - 1) * beta));
Fx(eqn, lambda) = p_lambda;
Fxp(eqn, lambda) = RHS_flag * (m * xi_c_n * tau_n_lambda + (m + (delta - 1) * beta) * (tau_k_n * tau_n_lambda + s * tau_k_lambda));
Fx(eqn, Z) = p_z;
Fxp(eqn, Z) = RHS_flag * (m * xi_c_n * tau_n_z + (m + (delta - 1) * beta) * (tau_k_n * tau_n_z + s * tau_k_z));
Fx(eqn, K) = p_k;
Fxp(eqn, K) = RHS_flag * (m * xi_c_n * tau_n_k + (m + (delta - 1) * beta) * (tau_k_n * tau_n_k + s * tau_k_k));
end

C_h = var_offset(1, vars('C'));
C_f = var_offset(2, vars('C'));
K_h = var_offset(1, vars('K'));
K_f = var_offset(2, vars('K'));
lambda_h = var_offset(1, vars('lambda'));
lambda_f = var_offset(2, vars('lambda'));
Z_h = var_offset(1, vars('Z'));
Z_f = var_offset(2, vars('Z'));

eqn = 5;
Fxp(eqn,lambda_h)  = RHS_flag;
Fx(eqn,lambda_h)  = A(1, 1);
Fx(eqn,lambda_f) = A(1, 2);


eqn = 6;
Fxp(eqn,lambda_f) = RHS_flag;
Fx(eqn,lambda_h) = A(2, 1);
Fx(eqn,lambda_f) = A(2, 2);

%1. Marginal Utility eq
eqn = 7;
Fy(eqn,C_h) = q_c;
Fx(eqn,lambda_h) = p_lambda;
Fx(eqn,K_h) = p_k;
Fx(eqn,Z_h) = p_z;
Fy(eqn,C_f) = RHS_flag * q_c;
Fx(eqn,lambda_f) = RHS_flag * p_lambda;
Fx(eqn,K_f) = RHS_flag * p_k;
Fx(eqn,Z_f) = RHS_flag * p_z;

%Resource constraint
eqn = 8;
Fxp(eqn,K_h) = RHS_flag * s * kss;
Fxp(eqn,K_f) = RHS_flag * s * kss;
Fy(eqn,C_h) = (zita_n * (xi_c_c - xi_n_c) * yss - s * css);
Fy(eqn,C_f) = (zita_n * (xi_c_c - xi_n_c) * yss - s * css);
Fx(eqn,K_h) = yss * (zita_k * s + zita_n * tau_n_k) + (1 - delta) * s * kss;
Fx(eqn,K_f) = yss * (zita_k * s + zita_n * tau_n_k) + (1 - delta) * s * kss;
Fx(eqn, lambda_h) = yss * (zita_lambda * s + zita_n * tau_n_lambda);
Fx(eqn, lambda_f) = yss * (zita_lambda * s + zita_n * tau_n_lambda);
Fx(eqn,Z_h) = yss * (zita_z * s + zita_n * tau_n_z);
Fx(eqn,Z_f) = yss * (zita_z * s + zita_n * tau_n_z);

At = [-Fxp -Fyp]; Bt = [Fx Fy];
[H,G]=solab(At,Bt,size(Fx,2));

%% Impulse responses
 shock(:,1) = [1;1;1;1;1;1];
 for i=1:100
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
z_h = shock(:,2);
z_f = shock(:,5);

i_h = zeros(100);
i_f = zeros(100);
for i=1:99
    i_h = k_h(i+1) + (1 - delta) * k_h(i);
    i_f = k_f(i+1) + (1 - delta) * k_f(i);
end

%y_h = ...
%nx_h = y_h - c..

 % Plot Impulse Responses
plot(c_h,'-d','MarkerSize',3,'Color',[0,0,0])
hold on 
plot(i_h,'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
hold off
ylabel('Percent Deviations')
xlabel('Quarters')
legend('c', 'i');
saveas(gcf,'home','psc2')

plot(c_f,'-d','MarkerSize',3,'Color',[0.9, 0, 0])
hold on
plot(i_f,'-^','MarkerSize',3)
hold on
hold off
ylabel('Percent Deviations')
xlabel('Quarters')
legend('c', 'i');
saveas(gcf,'foreign','psc2')