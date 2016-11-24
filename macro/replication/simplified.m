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
 rho = sum(A(1,:));
 
 % Steady state
 yss = 1;    % Output is normalized to 1
 zss = yss * ((1 - beta) / (beta * sigma))^(1/(-v - 1));
 kss_denom = 1 + (delta - 1)*beta
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
 zita_lambda = 1 - sigma * (yss)^v * (zss)^(-v);
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
K=1; lambda=2; 
C=1; Z=2;

RHS_flag = -1;

Fy=zeros(4,2);  Fx=zeros(4,2); 
Fyp=zeros(4,2); Fxp=zeros(4,2);

%1. Resource Constraint (equation 4)
eqn = 1;
Fy(eqn,C) = (zita_n * (xi_c_c - xi_n_c) * yss - s * css);
Fx(eqn,K) = yss * (zita_k * s + zita_n * tau_n_k) + (1 - delta) * s * kss;
Fxp(eqn,K) = RHS_flag * s * kss;
Fx(eqn, lambda) = yss * (zita_lambda * s + zita_n * tau_n_lambda);
Fx(eqn,Z) = yss * (zita_z * s + zita_n * tau_n_z);


%2. Inventory FOC
Fy(eqn,C) = q_c;
Fyp(eqn,C) = RHS_flag * (q_c + (1 - beta) * tau_z_n * (xi_c_c - xi_n_c));
Fx(eqn,lambda) = p_lambda;
Fx(eqn,K) = p_k;
Fx(eqn,Z) = p_z;
Fxp(eqn,lambda) = RHS_flag * (xi_c_n * tau_n_lambda + (1 - beta) * (tau_z_n * tau_n_lambda + s * tau_z_lambda));
Fxp(eqn, K) = RHS_flag * (xi_c_n * tau_n_k + (1 - beta) * (tau_z_n * tau_n_k + s * tau_z_k));
Fxp(eqn, Z) = RHS_flag * (xi_c_n * tau_n_z + (1 - beta) * (tau_z_n * tau_n_z + s * tau_z_z));
 
%3. Euler
eqn = 3;
Fy(eqn, C) = p_c;
Fyp(eqn, C) = RHS_flag * (m * q_c + tau_k_n * (xi_c_c - xi_n_c) * (m + (delta - 1) * beta));
Fx(eqn, lambda) = p_lambda;
Fyp(eqn, lambda) = RHS_flag * (m * xi_c_n * tau_n_lambda + (m + (delta - 1) * beta) * (tau_k_n * tau_n_lambda + s * tau_k_lambda));
Fx(eqn, Z) = p_z;
Fxp(eqn, Z) = RHS_flag * (m * xi_c_n * tau_n_z + (m + (delta - 1) * beta) * (tau_k_n * tau_n_z + s * tau_k_z));
Fx(eqn, K) = p_k;
Fxp(eqn, K) = RHS_flag * (m * xi_c_n * tau_n_k + (m + (delta - 1) * beta) * (tau_k_n * tau_n_k + s * tau_k_k));

%6. Technology process
eqn = 4;
Fxp(eqn,lambda)  = RHS_flag;
Fx(eqn,lambda)  = rho;


At = [-Fxp -Fyp]; Bt = [Fx Fy];
[H,G]=solab(At,Bt,size(Fx,2));
