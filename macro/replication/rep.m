 % Based on code by Karel Mertens for solving RBC models
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
 
 yss = 1;    % Output is normalized to 1
 zss = yss * ((1 - beta) / (beta * sigma))^(1/(-v - 1));
 kss_denom = 1 + delta * beta + delta * beta^2 + delta * beta^3 + (delta - 1) * beta^4;
 kss = (theta * (4 * beta^4) / kss_denom) * yss^(v+1) * ((yss)^(-v) - sigma * (zss)^(-v));
 nss_first_term = ((1 - theta) * (yss)^(v+1)*((yss)^(-v) - sigma*(zss)^(-v)))^(-1);
 nss_second_term = (1/mu - 1) * (yss - delta * kss);
 nss = (nss_first_term * nss_second_term + 1)^-1;
 css = yss - delta * kss;
 
 triangle = (delta * (1 + beta + beta^2 + beta^3) / (4 * beta^4)) - (1 - delta);
 D = (4 * triangle * beta^4 + (1 - delta) * (4 * beta^4 + 1))^-1;
 sc = css / yss;
 square = ((1 - beta) / (beta * sigma))^(v + v^2);
 B = (square - sigma)^-1;
 E = (B^-1 + sigma)^(-1/v);
 G = nss / (1 - nss);
 
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
 
 
 H = -D*xi_c_l*G;
 M = tau_k_n * triangle / (triangle + 1 - delta) - xi_c_l * G;

 RHS_flag = -1;
 % Model Solution
 % The Linearized System is of the form
 % 
 % Fyp*Ey(t+1)+Fxp x(t+1)+Fy*y(t)+Fx*x(t)=0
 % 
 % where y are the controls and x are the states
Kt=1; Kt1=2; Kt2=3; Kt3=4;
lambdat=5; lambdat1=6; lambdat2=7; lambdat3=8;

Ct=1;Ct1=2;Ct2=3;Ct3=4;
Zt=5;Zt1=6;Zt2=7;Zt3=8;


% Fx(eqn, Kt) has same index as Fxp(eqn, Kt) which corresponds to Kt and
% Kt+1
%matrix with equations from ou's version
Fy=zeros(16,8);  Fx=zeros(16,8);
Fyp=zeros(16,8); Fxp=zeros(16,8);


% matrix with equations from james' version
Fyy=zeros(16,8);  Fxx=zeros(16,8);
Fyyp=zeros(16,8); Fxxp=zeros(16,8);

m = 1 + delta * beta + delta * beta^2 + delta * beta^3;
s = xi_n_n - xi_c_n - tau_n_n;
p_c = xi_c_c + xi_c_n * (xi_c_c - xi_n_c);
q_c = xi_c_c * (xi_n_n - tau_n_n) - xi_c_n * xi_n_c;
p_lambda = xi_c_n * tau_n_lambda;
p_k = xi_c_n * tau_n_k;
p_z = xi_c_n * tau_n_z;

%1. Resource Constraint (equation 4)
eqn = 1;
Fy(eqn,Ct) = sc - zita_n*(1 / (tau_n_n - 1)); 
Fx(eqn,Kt) = 0.25*(1 - delta) * (1 - sc) / delta - zita_k + zita_n*(tau_n_k / (tau_n_n - 1));
Fx(eqn,Kt1)  = 0.25*(1- sc);
Fx(eqn,Kt2)  = 0.25*(1-sc);
Fx(eqn,Kt3)  = 0.25*(1-sc);
Fxp(eqn,Kt3) = RHS_flag * -0.25 * (1-sc) / delta;
Fx(eqn,lambdat) = RHS_flag * (zita_lambda + zita_n * tau_n_lambda / (tau_n_n - 1));
Fx(eqn,Zt) = zita_n * tau_n_z / (tau_n_n - 1) - zita_z;

Fyy(eqn,Ct) = 4 * (zita_n * (xi_c_c - xi_n_c) * yss - s * css);
Fxx(eqn,Kt) = 4 * yss * (zita_k * s + zita_n * tau_n_k) + (1 - delta) * s * kss;
Fxxp(eqn,Kt) = RHS_flag * delta * s * kss;
Fxxp(eqn,Kt1) = RHS_flag * delta * s * kss;
Fxxp(eqn,Kt2) = RHS_flag * delta * s * kss;
Fxxp(eqn,Kt3) = RHS_flag * s * kss;
Fxx(eqn, lambdat) = 4 * yss * (zita_lambda * s + zita_n * tau_n_lambda);
Fxx(eqn,Zt) = 4 * yss * (zita_z * s + zita_n * tau_n_z);


%2. Inventory equation (equation 3)
eqn = 2;
tmp = (1 - beta) * tau_z_n - xi_c_l * G;
Fy(eqn,Ct) = xi_c_c - xi_c_l * G * 1 / (tau_n_n - 1);
Fyp(eqn,Ct) = RHS_flag * (xi_c_c + (tmp * (tau_n_n - 1)));
Fx(eqn,Kt) = xi_c_l * G * tau_n_k / (tau_n_n - 1) + tmp * tau_n_k / (tau_n_n - 1);
Fxp(eqn,Kt) = RHS_flag * (1 - beta) * tau_z_k;
Fy(eqn,Zt) = xi_c_l * G * tau_n_z / (tau_n_n - 1);
Fyp(eqn,Zt) = RHS_flag * ((tmp * tau_n_z / (tau_n_n - 1)) + ((1 - beta) * tau_z_z));
Fx(eqn, lambdat) = xi_c_l * G * (tau_n_lambda / (tau_n_n - 1));
Fxp(eqn, lambdat) =  RHS_flag * ((1 - beta) * tau_z_lambda - tmp * tau_n_lambda / (tau_n_n - 1));

Fyy(eqn,Ct) = q_c;
Fyyp(eqn,Ct) = RHS_flag * (q_c + (1 - beta) * tau_z_n * (xi_c_c - xi_n_c));
Fxx(eqn,lambdat) = p_lambda;
Fxx(eqn,Kt) = p_k;
Fyy(eqn,Zt) = p_z;
Fxxp(eqn,lambdat) = RHS_flag * (xi_c_n * tau_n_lambda + (1 - beta) * (tau_z_n * tau_n_lambda + s * tau_z_lambda));
Fxxp(eqn, Kt) = RHS_flag * (xi_c_n * tau_n_k + (1 - beta) * (tau_z_n * tau_n_k + s * tau_z_k));
Fyyp(eqn, Zt) = RHS_flag * (xi_c_n * tau_n_z + (1 - beta) * (tau_z_n * tau_n_z + s * tau_z_z));

% Euler equation (eq. 1)
eqn = 3;
Fy(eqn,Ct) = (H / (tau_n_n - 1)) + D * xi_c_c;
Fy(eqn,Ct1) = (H * beta * delta / (tau_n_n - 1)) + D * xi_c_c * beta * delta;
Fy(eqn,Ct2) = (H * beta^2 * delta / (tau_n_n - 1)) + D * xi_c_c * beta^2 * delta;
Fy(eqn,Ct3) = (H * beta^3 * delta / (tau_n_n - 1)) + D * xi_c_c * beta^3 * delta;
Fyp(eqn,Ct3) = RHS_flag * (M / (tau_n_n - 1) + xi_c_c);
Fy(eqn,Zt) =  (H * tau_n_z / (tau_n_n - 1));
Fy(eqn,Zt1) =  (H * beta * delta * tau_n_z / (tau_n_n - 1));
Fy(eqn,Zt2) =  (H * beta^2 * delta * tau_n_z / (tau_n_n - 1));
Fy(eqn,Zt3) =  (H * beta^3 * delta * tau_n_z / (tau_n_n - 1));
Fyp(eqn,Zt3) = RHS_flag * (triangle / (triangle + 1 - delta) * tau_k_z - (M*tau_n_z) / (tau_n_n - 1));
Fx(eqn,lambdat) = (H * tau_n_lambda / (tau_n_n - 1));
Fx(eqn,lambdat1) = (H * beta * delta * tau_n_lambda / (tau_n_n - 1));
Fx(eqn,lambdat2) = (H * beta^2 * delta * tau_n_lambda / (tau_n_n - 1));
Fx(eqn,lambdat3) = (H * beta^3 * delta * tau_n_lambda / (tau_n_n - 1));
Fxp(eqn,lambdat3) = RHS_flag * (triangle / (triangle + 1 - delta) * tau_k_lambda - ((M * tau_n_lambda) / (tau_n_n - 1)));
Fx(eqn,Kt) = (H * tau_n_k / (tau_n_n - 1));
Fx(eqn,Kt1) = (H * beta * delta * tau_n_k / (tau_n_n - 1));
Fx(eqn,Kt2) = (H * beta^2 * delta * tau_n_k / (tau_n_n - 1));
Fx(eqn,Kt3) = (H * beta^3 * delta * tau_n_k / (tau_n_n - 1));
Fxp(eqn,Kt3) = RHS_flag * (triangle / (triangle + 1 - delta) * tau_k_k - ((M * tau_n_k) / (tau_n_n - 1)));

Fyy(eqn, Ct) = p_c;
Fyy(eqn, Ct1) = delta * beta * p_c;
Fyy(eqn, Ct2) = delta * beta^2 * p_c;
Fyy(eqn, Ct3) = delta * beta^3 * p_c;
Fyyp(eqn, Ct3) = RHS_flag * (m * q_c + tau_k_n * (xi_c_c - xi_n_c) * (m + (delta - 1) * beta^4));
Fxx(eqn, lambdat) = p_lambda;
Fxx(eqn, lambdat1) = delta * beta * p_lambda;
Fxx(eqn, lambdat2) = delta * beta^2 * p_lambda;
Fxx(eqn, lambdat3) = delta * beta^3 * p_lambda;
Fxxp(eqn, lambdat3) = RHS_flag * (m * xi_c_n * tau_n_lambda + (m + (delta - 1) * beta^4) * (tau_k_n * tau_n_lambda + s * tau_k_lambda));
Fyy(eqn, Zt) = p_z;
Fyy(eqn, Zt1) = delta * beta * p_z;
Fyy(eqn, Zt2) = delta * beta^2 * p_z;
Fyy(eqn, Zt3) = delta * beta^3 * p_z;
Fyyp(eqn, Zt3) = RHS_flag * (m * xi_c_n * tau_n_z + (m + (delta - 1) * beta^4) * (tau_k_n * tau_n_z + s * tau_k_z));
Fxx(eqn, Kt) = p_k;
Fxx(eqn, Kt1) = delta * beta * p_k;
Fxx(eqn, Kt2) = delta * beta^2 * p_k;
Fxx(eqn, Kt3) = delta * beta^3 * p_k;
Fxxp(eqn, Kt3) = RHS_flag * (m * xi_c_n * tau_n_k + (m + (delta - 1) * beta^4) * (tau_k_n * tau_n_k + s * tau_k_k));

% Technology process
eqn = 4;
Fxp(eqn,lambdat)  = -1;
Fx(eqn,lambdat)  = rho;

eqn = 5;
Fx(eqn,Kt1) = 1;
Fxp(eqn,Kt) = RHS_flag;

eqn = 6;
Fx(eqn,Kt2) = 1;
Fxp(eqn,Kt1) = RHS_flag;

eqn = 7;
Fx(eqn,Kt3) = 1;
Fxp(eqn,Kt2) = RHS_flag;

eqn = 8;
Fx(eqn, lambdat1) = 1;
Fxp(eqn, lambdat) = RHS_flag;

eqn = 9;
Fx(eqn, lambdat2) = 1;
Fxp(eqn, lambdat1) = RHS_flag;

eqn = 10;
Fx(eqn, lambdat3) = 1;
Fxp(eqn, lambdat2) = RHS_flag;

eqn = 11;
Fy(eqn, Zt1) = 1;
Fyp(eqn,Zt) = RHS_flag;

eqn = 12;
Fy(eqn, Zt2) = 1;
Fyp(eqn,Zt1) = RHS_flag;

eqn = 13;
Fy(eqn, Zt3) = 1;
Fyp(eqn,Zt2) = RHS_flag;

eqn = 14;
Fy(eqn, Ct1) = 1;
Fyp(eqn, Ct) = RHS_flag;

eqn = 15;
Fy(eqn, Ct2) = 1;
Fyp(eqn, Ct1) = RHS_flag;

eqn = 16;
Fy(eqn, Ct3) = 1;
Fyp(eqn, Ct2) = RHS_flag;


At = [-Fxp -Fyp]; Bt = [Fx Fy];
[H,G]=solab(At,Bt,size(Fx,2)); 
