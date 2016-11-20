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
 triangle = (delta * (1 + beta + beta^2 + beta^3) / 4 * beta^4) - (1 - delta);
 D = (4 * triangle * beta^4 + (1 - delta) * (4 * beta^4 + 1))^-1;
 sc = css / yss;
 square = ((1 - beta) / (beta * sigma))^(v + v^2);
 B = (square - sigma)^-1;
 E = (B^-1 + sigma)^(-1/v);
 G = nss / 1 - nss;
 
 tmp = sigma * yss^v * zss^-v;
 xi_c_c = gamma * (1 - mu);
 xi_c_l = gamma - gamma * mu;
 xi_l_c = gamma * mu;
 xi_l_l = -1 + gamma - gamma * mu;
 zita_k = theta * (1 - tmp);
 zita_n = (1 - theta) * (1 - tmp);
 zita_z = tmp;
 zita_lambda = B / (1 + B*sigma);
 tau_k_k = (v + 1) * theta * (1 - tmp) - v * theta - 1;
 tau_k_n = (v + 1) * (1 - theta) * (1 - tmp) + v * (theta - 1);
 tau_k_z = (v + 1) * tmp;
 tau_n_z = tau_k_z;
 tau_n_k = (v + 1) * theta * (1 - tmp) - v * theta;
 tau_n_n = (v + 1) * (1 - theta) * (1 - tmp) + v * (theta - 1) - 1;
 tau_z_k = (v + 1) * theta * (1 - tmp);
 tau_z_n = (v + 1) * (1 - theta) * (1 - tmp);
 tau_z_z = (v + 1) * tmp - v - 1;
 tau_z_lambda = (1 + v) / (1 + sigma);
 tau_n_lambda = (1 - B * v * sigma) / (1 + B * sigma);
 tau_k_lambda = (1 - B * v * sigma) / (1 - B * sigma);
 
 
 H = -D*xi_c_l*G;
 M = tau_k_n * triangle / (triangle + 1 - delta) - xi_c_l * G;

 
 % Model Solution
 % The Linearized System is of the form
 % 
 % Fyp*Ey(t+1)+Fxp x(t+1)+Fy*y(t)+Fx*x(t)=0
 % 
 % where y are the controls and x are the states
Kt=1; Kt1=2; Kt2=3; Kt3=4; lambdat=5; lambdat1=6; lambdat2=7; lambdat3=8;
Zt=9;Zt1=10;Zt2=11;Zt3=12;
Ct=1;Ct1=2;Ct2=3;Ct3=4;
Nt=5;Nt1=6;Nt2=7;Nt3=8;

% Fx(eqn, Kt) has same index as Fxp(eqn, Kt) which corresponds to Kt and
% Kt+1
Fy=zeros(5,8);  Fx=zeros(5,12); 
Fyp=zeros(5,8); Fxp=zeros(5,12);

%1. Resource Constraint (equation 4)
eqn = 1;
Fy(eqn,Ct) = sc - zita_n*(1 / (tau_n_n - 1)); 
Fx(eqn,Kt) = -1 * 0.25*(1 - delta) * (1 - sc) / delta + zita_n*(tau_n_k / (tau_n_n - 1));
Fx(eqn,Kt1)  = 0.25*(1- sc);
Fx(eqn,Kt2)  = 0.25*(1-sc);
Fx(eqn,Kt3)  = 0.25*(1-sc);
Fxp(eqn,Kt3) = 0.25 * (1-sc) / delta;
Fx(eqn,lambdat) = -1 * zita_lambda + zita_n * tau_n_lambda / (tau_n_n - 1); 
Fx(eqn,Zt) = zita_n * tau_n_z / (tau_n_n - 1) - zita_z;

%2. Inventory equation (equation 3)
eqn = 2;
tmp = (1 - beta) * tau_z_n - xi_c_l * G;
Fy(eqn,Ct) = xi_c_c - xi_c_l * G * 1 / (tau_n_n - 1);
Fyp(eqn,Ct) = -1 * (xi_c_c - (tmp * (tau_n_n - 1)));
Fx(eqn,Kt) = xi_c_l * G * tau_n_k / (tau_n_n - 1) + tmp * tau_n_k / (tau_n_n - 1);
Fxp(eqn,Kt) =  -1 * (1 - beta) * tau_z_k;
Fx(eqn,Zt) = xi_c_l * G * tau_n_z / (tau_n_n - 1);
Fxp(eqn,Zt) = tmp * tau_n_z * (tau_n_n - 1) - (1 - beta) * tau_z_z;
Fx(eqn, lambdat) = xi_c_l * G * (tau_n_lambda / (tau_n_n - 1));
Fxp(eqn, lambdat) =  -1 * (1 - beta) * tau_z_lambda + tmp * tau_n_lambda / (tau_n_n - 1);

%3. Equation for n_hat
eqn = 3;
Fy(eqn,Ct) = -1;
Fy(eqn,Nt) = tau_n_n - 1;
Fx(eqn,Zt) = tau_n_z;
Fx(eqn,lambdat) = tau_n_lambda;
Fx(eqn,Kt) = tau_n_k;

% Euler equation (eq. 1)
eqn = 4;
Fy(eqn,Nt)=H;
Fy(eqn,Nt1) = H * beta * delta;
Fy(eqn,Nt2) = H * beta^2 * delta;
Fy(eqn,Nt3) = H * beta^3 * delta;
Fy(eqn,Ct) = D * xi_c_c;
Fy(eqn,Ct1) = D * xi_c_c * beta * delta;
Fy(eqn,Ct2) = D * xi_c_c * beta^2 * delta;
Fy(eqn,Ct3) = D * xi_c_c * beta^3 * delta;
Fyp(eqn,Ct3) = -1* (M / (tau_n_n - 1) + xi_c_c);
Fxp(eqn,Zt3) = -1 * (triangle / (triangle + 1 - delta) * tau_k_z - (M*tau_n_z) / (tau_n_n - 1));
Fxp(eqn,lambdat3) = -1 * (triangle / (triangle + 1 - delta) * tau_k_lambda - ((M * tau_n_lambda) / (tau_n_n - 1)));
Fxp(eqn,Kt3) = -1 * (triangle / (triangle + 1 - delta) * tau_k_k - ((M * tau_n_k) / (tau_n_n - 1)));

% Technology process
eqn = 5;
Fxp(eqn,lambdat)  = -1;
Fx(eqn,lambdat)  = rho;


At = [-Fxp -Fyp]; Bt = [Fx Fy];
[H,G]=solab(At,Bt,size(Fx,2));


 
 