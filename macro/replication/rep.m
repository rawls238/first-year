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
 yss = 1;    % Output is normalized to 1
 zss = yss * ((1 - beta) / (beta * sigma))^(1/(-v - 1));
 kss_denom = 1 + delta * beta + delta * beta^2 + delta * beta^3 + (delta - 1) * beta^4;
 kss = ((4 * beta^4) / kss_denom) * yss^(v+1) * ((yss)^(-v) - sigma * (zss)^(-v));
 nss_first_term = ((1 - theta) * (yss)^(v+1)*((yss)^(-v) - sigma*(zss)^(-v)))^(-1);
 nss_second_term = (1 - 1/mu) * (yss - delta * kss);
 nss = (nss_first_term * nss_second_term + 1)^-1;
 css = yss - delta * kss;
 
 tmp = sigma * yss^v * zss^-v;
 pi_c_c = gamma * mu - 1;
 pi_c_n = gamma * (mu - 1) * nss / (nss - 1);
 pi_n_c = gamma * mu;
 pi_n_n = (nss / (nss - 1)) * (gamma * (1 - mu) - 1);
 zita_k = theta * (1 - tmp);
 zita_n = (1 - theta) * (1 - tmp);
 zita_z = tmp;
 tau_k_k = (v + 1) * theta * (1 - tmp) - v * theta - 1;
 tau_k_n = (v + 1) * (1 - theta) * (1 - tmp) + v * (theta - 1);
 tau_k_z = (v + 1) * tmp;
 tau_n_z = tau_k_z;
 tau_n_k = (v + 1) * theta * (1 - tmp) - v * theta;
 tau_n_n = (v + 1) * (1 - theta) * (1 - tmp) + v * (theta - 1) - 1;
 tau_z_k = (v + 1) * theta * (1 - tmp);
 tau_z_n = (v + 1) * (1 - theta) * (1 - tmp);
 tau_z_z = (v + 1) * tmp - v - 1;
 
 % Model Solution
 % The Linearized System is of the form
 % 
 % Fyp*Ey(t+1)+Fxp x(t+1)+Fy*y(t)+Fx*x(t)=0
 % 
 % where y are the controls and x are the states
K=1; a=2; 
Y=1; C=2; I=3; N=4; 

Fy=zeros(6,4);  Fx=zeros(6,2); 
Fyp=zeros(6,4); Fxp=zeros(6,2);

%1. Consumption
eqn = 1;
Fy(eqn,C)=1; 
Fyp(eqn,C)=-1;
Fyp(eqn,N) = alfa*(1-betas/gammax*(1-delta));
Fxp(eqn,a)=(1-betas/gammax*(1-delta));
Fxp(eqn,K) = -alfa*(1-betas/gammax*(1-delta));

%2. Labor euler equation
eqn = 2;
Fy(eqn,N) = -(1-alfa)-nss/(1-nss);
Fy(eqn,C) = -1;
Fx(eqn,a) = 1;
Fx(eqn,K) = (1-alfa);
 
%3. Resource constraint
eqn = 3;
Fy(eqn,C)  = css/yss;
Fy(eqn,N)  = - alfa;
Fx(eqn,K)  = - 1+alfa-iss/yss*(1-delta)/(delta+gammax-1);
Fx(eqn,a)  = - 1;
Fxp(eqn,K) = iss/yss*gammax/(delta+gammax-1);

 
%4. Production function
eqn = 4;
Fy(eqn,Y)  = -1;
Fy(eqn,N)  = alfa;
Fx(eqn,K)  = (1-alfa);
Fx(eqn,a)  = 1;

%5. Capital Accumulation
eqn = 5;
Fy(eqn,I)  = -1;
Fx(eqn,K)  = - (1-delta)/(delta+gammax-1);
Fxp(eqn,K) = gammax/(delta+gammax-1);

%6. Technology process
eqn = 6;
Fxp(eqn,a)  = -1;
Fx(eqn,a)  = rho;


At = [-Fxp -Fyp]; Bt = [Fx Fy];
[H,G]=solab(At,Bt,size(Fx,2));


 
 