 % CH2 Solving the real business cycle model
 %
 % ECON 6140, Karel Mertens, Cornell University
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear all;
 close all;
 
 % Utility function: u(C,L)=log(C)+thetal*log (L)
 % Production function Y = AK^(1-alpha) (XN)^alpha

 % The Euler conditions are:
 
 %  1/C(t) = beta*gammax^(1-sigma)E_t[1 /C(t+1) * ((1-alpha)* A(t+1)(N(t+1)/K(t+1))^alpha+1-\delta)] 
 %  alpha * A(t)(K(t)/N(t))^(1-alpha) / C(t) = thetal/(1-N(t))
 %  C(t) + gammax*K(t+1)-(1-delta)K(t) = A(t)*K(t)^(1-alpha)*N(t)^(alpha)
  
 % Parameter values
 alfa   = 0.58;
 betas  = 0.988;
 delta  = 0.025;
 gammax = 1.004;
 yss    = 1;    % Output is normalized to 1
 nss    = 0.20; % Steady state hours are 0.20
 rho    = 0.9;
 
 % Steady state capital, consumption and investment:
 kss    = (1-alfa)*yss/(gammax/betas-1+delta);
 iss    = (delta+gammax-1)*kss;
 css    = yss-iss;
 
 % Note: the implied values of A and thetal are
 A      = yss/kss^(1-alfa)/nss^alfa;
 thetal = alfa * yss/nss/css*(1-nss);
 
 % Model Solution
 W1 = [(1-betas/gammax*(1-delta))*(alfa*(1-alfa)/(nss/(1-nss)+1-alfa)-alfa)...
      (1-betas/gammax*(1-delta))*(1+alfa/(nss/(1-nss)+1-alfa))...
      -1-(1-betas/gammax*(1-delta))*alfa/(nss/(1-nss)+1-alfa);
      0 1 0;
      -iss/yss*gammax/(delta+gammax-1) 0 0];
 W2 = [0 0 -1;
       0 rho 0;
       -(1-alfa+iss/yss*(1-delta)/(delta+gammax-1))-alfa*(1-alfa)/(nss/(1-nss)+1-alfa)...
       -1- alfa/(nss/(1-nss)+1-alfa)...
       css/yss+ alfa/(nss/(1-nss)+1-alfa)];
 W = inv(W1)*W2;
 
 phi21 = (-(W(1,1)-W(3,3))-sqrt((W(1,1)-W(3,3))^2+4*W(1,3)*W(3,1)))/2/W(1,3);
 phi11 = W(1,1) + W(1,3)*phi21;
 phi22 = (W(3,2)-phi21*W(1,2))/(phi21*W(1,3)+rho-W(3,3));
 phi12 = W(1,2)+W(1,3)*phi22;
 
 % Let s = [k a] and z = [y c i n]
 
 G = [phi11 phi12;
      0     rho];
  
 H = [1-alfa+alfa*(nss/(1-nss)+1-alfa)^(-1)*(-phi21+(1-alfa))...
      1 + alfa*(nss/(1-nss)+1-alfa)^(-1)*(1-phi22);
     phi21          phi22;
     gammax/(delta+gammax-1)*phi11-(1-delta)/(delta+gammax-1) gammax/(delta+gammax-1)*phi12;
     (nss/(1-nss)+1-alfa)^(-1)*(-phi21+(1-alfa)) (nss/(1-nss)+1-alfa)^(-1)*(1-phi22)];
 
%% Impulse responses
 s(:,1) = [0;1];
 x(1,1) =  1;
 for i=1:20
 s(:,i+1) = G*s(:,i);
 x(i+1,1) = gammax*x(i,1);
 end
 
 z = (H*s)';
 s = s';
 
 % Plot Impulse Responses
 plot(s(:,1),'-s','MarkerSize',3,'Color',[0,0.5,0])
 hold on 
 plot(z(:,2),'-+','MarkerSize',3,'Color',[0,0,0.5])
 hold on 
 plot(z(:,1),'-d','MarkerSize',3,'Color',[0,0,0])
 hold on 
 plot(z(:,4),'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
 hold on 
 plot(z(:,3),'-^','MarkerSize',3,'Color',[0.9, 0, 0])
 hold on
 plot(s(:,2),'--','MarkerSize',3)
 hold on 
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('k','c','y','n','i','a')
 saveas(gcf,'rbcmodelgraph1_high_rho','psc2')
%% Simulated Time Series
 clear s z
 sigmae = 0.01;
 e      = sigmae*randn(10000,1);
 s(:,1) = [0 e(1)];
 x(1,1) =  1;
 for i=1:length(e)
 s(:,i+1) = G*s(:,i)+[0;1]*e(i);
 x(i+1,1) = gammax*x(i,1);
 end
 z = (H*s)';
 s = s';
 
 sim_y_series = yss*exp(z(1:200,1)).*x(1:200);
 sim_c_series = css*exp(z(1:200,2)).*x(1:200);
 sim_n_series = nss*exp(z(1:200,4));
 sim_wage = sim_c_series + nss/(1-nss) * sim_n_series;
 sim_a_series = A * exp(s(1:200,2));
 solow_residual = log(sim_a_series) + alfa*log(x(1:200));
 %% Plot simulated series
 plot(sim_c_series,'MarkerSize',3,'Color',[0,0,0.5])
 hold on 
 plot(sim_y_series,'MarkerSize',3,'Color',[0,0,0])
 hold on 
 plot(sim_n_series,'MarkerSize',3,'Color',[0.5, 0.5, 0])
 hold on 
 plot(iss*exp(z(1:200,3)).*x(1:200),'MarkerSize',3,'Color',[0.9, 0, 0])
 ylabel('Levels')
 xlabel('Quarters')
 hold off
 legend('c','y','n','i')
 saveas(gcf,'rbcmodelgraph2','psc2')
 
%% Compute Population Moments Analytically
 [m,m] = size(G);
 Sigmae = [0 0; 0 sigmae^2];
 % Var matrix of the states
 Sigmas = reshape(inv(eye(m^2)-kron(G,G))*Sigmae(:),2,2);
 % Var matrix of the variables of interest
Sigmaz = H*Sigmas*H';

sdzLT = sqrt(diag(Sigmaz));     % standard deviations
corLT =  Sigmaz./(sqrt(diag(Sigmaz))*sqrt(diag(Sigmaz))'); % cross-correlations
autocovLT =   H*G*Sigmas*H';
autocorrLT =   autocovLT./(sqrt(diag(Sigmaz))*sqrt(diag(Sigmaz))'); %auto-correlations

%% Compute Population Moments by simulation and hp-filtering
for j = 1:10000/200
 v = cov([(log(yss*exp(z(1+(j-1)*200:j*200,1)).*x(1+(j-1)*200:j*200,:))...
     -hpfilter(log(yss*exp(z(1+(j-1)*200:j*200,1)).*x(1+(j-1)*200:j*200,:)),1600)),...
    (log(css*exp(z(1+(j-1)*200:j*200,2)).*x(1+(j-1)*200:j*200,:))...
    -hpfilter(log(css*exp(z(1+(j-1)*200:j*200,2)).*x(1+(j-1)*200:j*200,:)),1600)),...
    (log(nss*exp(z(1+(j-1)*200:j*200,4)))-hpfilter(log(nss*exp(z(1+(j-1)*200:j*200,4))),1600)),...
    (log(iss*exp(z(1+(j-1)*200:j*200,3)).*x(1+(j-1)*200:j*200,:))...
    -hpfilter(log(iss*exp(z(1+(j-1)*200:j*200,3)).*x(1+(j-1)*200:j*200,:)),1600))]);
sdzhp(:,j) = v(:);
z_1 = log(yss*exp(z(2+(j-1)*200:j*200,1)).*x(2+(j-1)*200:j*200,:));
z_2 = log(yss*exp(z(1+(j-1)*200:j*200-1,1)).*x(1+(j-1)*200:j*200-1,:));
acy(:,j) = corr((z_1-hpfilter(z_1, 1600)), (z_2 - hpfilter(z_2, 1600)));
end
V = reshape(mean(sdzhp')',4,4);

sdzHP=sqrt(diag(V)); % standard deviations
corHP =  V./(sqrt(diag(V))*sqrt(diag(V))'); % cross-correlations
autocorrHPy = mean(acy); % auto-correlation of output