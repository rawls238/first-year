 % CH2 Transition dynamics in the neoclassical growth model
 %
 % ECON 6140, Karel Mertens, Cornell University
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear all;
 close all;
 %%%%%%%%%%%%%%%%%%%%%%
 % 1) No-growth Model
 %%%%%%%%%%%%%%%%%%%%%%
 % Utility function: u(C,L)=log(C)+thetal * log (L)
 % Production function Y = K^(1-alpha) N^alpha

 % The Euler conditions are:
 
 %  1/C(t) = beta /C(t+1) * ((1-alpha)* (N(t+1)/K(t+1))^alpha+1-\delta) 
 %  alpha * (K(t)/N(t))^(1-alpha) / C(t) = thetal/(1-N(t))
 %  C(t) + K(t+1)-(1-delta)K(t) = A*K(t)^(1-alpha)*N(t)^(alpha)
 
 % Parameter values
 alfa   = 0.58;
 beta   = 0.988;
 delta  = 0.025;
 yss    = 1;    % Output is normalized to 1
 nss    = 0.20; % Steady state hours are 0.20
  
 % Steady state capital, consumption and investment:
 kss    = (1-alfa)*yss/(1/beta-1+delta);
 iss    = delta*kss;
 css    = yss-iss;
 
 
 % Note: the implied values of A and thetal are
 A      = yss/kss^(1-alfa)/nss^alfa;
 thetal = alfa * yss/nss/css*(1-nss);
 
 % Model Solution
 W1 = [(1-beta*(1-delta))*(alfa*(1-alfa)/(nss/(1-nss)+1-alfa) - alfa)...
     -1-(1-beta*(1-delta))*alfa/(nss/(1-nss)+1-alfa);
      -iss/yss/delta 0];
 W2 = [0 -1;
      -(1-alfa+iss/yss*(1-delta)/delta)-alfa*(1-alfa)/(nss/(1-nss)+1-alfa)...
      css/yss+ alfa/(nss/(1-nss)+1-alfa)];
 W = inv(W1)*W2;
 
 phi2 = (-(W(1,1)-W(2,2))-sqrt((W(1,1)-W(2,2))^2+4*W(1,2)*W(2,1)))/2/W(1,2);
 phi1 = W(1,1) + W(1,2)*phi2;
 
 % Compute the transition path
 k(1,1) = -1;
 for i=1:200
 k(i+1,1) = phi1*k(i,1);
 end
 c = phi2*k;
 n = (nss/(1-nss)+1-alfa)^(-1)*(-phi2+(1-alfa))*k;
 y = (1-alfa)*k+alfa*n;
 i = (1/delta*phi1-(1-delta)/delta)*k;
 % Plot Transition Paths
 plot(k,'-s','MarkerSize',3,'Color',[0,0.5,0])
 hold on 
 plot(c,'-+','MarkerSize',3,'Color',[0,0,0.5])
 hold on 
 plot(y,'-d','MarkerSize',3,'Color',[0,0,0])
 hold on 
 plot(n,'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
 hold on 
 plot(i,'-^','MarkerSize',3,'Color',[0.9, 0, 0])
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('k','c','y','n','i')
 saveas(gcf,'growthmodelgraph1','psc2') 

 % Plot Transition Paths 2
 %  plot(kss*exp(k).*x,'k-s')
 %  hold on 
 plot(css*exp(c),'-+','MarkerSize',3,'Color',[0,0,0.5])
 hold on 
 plot(yss*exp(y),'-d','MarkerSize',3,'Color',[0,0,0])
 hold on 
 plot(nss*exp(n),'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
 hold on
 plot(iss*exp(i),'-^','MarkerSize',3,'Color',[0.9, 0, 0])
 ylabel('Levels')
 xlabel('Quarters')
 hold off
 legend('c','y','n','i')
 saveas(gcf,'growthmodelgraph2','psc2')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 2) Model with Steady State Growth
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear all
 % Utility function: u(C,L)=log(C)+thetal * log (L)
 % Production function Y = A K^(1-alpha) (XN)^alpha

 % The Euler conditions are:
 
 %  1/c(t) = beta*gammax^(1-sigma)/c(t+1) * ((1-alpha)* A(N(t+1)/k(t+1))^alpha+1-\delta) 
 %  alpha * A(k(t)/N(t))^(1-alpha) / c(t) = thetal/(1-N(t))
 %  c(t) + gammax*k(t+1)-(1-delta)K(t) = A*k(t)^(1-alpha)*N(t)^(alpha)
 
 % Parameter values
 alfa   = 0.58;
 betas  = 0.988;
 delta  = 0.025;
 gammax = 1.004;
 yss    = 1;    % Output is normalized to 1
 nss    = 0.20; % Steady state hours are 0.20
  
 % Steady state capital, consumption and investment:
 kss    = (1-alfa)*yss/(gammax/betas-1+delta);
 iss    = (delta+gammax-1)*kss;
 css    = yss-iss;
 
 % Note: the implied values of A and thetal are
 A      = yss/kss^(1-alfa)/nss^alfa;
 thetal = alfa * yss/nss/css*(1-nss);
 
 % Model Solution
 W1 = [(1-betas/gammax*(1-delta))*(alfa*(1-alfa)/(nss/(1-nss)+1-alfa)-alfa)...
     -1-(1-betas/gammax*(1-delta))*alfa/(nss/(1-nss)+1-alfa);
      -iss/yss*gammax/(delta+gammax-1) 0 ];
 W2 = [0 -1;
       -(1-alfa+iss/yss*(1-delta)/(delta+gammax-1))-alfa*(1-alfa)/(nss/(1-nss)+1-alfa)...
       css/yss+ alfa/(nss/(1-nss)+1-alfa)];
 W = inv(W1)*W2;
 
 phi2 = (-(W(1,1)-W(2,2))-sqrt((W(1,1)-W(2,2))^2+4*W(1,2)*W(2,1)))/2/W(1,2);
 phi1 = W(1,1) + W(1,2)*phi2;
 
 % Compute the transition path
 k(1,1) = -1;
 x(1,1) =  1;
 for i=1:200
 k(i+1,1) = phi1*k(i,1);
 x(i+1,1) = gammax*x(i,1);
 end
 c = phi2*k;
 n = (nss/(1-nss)+1-alfa)^(-1)*(-phi2+(1-alfa))*k;
 y = (1-alfa)*k+alfa*n;
 i = (gammax/(delta+gammax-1)*phi1-(1-delta)/(delta+gammax-1))*k;
 
 % Plot Transition Paths 1
 plot(k,'-s','MarkerSize',3,'Color',[0,0.5,0])
 hold on 
 plot(c,'-+','MarkerSize',3,'Color',[0,0,0.5])
 hold on 
 plot(y,'-d','MarkerSize',3,'Color',[0,0,0])
 hold on 
 plot(n,'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
 hold on 
 plot(i,'-^','MarkerSize',3,'Color',[0.9, 0, 0])
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('k','c','y','n','i')
 saveas(gcf,'growthmodelgraph3','psc2')
 
 
 
 % Plot Transition Paths 2
 %  plot(kss*exp(k).*x,'k-s')
 %  hold on 
 plot(css*exp(c).*x,'-+','MarkerSize',3,'Color',[0,0,0.5])
 hold on 
 plot(yss*exp(y).*x,'-d','MarkerSize',3,'Color',[0,0,0])
 hold on 
 plot(nss*exp(n),'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
 hold on 
 plot(iss*exp(i).*x,'-^','MarkerSize',3,'Color',[0.9, 0, 0])
 ylabel('Levels')
 xlabel('Quarters')
 hold off
 legend('c','y','n','i')
 saveas(gcf,'growthmodelgraph4','psc2')
 