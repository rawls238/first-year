 % CH2 Transition dynamics in the neoclassical growth model
 %
 % ECON 6140, Karel Mertens, Cornell University
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear all;
 close all;
 
 % Parameter values
 alfa   = 0.58;
 betas  = 0.988;
 delta  = 0.025;
 gammaa = 1.001;
 gammav = 1.002;
 sigma = 2.0;
 theta = 0.25;
 gammas = gammaa^(1/alfa) * gammav^((1-alfa)/alfa);
 yss    = 1;    % Output is normalized to 1
 nss    = 0.20; % Steady state hours are 0.20
  
 % Steady state capital, consumption and investment:
 beta_tilda = (1/betas) * gammas^(theta*(1-alfa));
 delta_tilda = (1 - delta) / (gammav * gammas);
 kss    = (1-alfa)*yss/(beta_tilda - delta_tilda);
 css    = yss - kss + kss*(1-delta)/(gammav*gammas);
 
 % Model Solution
 first = yss / (yss + kss);
 second = kss / (kss + yss);
 third = css / (css + (kss/(gammas * gammav)));
 fourth = (kss/(gammas * gammav)) / (css + (kss/(gammas * gammav)));
 k_hat_coeff = first*(1-alfa) - third*(1-alfa) - fourth;
 n_hat_coeff = alfa*(third - first);
 first_1 = (1 - sigma)*theta - 1;
 second_1 = (1 - sigma)*(1 - theta) * nss / (1 - nss);
 third_1 = (1 - alfa) * (yss / kss) / ((1-alfa)*(yss / kss) - (1-delta)/(gammas*gammav));
 k_hat = -1 * (first_1 * alfa + third_1 * alfa);
 nt_1 = (alfa - 1) + 1/(1-nss) + second_1;
 nt = first_1*(1 - alfa - (1 / (1 - nss))) - third_1 * alfa - second_1;
 W1 = [-1*second 0; 
       -1*alfa*(third_1 - first_1) nt_1];
 W2 = [k_hat_coeff n_hat_coeff;
       k_hat nt];
 W = inv(W1)*W2;
 
 phi2 = (-(W(1,1)-W(2,2))-sqrt((W(1,1)-W(2,2))^2+4*W(1,2)*W(2,1)))/2/W(1,2);
 phi1 = W(1,1) + W(1,2)*phi2;
 
 % Compute the transition path
 k(1,1) = -1;
 S(1,1) =  1;
 for i=1:500
 k(i+1,1) = phi1*k(i,1);
 S(i+1,1) = gammas * S(i, 1);
 end
 n = phi2*k;
 y = (1-alfa)*k+alfa*n;
 c = y - 1/(1-nss)*n;
 x = y - c;
 
 % Plot Transition Paths 1
 plot(k,'-s','MarkerSize',3,'Color',[0,0.5,0])
 hold on 
 plot(c,'-+','MarkerSize',3,'Color',[0,0,0.5])
 hold on 
 plot(y,'-d','MarkerSize',3,'Color',[0,0,0])
 hold on 
 plot(n,'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
 hold on
 plot(x,'-^','MarkerSize',3,'Color',[0.9, 0, 0])
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('k','c','y','n','x')
 saveas(gcf,'growthmodelgraph3','psc2')
 
 
 
 xss = yss - css;
 % Plot Transition Paths 2
 plot(css*exp(c).*S,'-+','MarkerSize',3,'Color',[0,0,0.5])
 hold on 
 plot(yss*exp(y).*S,'-d','MarkerSize',3,'Color',[0,0,0])
 hold on 
 plot(nss*exp(n),'-^','MarkerSize',3,'Color',[0.5, 0.5, 0])
 hold on 
 plot(xss*exp(x).*S,'-^','MarkerSize',3,'Color',[0.9, 0, 0])
 ylabel('Levels')
 xlabel('Quarters')
 hold off
 legend('c','y','n','x')
 saveas(gcf,'growthmodelgraph4','psc2')
 