function calvo_model
%Produce the symbolic first- order approximation to the equilibrium 
%conditions of the Calvo model with staggered price setting  in the nontraded sector as  developed in the chapter entitled  ``Nominal Rigidity, Exchange Rates,  And Unemployment''
%of the book ``Open Economy Macroeconomics,'' by Mart�n Uribe and Stephanie Schmitt-Groh�, Princeton University Press. 
%Run this program once  for each exchange-rate policy specification. 
%Set the variable optimal_policy to 0 or 1 for the case peg or optimal policy, respectively.  No need to re-run it 
%when parameter values change. 
%calvo_model.m computes a symbolic  log-linear approximation to the  function f, which defines  the DSGE model: 
%  E_t f(yp,y,xp,x) =0. 
%here, a p denotes next-period variables.  
%
%Output: Analytical expressions for f and its first derivatives as well as x and y. 
%The output is written to either calvo_peg_num_eval.m or calvo_opt_num_eval.m, depending on the value of optimal_policy, which can then be run for numerical evaluations
%
%Calls: anal_deriv.m and anal_deriv_print2f.m 
%
%� Mart�n Uribe and Stephanie Schmitt-Groh�, January 29, 2016 (BSGDAY)

optimal_policy = 2 %this indicator takes the value 1 if optimal exchange-rate policy and 0 if currency peg. Other policies can readily be considered. Here, we contemplate only these two. 

syms v vp cT cTp r rp rstar rstarp   p pp h hp la lap w wp pNtilde pNtildep  yN yNp s sp d dp  dy dyp paiN paiNp pvmc pvmcp pvmr pvmrp output outputp tby tbyp  c cp pai paip a2 a2p ry ryp intr intrp RY INT PAI eta etap
syms v_cons v_consp v_h v_hp devShock devShockp

syms a11 a12 a21 a22 eta11 eta12 eta13 eta21 eta22 eta23 eta31 eta32 eta33 

%%Note: all varaibles with a p at the end are dated t+1, excpt 
%s=s_{t-1}; sp = s_t
%p = p_{t-1}; pp = p_t

syms SIGG A XI ALFA VARPHI HBAR MU TAU THETA BETTA RSTAR PSSI DBAR 
syms yT yTp rstar rstarp epsi epsip 


%Utitlity function
U = (c^(1-SIGG) -1 ) /(1-SIGG) + VARPHI*log(HBAR - h) ;

ces =1; %this indicator takes the value 1 if the Armington aggregator is CES and 0 if it is Cobb-Dpuglas

if ces==1
c = (A * cT^(1-1/XI ) + (1-A) * yN^(1-1/XI))^(1/(1-1/XI));
elseif ces==0
c = cT^A*yN^(1-A);
else
error('The variable ces must take the value 0 or 1')
end

%Equilibrium conditions. The symbols e1, e2, ... denote equation 1, equation2, ...

e1 = yT + dp/(1+r) - cT - d;

e2 = -la +  diff(U, 'c') * diff(c, 'cT'); 

e3 = -pp + diff(c, 'yN')/ diff(c, 'cT'); 

e4 =  w*la   +  diff(U, 'h'); 

e5 = -yN +  (h/sp)^(ALFA); 

e6 = -sp + THETA * s * paiN^(MU/ALFA) +  (1-   THETA)* pNtilde^(-MU/ALFA); 

e7 = -1 + THETA* paiN^(MU-1) + (1-THETA)*pNtilde^(1-MU); 

e8  = -pp + p * paiN / epsi; 

e9 = pvmr - pvmc; 

e10 = - pvmr + (MU-1)/MU * yN * pp*pNtilde^(1-MU) + BETTA * THETA * lap /la * ( pNtilde/pNtildep / paiNp )^(1-MU) * pvmrp ; 

e11 = - pvmc + (1-TAU)/ ALFA * yN^(1/ALFA) * w * pNtilde^(-MU/ALFA) + BETTA * THETA * lap / la *(pNtilde/pNtildep / paiNp)^(-MU/ALFA) * pvmcp; 

e12 = -r + rstar + PSSI * (exp(dp-DBAR) -1);

e13 =  -la/(1+r) + BETTA * lap; 

if optimal_policy==1
e14 = paiN -1; 
filename = 'calvo_opt'
elseif optimal_policy==0
e14 = -epsi + devShock;
filename = 'calvo_peg'
elseif optimal_policy == 2
e14=-log(intr/INT)+1.5*log(pai/PAI)+0.125*log(ry/RY)+eta;
filename = 'calvo_taylor';
else
error('optimal_policy must take either the vlaue 1 or the value 0')
end
 %filename is a string. A file with this name and suffix _num_eval.m is created at the end of this program and contains symbolic expression for the function f and its derivatives fx fy fxp fyp 
e15 = [-output + yT + pp * yN];             % DEF real gdp in terms of tradable goods, pp relative price at t  

e16 = - tby  +  (yT - cT)/output;           % DEF tby

e17 = -a2p + diff(c,'yN');                  % DEF a2p: derivative of c wrt cN at t

e18 = -pai + paiN * a2 / a2p;               % DEF pai is CPI inflation (CPI=pN/a2), pai and paiN both GROSS, a2 at t-1

e19 = -[ log(yTp); log((1+rstarp)/(1+RSTAR))] + [a11 a12; a21 a22] * [log(yT); log((1+rstar)/(1+RSTAR))]; 

e20 = -dy + d/output/4;                     % DEF debt-to-ouput ratio 

e21 = [-v + U+BETTA*vp;
           -v_cons + (c^(1-SIGG) -1 ) /(1-SIGG)  + BETTA*v_consp;
           -v_h +  VARPHI*log(HBAR - h)   + BETTA*v_hp];  %value function

e22 = -(1+r)+intr/epsip;

e23 = etap;

e24 = -ry+(a2p/pp)*output;

e25 = devShockp;


%Create function f
f = eval([e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17;e18;e19;e20;e21;e22;e23;e24;e25]);

% Define the vector of controls in periods t and t+1, controlvar and controlvarp, and the vector of states in periods t and t+1, statevar and statevarp
%States to substitute from levels to logs
states_in_logs = [a2 s p yT];
states_in_logsp = [a2p sp pp yTp];
statevar = [d states_in_logs rstar eta devShock];
statevarp = [dp  states_in_logsp rstarp etap devShockp];
 
controls_in_logs = [cT h la w pNtilde yN paiN pvmc pvmr epsi pai output intr ry];
controls_in_logsp = [cTp hp lap wp pNtildep yNp paiNp pvmcp pvmrp epsip paip outputp intrp ryp];
controlvar = [controls_in_logs tby  r dy v v_cons v_h]; 
controlvarp = [controls_in_logsp tbyp rp  dyp vp v_consp v_hp]; 

%Number of states
ns = length(statevar);

%Make f a function of the logarithm of the state and control vector

%variables to substitute from levels to logs
variables_in_logs = transpose([states_in_logs, controls_in_logs, states_in_logsp, controls_in_logsp]);

%variable transformations
f = subs(f, variables_in_logs, exp(variables_in_logs));

approx = 1;

%Compute analytical derivatives of f
[fx,fxp,fy,fyp]=anal_deriv(f,statevar,controlvar,statevarp,controlvarp,approx);

%Make f and its derivatives a function of the level of its arguments rather than the log
f = subs(f, variables_in_logs, log(variables_in_logs));
fx = subs(fx, variables_in_logs, log(variables_in_logs));
fy = subs(fy, variables_in_logs, log(variables_in_logs));
fxp = subs(fxp, variables_in_logs, log(variables_in_logs));
fyp = subs(fyp, variables_in_logs, log(variables_in_logs));

%Symbolically evaluate f and its derivatives at the nonstochastic steady state (c=cp, etc.)
cu = transpose([statevar controlvar]);
cup = transpose([statevarp controlvarp]);

for subcb=1:2 %substitution must be run twice  in case the original system is a stochastic difference equation of order higher than one. For example, the current model features k_t, k_t+1, and k_t+2 and thus is a 2nd order difference equation

f = subs(f, cup,cu,0);
fx = subs(fx, cup,cu,0);
fy = subs(fy, cup,cu,0);
fxp = subs(fxp, cup,cu,0);
fyp = subs(fyp, cup,cu,0);
end

%Construct ETASHOCK matrix, which determines the var/cov of the forcing term of the system. Specifically, the state vector evolves over time according to 
%x_t+1 = hx x_t + ETASHOCK epsilon_t+1
ETASHOCK = sym(zeros(ns,4)); 
ETASHOCK(ns-3:end,:) =[eta11 eta12 eta13 0;eta21 eta22 eta23 0;eta31 eta32 eta33 0; 0 0 0 1];

%Print derivatives to file <filename>_num_eval.m'  for model evaluation
anal_deriv_print2f(filename,fx,fxp,fy,fyp,f,ETASHOCK);

eval(['save ' filename  '.mat'])