%calvo_ss.m
%Steady state of the Two-sector open economy model with Staggered price 
%setting in the nontraded sector developed  in the  chapter entitled ``Nominal Rigidity, Exchange Rates,  And Unemployment,''  of  
%the book ``Open Economy Macroeconomics,'' 
%by Martin Uribe and Stephanie Schmitt-
%Grohe, Princeton University Press, forthcoming. 
%© Martín Uribe and Stephanie Schmitt-
%Grohé, 2016.


%Calibration
%Time unit is a quarter

SIGG = 2; %
ALFA =0.75;% labor elasticity of nontradable production 
MU = 6;%elasticity of substitution across varieties of intermediate nontraded goods. It implies a markup of price over marginal cost of 20 percent 

RSTAR = 0.0316;%interest rate
BETTA = 1/(1+RSTAR); %subjective discount factor 
%THETA =  0.7;% Degree of price rigidity. 

A = 0.26; %share parameter of Armington aggregator
XI = 1/SIGG; %elasticity of substitution between tradables and nontradables (parameter of the Armington aggregator). 

h = 1; %hours worked in steady state
HBAR  = 3; %time endowment
VARPHI =  (HBAR*ALFA*(1-A) - ALFA * (1-A)*h) / h^(1-ALFA+ALFA/XI); %parameter of subutility function of leisure 

s = 1; % price dispersion measure
yN =  h^ALFA; %nontraded output
cN = yN; %nontraded consumption
DBAR =  2.9014; %mean value of debt in under optimal exchange-rate policy 
%in the model with a constant subjective discount factor studied in earlier sections 
%(see preamble to this program).  It implies a  debt-to-output ratio of 
%21.6 percent per annum in the deterministic steady sate of the present model. 

PSSI =  3.35e-05; %parameter of debt-elastic interest rate premium function. Produced by running calibrate_pssi.m in c:\data\uribe\book\dnwr\calvo 
%this value matches the unconditional standard deviation of log(cT) of 18.5 percent implied by  the version of this model without a stationarity inducing device characterized in the section entitled ``External Crises and Exchange-Rate Policy: A Quantitative Analysis.''

d = DBAR;
yT = 1;
r = RSTAR;
rstar = RSTAR; 
cT = yT - r/(1+r) * d;
p = (1-A)/A * (cT/cN)^(1/XI);
paiN =1;% gross rate of inflation of nontradable goods

pai = paiN; %CPI inflation

la = A*cT^(-1/XI); 
w = VARPHI*(HBAR -h )^(-1) / la; 
pNtilde = 1; 

TAU = 1 /MU; % Labor Subsidy

epsi = paiN; 

pvmc = (1-TAU)/ALFA * w*h /(1-BETTA*THETA); 

pvmr = (MU-1)/MU*p*yN/(1-BETTA*THETA);

c = 1/(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1)); %consumption

v = (c^(1-SIGG)-1)/(1-SIGG) + VARPHI * log(HBAR-h);
v = v/(1-BETTA); 

%utility function split into consumption and labor part: 
v_cons = (c^(1-SIGG)-1)/(1-SIGG)/(1-BETTA);
v_h =  VARPHI * log(HBAR-h) / (1-BETTA); 

a2 = (1-A)*(yN/c)^(-1/XI); %partial derivative of Armington aggregator with respect to cN

output = yT+p*yN; %gdp in terms of tradables

dy = d/output/4;

tb = yT - cT;
tby = tb/output;

load tpm b omega
%produced by running tpm.m 
%in
%c:\data\uribe\default\mfiles\one_bond\two_shocks\

b = b';
a11 = b(1,1);
a12 = b(1,2);
a21 = b(2,1);
a22 = b(2,2);

eta11 = omega(1,1);
eta12 = omega(1,2);
eta21 = omega(2,1);
eta22 = omega(2,2);