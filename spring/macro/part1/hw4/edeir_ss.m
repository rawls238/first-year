%edeir_ss.m
% computes the deterministic steady state  of the EDEIR open economy model  
% presented in Chapter 4 of ``Open Economy Macroeconomics,'' 
% by Martín Uribe and Stephanie Schmitt-Grohé (Princeton University Press 2017). 
% ©Martín Uribe and Stephanie  Stephanie Schmitt-Grohé, January 2013.

%Calibration 
%Time unit is a year
SIGG = 2; %mENDOZA
DELTA = 0.1; %depreciation rate
RSTAR = 0.04; %long-run interest rate
ALFA = 0.32; %F(k,h) = k^ALFA h^(1-ALFA)
OMEGA = 1.020894482180640; %Frisch ela st. from Mendoza 1991
DBAR =  0.867036295078822; %debt
PSSI = 0.000024591310326; %debt elasticity of interest rate
PHI = 0.034704484826501; %capital adjustment cost
RHO = 0.514357598235313; %persistence of TFP shock
ETATILDE = 0.004346886880154; %standard deviation of innovation to TFP shock

%OMEGA= 1.397907023457787;
%DBAR =   0.744410937075253;
%PSSI =   0.000098805777619;
%PHI =   0.027657543538344;
%RHO=   0.546977455306480;
%ETATILDE=   0.009228754480410;
%Implied parameters and steady state of endogenous variables
BETTA = 1/(1+RSTAR); %subjective discount factor

r = RSTAR; %interest rate

d = DBAR; %debt

KAPA = ((1/BETTA - (1-DELTA)) / ALFA)^(1/(ALFA-1)); %k/h

h = ((1-ALFA)*KAPA^ALFA)^(1/(OMEGA -1)); 

k = KAPA * h; %capital

kfu = k; 

output = KAPA^ALFA * h; %output

c = output-DELTA*k-RSTAR*DBAR;

ivv = DELTA * k; %investment

tb = output - ivv - c; %trade balance

tby = tb/output;

ca = -r*d+tb;

cay = ca/output;

a = 1; %technological factor

la = ((c - h^OMEGA/OMEGA))^(-SIGG); %marginal utility of wealth


