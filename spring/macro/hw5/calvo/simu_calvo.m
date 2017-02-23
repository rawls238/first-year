%simu_calvo.m
%Produce artificial time series implied by the Calvo model with staggered price setting  in the nontraded sector  developed in the chapter entitled  ``Nominal Rigidity, Exchange Rates,  And Unemployment''
%of the book  ``Open Economy Macroeconomics,'' by Martín Uribe and Stephanie Schmitt-Grohé, 2016, Princeton University Press, forthcoming. 
%
%set variable filename to `calvo_opt' or  'calvo_peg' 
%
%Output: simu_calvo_opt.mat   or simu_calvo_peg.mat
%
%Calls: calvo_ss.m, calvo_opt_num_eval.m gx_hx.m
%
%Requirement: Matlab Toolbox for First-Order and Second-Order Accurate Solution to DSGE Models,  S. Schmitt-Grohe and M. Uribe <http://www.columbia.edu/~mu2166/>
%© Martín Uribe and Stephanie Schmitt-Grohé, January 29, 2016.
  
clear all
format compact

filename = 'calvo_peg'

T = 3e6; %length of model simulation 

Tburn = round(0.01*T);

calvo_ss %read parameter values and deterministic steady state

eval([filename '_num_eval']) %this .m script was created by running calvo_model.m.  
%in c:\data\uribe\book\dnwr\calvo

%The linearized equilibrium system is of the form
%y_t=gx x_t
%x_t+1 = hx x_t + ETASHOCK epsilon_t+1
[gx, hx, exitflag] = gx_hx(nfy, nfx, nfyp, nfxp);

nx = size(hx,1); %number of states

%Variance/Covariance matrix of innovation to state vector x_t
varshock = nETASHOCK*nETASHOCK';

%Position of variables in the control vector
eval(['load ' filename '.mat controlvar statevar'])
ncontrolvar = numel(controlvar);
for j=1:ncontrolvar
eval(['n' char(controlvar(j)) '='  num2str(find(controlvar==controlvar(j))) ';'])
end

%Position of variables in the state vector
nstatevar = numel(statevar);
for j=1:nstatevar
eval(['n' char(statevar(j)) '='  num2str(find(statevar==statevar(j))) ';'])
end

%create policy funciton for the vector 
%yT rstar h w epsi p pai cT tby dy
GX = zeros(10,nx);
GX(1:2,end-1:end) = eye(2);
GX(3,:) = gx(nh,:);
GX(4,:) = gx(nw,:);
GX(5,:) = gx(nepsi,:);
GX(6,np) = 1;
GX(7,:) = gx(npai,:);
GX(8,:) = gx(ncT,:);
GX(9,:) = gx(ntby,:);
GX(10,:) = gx(ndy,:);
GX(11,:) = gx(nrer,:);

x0 = zeros(nx,1); %initial condition

rng('default')
e = [randn(T+Tburn,1) zeros(T+Tburn, 1)];
[Y,X, e] = simu_1st(GX, hx, nETASHOCK, T+Tburn,x0,e);
Y = Y(Tburn+1:end,:);
X = X(Tburn+1:end,:);

YT = Y(:,1)*100;

RS = Y(:,2);
RS = exp(RS)*(1+RSTAR);
RS = (RS.^4-1)*100;
RS = RS-mean(RS);

H = Y(:,3)*100;

W = Y(:,4)*100;

EPSI = exp(Y(:,5));
EPSI = (EPSI.^4-1)*100;

%Note that Y(:,6) is p_{t-1}. Shift it one period forward to obtain p_t: 
P = Y(:,6)*100;
P=[P(2:end); P(end)]; 

PAI = ( exp(Y(:,7)).^4-1)*100;

CT = Y(:,8)*100;

TBY = Y(:,9)*100;

DY = Y(:,10) *100;

EPSIRER = exp(Y(:,11));
EPSIRER = (EPSIRER.^4-1)*100;

eval(['save simu_' filename '.mat']) 