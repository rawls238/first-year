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

filename = 'calvo_taylor';

T = 3e6; %length of model simulation 

Tburn = round(0.01*T);

calvo_taylor_ss %read parameter values and deterministic steady state

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

%% Compute Impulse Responses
T = 20; %number of periods for impulse responses
%Give a unit innovation to TFP
x0 = zeros(nstatevar,1);
x0(end-1) = -0.01;
%Compute Impulse Response
[IR IRy IRx]=ir(gx,hx,x0,T);

%Plot Impulse responses
t=(0:T-1)';

subplot(3,2,1)
plot(t,IRy(:,nepsi)*100)
title('rate')
subplot(3,2,2)
plot(t,IRy(:,nv)*100)
title('welfare value function')

x0 = zeros(nstatevar,1);
x0(end) = 0.01;

filename = 'calvo_opt';

calvo_ss %read parameter values and deterministic steady state

eval([filename '_num_eval']) %this .m script was created by running calvo_model.m.  
%in c:\data\uribe\book\dnwr\calvo

%The linearized equilibrium system is of the form
%y_t=gx x_t
%x_t+1 = hx x_t + ETASHOCK epsilon_t+1
[gx_opt, hx_opt, exitflag] = gx_hx(nfy, nfx, nfyp, nfxp);


filename = 'calvo_peg';
eval([filename '_num_eval']) %this .m script was created by running calvo_model.m.  

[gx_peg, hx_peg, exitflag] = gx_hx(nfy, nfx, nfyp, nfxp);

T=10;
x0=x0(:);
pd=length(x0);
T1X = [gx_peg; eye(pd)];
T2X = [gx_opt; eye(pd)];
IR=[];
x=x0;
for t=1:T
    if t == 1
        IR(t,:) = (T1X*x);
        x = hx_peg * x;
    else
        IR(t,:) = (T2X*x);
        x = hx_opt * x;
    end
end

IRx = IR(:,end-pd+1:end);
IRy = IR(:,1:end-pd);

t=(0:T-1)';
f = figure();
subplot(3,2,1)
plot(t,IRy(:,nepsi)*100)
title('rate')
subplot(3,2,2)
plot(t,IRy(:,nv)*100)
title('welfare value cost')


