%EDEIR_RUN.M
%Compute numerically a first-order approximation, 
%second moments, and impulse responses  implied by  the 
%Small Open Economy Model With An External Debt-Elastic Interest Rate 
%as presented in chapter 4 of ``Open Economy Macroeconomics,'' 
%by Martin Uribe and Stephanie Schmitt Grohe, Princeton University Press 2017.
%© Martín Uribe and Stephanie  Stephanie Schmitt-Grohé, January 2013.

clear all

% Create the symbolic expressions for the model. 
%You need to run this only once. 
edeir_model %This program generates the file edeir_model_num_eval.m

%Compute the steady state by running edeir_ss.m
edeir_ss

%Evaluate f and its derivates at the steady state by running edeir_model_num_eval.m
edeir_model_num_eval; 
%this is a .M file produced by running 
%edeir_model.m

%First-stage policy functions
[gx,hx,exitflag]=gx_hx(nfy,nfx,nfyp,nfxp); 
%the function gx_hx.m is available online at http://www.columbia.edu/~mu2166/1st_order/1st_order.htm

%Variance/Covariance matrix of innovation to state vector x_t
varshock = nETASHOCK*nETASHOCK';

%standard deviations
[sigy0,sigx0]=mom(gx,hx,varshock);
stds = sqrt(diag(sigy0));

%correlations with output
corr_xy = sigy0(:,noutput)./stds(noutput)./stds;

%serial correlations
[sigy1,sigx1]=mom(gx,hx,varshock,1);
scorr = diag(sigy1)./diag(sigy0);

%make a table containing second moments
num_table = [stds*100  scorr corr_xy];

%From this table, select variables of interest (output, c, ivv, h, tby, cay)
disp('In This table:');
disp('Rows: y,c,i,h, tb/y,ca/y');
disp('Columns: std,  serial corr., corr with y,');
num_table1 = num_table([noutput nc nivv nh ntby  ncay],:);
disp(num_table1);

edeir_model_2;

edeir_ss_1;

edeir_model_1_num_eval;

[gx1,hx1,exitflag]=gx_hx(nfy,nfx,nfyp,nfxp); 

%Compute Impulse Responses
T = 11; %number of periods for impulse responses
%Give a unit innovation to TFP
x0 = zeros(nstate,1);
x0(end) = 0.01;
%Compute Impulse Response
[IR IRy IRx]=ir(gx,hx,x0,T);
[IR1 IR1y IR1x]=ir(gx1,hx1,x0,T);

%Plot Impulse responses
t=(0:T-1)';

subplot(3,2,1)
plot(t,IR(:,noutput)*100,t,IR1(:,noutput)*100)
title('Output')

subplot(3,2,2)
plot(t,IR(:,nc)*100,t,IR1(:,nc)*100)
title('Consumption')

subplot(3,2,3)
plot(t,IR(:,nivv)*100,t,IR1(:,nivv)*100)
title('Investment')

subplot(3,2,4)
plot(t,IR(:,nh)*100,t,IR1(:,nh)*100)
title('Hours')

subplot(3,2,5)
plot(t,IR(:,ntby)*100,t,IR1(:,ntby)*100)
title('Trade Balance / Output')

subplot(3,2,6)
plot(t,IRx(:,na)*100,t,IR1x(:,na)*100)
title('TFP Shock')
shg