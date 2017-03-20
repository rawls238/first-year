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
edeir_model_3 %This program generates the file edeir_model_num_eval.m

%Compute the steady state by running edeir_ss.m
edeir_ss_3

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

%correlations with investment
corr_xi = sigy0(:,nivv)./stds(nivv)./stds;
investment_corr = corr_xi(11)
