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
  

filename = 'calvo_opt'

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

[gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx) ;
[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,nETASHOCK);

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


%Perturbation parameter
sig=1;

%Second derivative of V^{opt} with respect to sigma evaluated at
%(x_0,\sigma) = (x^{ss}, 0): 
V_opt_ss  =   gss(nv)  ;

[Ey,Ex] = unconditional_mean(gx, hx, gxx, hxx, gss, hss, nETASHOCK, 1);

E_of_V_opt = Ey(nv) ;%this is in deviations from steady state 


%Now consider the peg economy: 

filename = 'calvo_peg'

calvo_ss %read parameter values and deterministic steady state

eval([filename '_num_eval']) %this .m script was created by running calvo_model.m.  
%in c:\data\uribe\book\dnwr\calvo

%The linearized equilibrium system is of the form
%y_t=gx x_t
%x_t+1 = hx x_t + ETASHOCK epsilon_t+1
[gx, hx, exitflag] = gx_hx(nfy, nfx, nfyp, nfxp);

decomp = variance_decomposition(gx, hx, 1);
decomp(:,nc);

nx = size(hx,1); %number of states

%Variance/Covariance matrix of innovation to state vector x_t
varshock = nETASHOCK*nETASHOCK';

[gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx) ;
[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,nETASHOCK);

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




%Perturbation parameter
sig=1;
%Second derivative of V^{opt} with respect to sigma evaluated at
%(x_0,\sigma) = (x^{ss}, 0): 
V_peg_ss  =   gss(nv);%Second derivative
V_cons_peg_ss  =   gss(nv_cons) ; 
V_h_peg_ss  =   gss(nv_h) ;  

cwelfare_cost = (V_opt_ss  - V_cons_peg_ss - V_h_peg_ss  )/ ((1-SIGG)* v_cons + (1-BETTA)^(-1)) ;
cwelfare_cost = cwelfare_cost *sig^2/2*100

%Unconditional Welfare Cost: 
[Ey,Ex] = unconditional_mean(gx, hx, gxx, hxx, gss, hss, nETASHOCK, 1);

E_of_V_cons_peg = Ey(nv_cons) ;%this is in deviations from steady state 
E_of_V_h_peg = Ey(nv_h) ;%this is in deviations from steady state 
uwelfare_cost = ((E_of_V_opt-v)  - (E_of_V_cons_peg-v_cons) - (E_of_V_h_peg-v_h)  )/ ((1-SIGG)* v_cons + (1-BETTA)^(-1))  *100



