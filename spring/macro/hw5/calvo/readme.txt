Simulating a Typical Crisis in an Two Sector Open Economy Model with Calvo Price Setting
%Chapter ``Nominal Rigidity, Exchange Rates,  And Unemployment''
%Book ``Open Economy Macroeconomics''
%Authors: © Martín Uribe and Stephanie Schmitt-Grohé
%Date:February 5, 2016
(1) This set of programs requires the matlab toolbox for first-order <http://www.columbia.edu/~mu2166/1st_order/1st_order.htm> and second-order <http://www.columbia.edu/~mu2166/2nd_order.htm> approximation of DSGE models produced by Schmitt-Grohe and Uribe (JEDC, 2004). 

(2) run the program calvo_model.m once for the peg economy and once for  the optimal exchange-rate economy. Inside that program, set the variable optimal_policy to 0 or 1 to indicate the  peg economy or optimal exchange rate economy. This will produce the files calvo_opt.mat, calvo_opt_num_eval.m, calvo_peg.mat, and calvo_peg_num_eval.m

(3) run the program simu_calvo.m. Set the variable filename to calvo_peg or to calvo_opt. 
This program produces artificial time series in the Calvo model. It produces the files simu_calvo_peg.mat or simu_calvo_opt.mat

(4)  Run the program typical_crisis.m. It produces the figure entitled ``Crisis Dynamics in the Calvo Model: The Role Of Exchange-Rate Policy.''