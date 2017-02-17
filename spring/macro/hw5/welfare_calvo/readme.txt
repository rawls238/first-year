%Purpose: Compute the welfare cost of a currency peg vis-a-vis the optimal exchange rate policy in  a two sector open economy model with Calvo-type price stickiness in the nontraded sector. 
%Chapter ``Nominal Rigidity, Exchange Rates,  And Unemployment''
%Book ``Open Economy Macroeconomics''
%Authors: Martín Uribe and Stephanie Schmitt-Grohé
%© Martín Uribe and Stephanie Schmitt-Grohé, February  2016 

(1) This set of programs requires the matlab toolbox for first-order <http://www.columbia.edu/~mu2166/1st_order/1st_order.htm> and second-order <http://www.columbia.edu/~mu2166/2nd_order.htm> approximation of DSGE models produced by Schmitt-Grohe and Uribe (JEDC, 2004). 

(2) run the program calvo_model.m once for the peg economy and once for  the optimal exchange-rate economy. To this end, inside that program, set the variable optimal_policy to 0 or 1 to indicate the  peg economy or optimal exchange rate economy. This will produce the files calvo_opt.mat, calvo_opt_num_eval.m, calvo_peg.mat, and calvo_peg_num_eval.m. NOTE: calvo_model.m is different from the file with the same name in calvo.zip, because the present model includes welfare as an additional variable. 

(3) Run the program welfare_calvo_run. This script produces a table (called table)  with conditional and unconditional welfare costs as a function of the degree of price stickiness. 
