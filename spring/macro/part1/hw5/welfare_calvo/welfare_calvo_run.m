%welfare_calvo_run.m
%Compute the welfare cost of a currency peg vis-a-vis the optimal exchange-rate 
%policy in the model  with staggered price setting a la Calvo  in the nontraded sector   developed in the chapter entitled  ``Nominal Rigidity, Exchange Rates,  And Unemployment''
%of the book ``Open Economy Macroeconomics,'' by Martín Uribe and Stephanie Schmitt-Grohé, Princeton University Press. 
%Output: table (8x4 matrix)
%column1: Average duration of prices, in quarters.
%column 2: THETA, probability of not being allowed to change the price in a given quarter.
%column 3: Unconditional welfare cost of a currency peg, as percentage of consumption each period
%column 4: Conditional welfare cost of a currency peg, as percentage of consumption each period
%
%© Martín Uribe and Stephanie Schmitt-Grohé, January 29, 2016 . 
table = zeros(8,4)
jjj = 1;
for duration = [1 2 1/0.3 4 5 6 7  8]

THETA = 1-1/duration;


welfare_calvo

table(jjj,1:4) = [duration 1-1/duration uwelfare_cost cwelfare_cost];
jjj = jjj+1;
end
clc
table