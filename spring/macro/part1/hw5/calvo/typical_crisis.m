%typical_crisis.m
%Simulate the dynamics of a typical external crisis  in a two-sector model with sticky prices a la Calvo-Yun   in the nontraded sector as  developed in the  chapter entitled  ``Nominal Rigidity, Exchange Rates,  And Unemployment''
%of the book ``Open Economy Macroeconomics,'' by Martín Uribe and Stephanie Schmitt-Grohé, Princeton University Press. 
%
%© Martín Uribe and Stephanie Schmitt-Grohé, January 29, 2016 (BSGDAY)

clear all
clf

orient tall

load simu_calvo_opt 
%produced by running simu_calvo.m in
%c:\data\uribe\book\dnwr\calvo
tipe = '--';


my = mean(YT);
stdy = std(YT);

X = lagg(YT,10);
s1 = find(X(:,11)>=my & X(:,1)<my-2*stdy);

s1 = s1 + 10; 

ww1= 10;
ww2 = 30; 

find(s1<=ww1|s1>=length(YT)-ww2);
s1 = setdiff(s1,s1(ans));
nw = numel(s1);

ww = (ww1+ww2)/2; 

WYT = zeros(nw,2*ww+1);
WRS = zeros(nw,2*ww+1);
WH = zeros(nw,2*ww+1);
WW = zeros(nw,2*ww+1);
WEPSI = zeros(nw,2*ww+1);
WP = zeros(nw,2*ww+1);
WPAI = zeros(nw,2*ww+1);
WCT = zeros(nw,2*ww+1);
WTBY = zeros(nw,2*ww+1);
WDY = zeros(nw,2*ww+1);

ww1=ww;
ww2=ww;
for i=1:nw
WYT(i,:) = YT(s1(i)-ww1:s1(i)+ww2)';
WRS(i,:) = RS(s1(i)-ww1:s1(i)+ww2)';
WCT(i,:) = CT(s1(i)-ww1:s1(i)+ww2)';

WH(i,:) = H(s1(i)-ww1:s1(i)+ww2)';
WW(i,:) = W(s1(i)-ww1:s1(i)+ww2)';
WEPSI(i,:) = EPSI(s1(i)-ww1:s1(i)+ww2)';
WP(i,:) = P(s1(i)-ww1:s1(i)+ww2)';
WPAI(i,:) = PAI(s1(i)-ww1:s1(i)+ww2)';
WTBY(i,:) = TBY(s1(i)-ww1:s1(i)+ww2)';
WDY(i,:) = DY(s1(i)-ww1:s1(i)+ww2)';
end
 
thick = 2;
rows = 4;
cols=2; 

ww1=10;
ww2=30;

%%%%
var_name = { 'Hours', 'Real Wage (in terms of tradables)', 'Annualized Devaluation Rate', 'Relative Price of Nontradables $(P^N_t/\mathcal{E}_t)$', 'Annual CPI Inflation Rate', 'Consumption of Tradables', 'Trade-Balance-To-Output Ratio', 'Debt-To-Output Ratio (Annual)'};
ylabels = { '\% dev. from trend', '\% dev. from trend',    '\% ', '\% dev. from trend',    '\% ', '\% dev. from trend', '\% ', '\%'};

W = zeros(nw, 41, 8); 
W(:,:,1) = WH; 
W(:,:,2) = WW;
W(:,:,3) = WEPSI; 
W(:,:,4) = WP;
W(:,:,5) = WPAI; 
W(:,:,6) = WCT;
W(:,:,7) = WTBY; 
W(:,:,8) = WDY;

for j=1:(rows*cols)
    
subplot(rows,cols,j)
x =  mean(W(:, :,j));
plot(-ww1:ww2,x,'linestyle', tipe, 'linewidth',thick);
title(var_name{j},'interpreter','LaTeX')
ylabel(ylabels{j},'interpreter','LaTeX')
xlabel('$t$', 'interpreter','LaTeX')
xlim([-ww1 ww2])
hold on 
end



%%%%%%%%%%%%%%
clear all
load simu_calvo_peg 
%produced by running simu_calvo.m in
%c:\data\uribe\book\dnwr\calvo
%you must modify the load command
tipe = '-'; 

my = mean(YT);
stdy = std(YT);

X = lagg(YT,10);
s1 = find(X(:,11)>=my & X(:,1)<my-2*stdy);

s1 = s1 + 10; 

ww1= 10;
ww2 = 30; 

find(s1<=ww1|s1>=length(YT)-ww2);
s1 = setdiff(s1,s1(ans));
nw = numel(s1);

ww = (ww1+ww2)/2; 

WYT = zeros(nw,2*ww+1);
WRS = zeros(nw,2*ww+1);
WH = zeros(nw,2*ww+1);
WW = zeros(nw,2*ww+1);
WEPSI = zeros(nw,2*ww+1);
WP = zeros(nw,2*ww+1);
WPAI = zeros(nw,2*ww+1);
WCT = zeros(nw,2*ww+1);
WTBY = zeros(nw,2*ww+1);
WDY = zeros(nw,2*ww+1);

ww1=ww;
ww2=ww;
for i=1:nw
WYT(i,:) = YT(s1(i)-ww1:s1(i)+ww2)';
WRS(i,:) = RS(s1(i)-ww1:s1(i)+ww2)';
WCT(i,:) = CT(s1(i)-ww1:s1(i)+ww2)';

WH(i,:) = H(s1(i)-ww1:s1(i)+ww2)';
WW(i,:) = W(s1(i)-ww1:s1(i)+ww2)';
WEPSI(i,:) = EPSI(s1(i)-ww1:s1(i)+ww2)';
WP(i,:) = P(s1(i)-ww1:s1(i)+ww2)';
WPAI(i,:) = PAI(s1(i)-ww1:s1(i)+ww2)';
WTBY(i,:) = TBY(s1(i)-ww1:s1(i)+ww2)';
WDY(i,:) = DY(s1(i)-ww1:s1(i)+ww2)';
end


thick = 2;
rows = 4;
cols=2;

ww1=10;
ww2=30;

%%%%
var_name = { 'Hours', 'Real Wage (in terms of tradables)', 'Annualized Devaluation Rate', 'Relative Price of Nontradables $(P^N_t/\mathcal{E}_t)$', 'Annual CPI Inflation Rate', 'Consumption of Tradables', 'Trade-Balance-To-Output Ratio', 'Debt-To-Output Ratio (Annual)'};
ylabels = { '\% dev. from trend', '\% dev. from trend',    '\% ', '\% dev. from trend',    '\% ', '\% dev. from trend', '\% ', '\%'};

W = zeros(nw, 41, 8);
W(:,:,1) = WH;
W(:,:,2) = WW;
W(:,:,3) = WEPSI; 
W(:,:,4) = WP;
W(:,:,5) = WPAI;
W(:,:,6) = WCT;
W(:,:,7) = WTBY;
W(:,:,8) = WDY;




for j=1:(rows*cols)
    
subplot(rows,cols,j)
x =  mean(W(:, :,j));
plot(-ww1:ww2,x,'linestyle', tipe, 'linewidth',thick);
title(var_name{j},'interpreter','LaTeX')
ylabel(ylabels{j},'interpreter','LaTeX')
xlim([-ww1 ww2])
hold off  
end

