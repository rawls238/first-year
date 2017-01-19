%%  Identifying Monetary Shocks
% Xiao Xu, xx2251, Macro PS3, Prof. Mertens
clear all;
clc;
filename = 'CEE1999data.xls';
data = xlsread(filename);

dates = data(:,end,1);
T = length(dates);

y = data(:,2);  % Output
p = data(:,3); % implicit price deflator
pc = data(:,4); % change in CPI
rff = data(:,5);  % effective federal funds rate
tr = data(:,7); % total reserve
nbr = data(:,6); % nonborrowed reserves
m1 = data(:,8); % m1
m2 = data(:,9); %m2
mb = data(:,10); %mb

VARn = 7;
VARp = 4;

% construct small z
z = [ y'; p'; pc';rff'; tr'; nbr'; m1'];
Z = [lagmatrix(z',1) lagmatrix(z',2) lagmatrix(z',3) lagmatrix(z',4)];

z = z(:,VARp+1:length(z));

Z = Z(VARp+1:length(Z),:)';
Z = [ones(1,length(Z));Z];

% GLS estimator
beta  = kron(inv(Z*Z')*Z,eye(VARn))*z(:);
Sigma = (T-VARn*VARp-1)^(-1)*z*(eye(length(z))-Z'*inv((Z*Z'))*Z)*z';

A0 = beta(1:2)
A1 = reshape(beta(3:6),2,2)
A2 = reshape(beta(7:10),2,2)
A3 = reshape(beta(11:14),2,2)
A4 = reshape(beta(15:18),2,2)

% Federal fund rate model
% D matrix with recursive assumption
D    = chol(Sigma, 'lower');           % Estimated contemporaneous impact matrix:

%Compute Impulse Response to RFF shock
IRdiffT(:,5)  = D*[0;0;0;1;0;0;0];
IRlevelT(:,5) = D*[0;0;0;1;0;0;0];
for i = 6:18     
IRdiffT(:,i) = A1*IRdiffT(:,i-1)+A2*IRdiffT(:,i-2)+A3*IRdiffT(:,i-3)+A4*IRdiffT(:,i-4);     
IRlevelT(:,i)= IRlevelT(:,i-1)+IRdiffT(:,i);
end


% Plot Results
h1 = figure
subplot(3,2,1)
plot((IRlevelT(1,:))*100)
title('Productivity: Technology Shock')
ylabel('Percent')
%xlabel('Quarters')
subplot(3,2,3)
plot((IRlevelT(1,:)+IRlevelT(2,:))*100)
title('Output: Technology Shock')
ylabel('Percent')
%xlabel('Quarters')
subplot(3,2,5)
plot(IRlevelT(2,:)*100)
xlabel('Quarters')
ylabel('Percent')
title('Hours: Technology Shock')

subplot(3,2,2)
plot((IRlevelNT(1,:))*100)
title('Productivity: Non-Technology Shock')
ylabel('Percent')
%xlabel('Quarters')
subplot(3,2,4)
plot((IRlevelNT(1,:)+IRlevelNT(2,:))*100)
title('Output: Non-Technology Shock')
ylabel('Percent')
%xlabel('Quarters')
subplot(3,2,6)
plot(IRlevelNT(2,:)*100)
xlabel('Quarters')
ylabel('Percent')
title('Hours: Non-Technology Shock')
hold off


