clear all;
close all;

[dat, headers] = xlsread('CEE1999data.xls');
vals = zeros(1, length(headers));
for i=1:length(vals)
    vals(i) = i;
end
vars = containers.Map(headers, vals);

T = length(dat(:,vars('ENTRY')));
VARn = 7;
VARp = 4;

P = 1; Pc = 2; Y = 3; RFF = 4; TR = 5; NBR = 6; M1 = 7;
z = [dat(:,vars('P'))'; dat(:,vars('Pc'))'; dat(:, vars('Y'))'; dat(:, vars('RFF'))'; dat(:, vars('TR'))'; dat(:, vars('NBR'))'; dat(:, vars('M1'))'];
Z = [lagmatrix(z',1) lagmatrix(z',2) lagmatrix(z',3) lagmatrix(z',4);];

z = z(:,VARp+1:length(z));

Z = Z(VARp+1:length(Z),:)';

% GLS estimator
beta  = kron(inv(Z*Z')*Z,eye(VARn))*z(:);
Sigma = (T-VARn*VARp-1)^(-1)*z*(eye(length(z))-Z'*inv((Z*Z'))*Z)*z';

A1 = reshape(beta(1:49),7,7);
A2 = reshape(beta(50:98),7,7);
A3 = reshape(beta(99:147),7,7);
A4 = reshape(beta(148:196),7,7);

D=chol(Sigma)';

IRdiffRFF(:,5)  = D*[0;0;0;1;0;0;0];
for i = 6:19     
    IRdiffRFF(:,i) = A1*IRdiffRFF(:,i-1)+A2*IRdiffRFF(:,i-2)+A3*IRdiffRFF(:,i-3)+A4*IRdiffRFF(:,i-4);     
end

b = figure();
subplot(4,2,2)
plot((IRdiffRFF(Y,5:19)))
title('Output: RFF')
ylabel('Percent')
xlabel('Quarters')
subplot(4,2,4)
plot(IRdiffRFF(Pc,5:19))
title('Pcom: RFF')
ylabel('Percent')
xlabel('Quarters')
subplot(4,2,6)
plot(IRdiffRFF(M1,5:19))
xlabel('Quarters')
ylabel('Percent')
title('M1: RFF')
subplot(4,2,8)
plot((IRdiffRFF(RFF,5:19)))
title('Fed funds: RFF')
ylabel('Percent')
hold off

P = 1; Pc = 2; Y = 3; NBR = 4; TR = 5; RFF = 6; M1 = 7;
z = [dat(:,vars('P'))'; dat(:,vars('Pc'))'; dat(:, vars('Y'))'; dat(:, vars('NBR'))'; dat(:, vars('TR'))'; dat(:, vars('RFF'))'; dat(:, vars('M1'))'];

Z = [lagmatrix(z',1) lagmatrix(z',2) lagmatrix(z',3) lagmatrix(z',4);];

z = z(:,VARp+1:length(z));

Z = Z(VARp+1:length(Z),:)';

% GLS estimator
beta  = kron(inv(Z*Z')*Z,eye(VARn))*z(:);
Sigma = (T-VARn*VARp-1)^(-1)*z*(eye(length(z))-Z'*inv((Z*Z'))*Z)*z';
D = chol(Sigma)';

IRdiffNBR(:,5)  = D*[0;0;0;1;0;0;0];
for i = 6:19    
    IRdiffNBR(:,i) = A1*IRdiffNBR(:,i-1)+A2*IRdiffNBR(:,i-2)+A3*IRdiffNBR(:,i-3)+A4*IRdiffNBR(:,i-4);     
end

c = figure();
subplot(4,2,2)
plot((IRdiffNBR(Y,5:19)))
title('Output: NBR')
ylabel('Percent')
xlabel('Quarters')
subplot(4,2,4)
plot(IRdiffNBR(Pc,5:19))
title('Pcom: NBR')
ylabel('Percent')
xlabel('Quarters')
subplot(4,2,6)
plot(IRdiffNBR(M1,5:19))
xlabel('Quarters')
ylabel('Percent')
title('M1: NBR')
subplot(4,2,8)
plot((IRdiffNBR(RFF,5:19)))
title('Fed funds: NBR')
ylabel('Percent')
hold off