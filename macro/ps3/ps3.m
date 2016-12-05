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

z = [dat(:,vars('P'))'; dat(:,vars('Pc'))'; dat(:,vars('Y'))'; dat(:,vars('M1'))'; dat(:,vars('TR'))'; dat(:,vars('RFF'))'; dat(:,vars('NBR'))'];
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

% Long-run impact matrices
D=chol(Sigma)';

IRdiffRFF(:,5)  = D*[0;0;0;0;1;0;0];
IRlevelRFF(:,5) = D*[0;0;0;0;1;0;0];
for i = 6:19     
    IRdiffRFF(:,i) = A1*IRdiffRFF(:,i-1)+A2*IRdiffRFF(:,i-2)+A3*IRdiffRFF(:,i-3)+A4*IRdiffRFF(:,i-4);     
    IRlevelRFF(:,i)= IRlevelRFF(:,i-1)+IRdiffRFF(:,i);
end

b = figure();
subplot(4,2,2)
plot((IRlevelRFF(5,5:19)))
title('Output: RFF')
ylabel('Percent')
%xlabel('Quarters')
subplot(4,2,4)
plot(IRlevelRFF(7,5:19))
title('Pcom: RFF')
ylabel('Percent')
%xlabel('Quarters')
subplot(4,2,6)
plot(IRlevelRFF(4,5:19)*100)
xlabel('Quarters')
ylabel('Percent')
title('M1: RFF')
subplot(4,2,8)
plot((IRlevelRFF(1,5:19)))
title('Fed funds: RFF')
ylabel('Percent')
hold off

IRdiffNBR(:,5)  = D*[0;0;0;0;0;1;0];
IRlevelNBR(:,5) = D*[0;0;0;0;0;1;0];
for i = 6:19    
    IRdiffNBR(:,i) = A1*IRdiffNBR(:,i-1)+A2*IRdiffNBR(:,i-2)+A3*IRdiffNBR(:,i-3)+A4*IRdiffNBR(:,i-4);     
    IRlevelNBR(:,i)= IRlevelNBR(:,i-1)+IRdiffNBR(:,i);
end

c = figure();
subplot(4,2,2)
plot((IRlevelNBR(5,5:19)))
title('Output: NBR')
ylabel('Percent')
%xlabel('Quarters')
subplot(4,2,4)
plot(IRlevelNBR(7,5:19))
title('Pcom: NBR')
ylabel('Percent')
%xlabel('Quarters')
subplot(4,2,6)
plot(IRlevelNBR(4,5:19)*100)
xlabel('Quarters')
ylabel('Percent')
title('M1: NBR')
subplot(4,2,8)
plot((IRlevelNBR(1,5:19)))
title('Fed funds: NBR')
ylabel('Percent')
hold off