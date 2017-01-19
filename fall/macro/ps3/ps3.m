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

%2.b
P = 1; Pc = 2; Y = 3; RFF = 4; TR = 5; NBR = 6; M1 = 7;
z = [dat(:,vars('P'))'; dat(:,vars('Pc'))'; dat(:, vars('Y'))'; dat(:, vars('RFF'))'; dat(:, vars('TR'))'; dat(:, vars('NBR'))'; dat(:, vars('M1'))'];
Z = [(1:146)' ones(146, 1) lagmatrix(z',1) lagmatrix(z',2) lagmatrix(z',3) lagmatrix(z',4)];

z = z(:,VARp+1:length(z));

Z = Z(VARp+1:length(Z),:)';

% GLS estimator
beta  = kron(inv(Z*Z')*Z,eye(VARn))*z(:);
Sigma = (T-VARn*VARp-1)^(-1)*z*(eye(length(z))-Z'*inv((Z*Z'))*Z)*z';

A1 = reshape(beta(15:63),7,7);
A2 = reshape(beta(64:112),7,7);
A3 = reshape(beta(113:161),7,7);
A4 = reshape(beta(162:210),7,7);

D=chol(Sigma)';
IRdiffRFF(:,5) = D*[0;0;0;1;0;0;0];
for i = 6:19
    IRdiffRFF(:,i) = A1*IRdiffRFF(:,i-1)+A2*IRdiffRFF(:,i-2)+A3*IRdiffRFF(:,i-3)+A4*IRdiffRFF(:,i-4);     
end

plotGr(IRdiffRFF, Y, Pc, M1, TR, NBR, RFF, P);

% 2.c
P = 1; Pc=2; Y = 3; NBR = 4; TR = 5; RFF = 6; M1 = 7;
z = [dat(:,vars('P'))'; dat(:,vars('Pc'))'; dat(:, vars('Y'))'; dat(:, vars('NBR'))'; dat(:, vars('TR'))'; dat(:, vars('RFF'))'; dat(:, vars('M1'))'];

Z = [(1:T)' ones(146,1) lagmatrix(z',1) lagmatrix(z',2) lagmatrix(z',3) lagmatrix(z',4)];

z = z(:,VARp+1:length(z));

Z = Z(VARp+1:length(Z),:)';

% GLS estimator
beta  = kron(inv(Z*Z')*Z,eye(VARn))*z(:);
Sigma = (T-VARn*VARp-1)^(-1)*z*(eye(length(z))-Z'*inv((Z*Z'))*Z)*z';
D = chol(Sigma)';

A1 = reshape(beta(15:63),7,7);
A2 = reshape(beta(64:112),7,7);
A3 = reshape(beta(113:161),7,7);
A4 = reshape(beta(162:210),7,7);

IRdiffNBR(:,5)  = D*[0;0;0;-1;0;0;0];
for i = 6:19
    IRdiffNBR(:,i) = A1*IRdiffNBR(:,i-1)+A2*IRdiffNBR(:,i-2)+A3*IRdiffNBR(:,i-3)+A4*IRdiffNBR(:,i-4);     
end

plotGr(IRdiffNBR, Y, Pc, M1, TR, NBR, RFF, P);


% 2.d
P = 1; Y = 2; RFF = 3; NBR = 4; TR = 5; M1 = 6;
z = [dat(:,vars('P'))'; dat(:, vars('Y'))'; dat(:, vars('RFF'))'; dat(:, vars('TR'))'; dat(:, vars('NBR'))'; dat(:, vars('M1'))'];

Z = [(1:T)' ones(146,1) lagmatrix(z',1) lagmatrix(z',2) lagmatrix(z',3) lagmatrix(z',4);];

z = z(:,VARp+1:length(z));

Z = Z(VARp+1:length(Z),:)';

VARn = 6;

% GLS estimator
beta  = kron(inv(Z*Z')*Z,eye(VARn))*z(:);
Sigma = (T-VARn*VARp-1)^(-1)*z*(eye(length(z))-Z'*inv((Z*Z'))*Z)*z';
D = chol(Sigma)';

A1 = reshape(beta(13:48),6,6);
A2 = reshape(beta(49:84),6,6);
A3 = reshape(beta(85:120),6,6);
A4 = reshape(beta(121:156),6,6);

IRdiffnoPc(:,5)  = D*[0;0;1;0;0;0];
for i = 6:19
    IRdiffnoPc(:,i) = A1*IRdiffnoPc(:,i-1)+A2*IRdiffnoPc(:,i-2)+A3*IRdiffnoPc(:,i-3)+A4*IRdiffnoPc(:,i-4);     
end

c = figure();
subplot(1,2,2)
plot(IRdiffnoPc(P,5:19))
title('Price: No PC')
ylabel('Percent')
xlabel('Quarters')
hold off

function c = plotGr(IR, Y, Pc, M1, TR, NBR, RFF, P)
    c = figure();
    subplot(4,2,1)
    plot(IR(Y,5:19))
    title('Output')
    ylabel('Percent')
    xlabel('Quarters')
    subplot(4,2,2)
    plot(IR(Pc,5:19))
    title('Pcom')
    ylabel('Percent')
    xlabel('Quarters')
    subplot(4,2,3)
    plot(IR(M1,5:19))
    xlabel('Quarters')
    ylabel('Percent')
    title('M1')
    subplot(4,2,4)
    plot(IR(TR,5:19))
    title('Total Reserves=')
    ylabel('Percent')
    subplot(4,2,5)
    plot(IR(NBR,5:19))
    title('Non-borrowed Reserves')
    ylabel('Percent')
    subplot(4,2,6)
    plot(IR(RFF,5:19))
    title('RFF')
    ylabel('Percent')
    subplot(4,2,7)
    plot(IR(P,5:19))
    title('P')
    ylabel('Percent')
    hold off
end

