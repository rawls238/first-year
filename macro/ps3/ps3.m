[dat, headers] = xlsread('CEE1999data.xls');
vals = zeros(1, length(headers));
for i=1:length(vals)
    vals(i) = i;
end
vars = containers.Map(headers, vals);

T = length(dat(:,vars('ENTRY')));
VARn = 7;
VARp = 4;

z = [dat(:,vars('RFF'))'; dat(:,vars('TR'))'; dat(:,vars('NBR'))'; dat(:,vars('M1'))'; dat(:,vars('Y'))'; dat(:,vars('P'))'; dat(:,vars('Pc'))'];
Z = [lagmatrix(z',1) lagmatrix(z',2) lagmatrix(z',3) lagmatrix(z',4)];

z = z(:,VARp+1:length(z));

Z = Z(VARp+1:length(Z),:)';
Z = [ones(1,length(Z));Z];

% GLS estimator
beta  = kron(inv(Z*Z')*Z,eye(VARn))*z(:);
Sigma = (T-VARn*VARp-1)^(-1)*z*(eye(length(z))-Z'*inv((Z*Z'))*Z)*z';

A0 = beta(1:7);
A1 = reshape(beta(8:56),7,7);
A2 = reshape(beta(57:105),7,7);
A3 = reshape(beta(106:154),7,7);
A4 = reshape(beta(155:203),7,7);

% Long-run impact matrices
Phi  = inv(eye(VARn)-A1-A2-A3-A4);
Phie = (chol(Phi*Sigma*Phi'))'; % Estimated identified long run impact matrix:
D    = inv(Phi)*Phie;           % Estimated contemporaneous impact matrix: