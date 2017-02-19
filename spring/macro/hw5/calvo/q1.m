load('simu_calvo_opt.mat')
[sigy0,sigx0]=mom(gx,hx,varshock);

stds = sqrt(diag(sigy0));

%serial correlations
[sigy1,sigx1]=mom(gx,hx,varshock,1);
scorr = diag(sigy1)./diag(sigy0);

e_rer = scorr(np);
P_1 = lagmatrix(P, 1);
d = P ./ P_1;
d = detrend(d(2:length(d)));
PAI = detrend(PAI);
PSSI = detrend(EPSI);
a = std(d);
b = std(PAI);
c = std(PSSI);


function [detrended] = detrend(dat)
dat = dat';
t = length(dat);
t1 = 1:t;
t2 = t1.^2;
y = dat';
X = [ones(1, t); t1; t2];
beta = inv(X * X') * (X * y);
proj = beta(1) + beta(2) * t1 + beta(3) * t2;
detrended = dat - proj;
end
