load('simu_calvo_peg.mat')

std(PAI)
std(EPSIRER)
std(EPSI)
ep = autocorr(EPSI, 1);
re = autocorr(EPSIRER, 1);
ep(2)
re(2)