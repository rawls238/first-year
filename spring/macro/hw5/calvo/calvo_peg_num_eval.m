% File name: calvo_peg_num_eval.m 
% File generated by anal_deriv_print2f.m Date: 05-Feb-2016

nfx=[[-1] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [1/(4*output)] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [paiN] [0] [0] [0] [0] [0] [0] [0] [0] [THETA*paiN^(MU/ALFA)*s] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [(p*paiN)/epsi] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [yT] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [yT] [yT/output] [0] [0] [a11] [a21] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [1] [0] [0] [0] [0] [0] [0] [a12/(rstar + 1)] [a22/(rstar + 1)] [0]];
nfxp=[[1/(r + 1)] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [PSSI*exp(d - DBAR)] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-a2] [-paiN] [0] [0] [0] [0] [0] [0] [0] [-(ALFA*h*(h/s)^(ALFA - 1))/s] [-s] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-p] [0] [0] [0] [0] [-p] [0] [(p*pNtilde^(1 - MU)*yN*(MU - 1))/MU] [0] [0] [0] [0] [p*yN] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-1] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-1/(rstar + 1)] [0]];
nfy=[[-cT] [(A^2*cT*(1/XI - 1)*(1/(1/XI - 1) + 1))/(cT^(2/XI)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1) + 2)*(1/(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1)))^SIGG) - (A*cT)/(XI*cT^(1/XI + 1)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1) + 1)*(1/(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1)))^SIGG) - (A^2*SIGG*cT)/(cT^(2/XI)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(2/(1/XI - 1) + 2)*(1/(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1)))^(SIGG + 1))] [-(cT*cT^(1/XI - 1)*(A - 1))/(A*XI*yN^(1/XI))] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-cT/output] [-(A*cT*(1/XI - 1)*(A - 1)*(1/(1/XI - 1) + 1))/(cT^(1/XI)*yN^(1/XI)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1) + 2))] [0] [0] [0] [0] [0] [0] [0] [-(VARPHI*h)/(HBAR - h)^2] [(ALFA*h*(h/s)^(ALFA - 1))/s] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-la] [0] [la*w] [0] [0] [0] [0] [0] [-BETTA*THETA*pvmr*(1/paiN)^(1 - MU)] [-(BETTA*THETA*pvmc)/(1/paiN)^(MU/ALFA)] [0] [-la/(r + 1)] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [la*w] [0] [0] [0] [0] [0] [0] [-(w*yN^(1/ALFA)*(TAU - 1))/(ALFA*pNtilde^(MU/ALFA))] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [(MU*pNtilde*(THETA - 1))/(ALFA*pNtilde^(MU/ALFA + 1))] [(pNtilde*(MU - 1)*(THETA - 1))/pNtilde^MU] [0] [0] [- (BETTA*THETA*pvmr*(MU - 1))/(paiN*(1/paiN)^MU) - (p*pNtilde*yN*(MU - 1)^2)/(MU*pNtilde^MU)] [(MU*pNtilde*w*yN^(1/ALFA)*(TAU - 1))/(ALFA^2*pNtilde^(MU/ALFA + 1)) - (BETTA*MU*THETA*pvmc)/(ALFA*paiN*(1/paiN)^(MU/ALFA + 1))] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [(A*SIGG*yN*(A - 1))/(cT^(1/XI)*yN^(1/XI)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(2/(1/XI - 1) + 2)*(1/(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1)))^(SIGG + 1)) - (A*yN*(1/XI - 1)*(A - 1)*(1/(1/XI - 1) + 1))/(cT^(1/XI)*yN^(1/XI)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1) + 2)*(1/(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1)))^SIGG)] [(cT^(1/XI)*yN*(A - 1))/(A*XI*yN^(1/XI + 1))] [0] [-yN] [0] [0] [0] [0] [(p*pNtilde^(1 - MU)*yN*(MU - 1))/MU] [-(w*yN*yN^(1/ALFA - 1)*(TAU - 1))/(ALFA^2*pNtilde^(MU/ALFA))] [0] [0] [0] [p*yN] [0] [(yN*(A - 1))/(XI*yN^(1/XI + 1)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1) + 1)) + (yN*(1/XI - 1)*(A - 1)^2*(1/(1/XI - 1) + 1))/(yN^(2/XI)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1) + 2))] [0] [0] [0] [0] [0] [0] [0] [0] [0] [(MU*THETA*paiN*paiN^(MU/ALFA - 1)*s)/ALFA] [THETA*paiN*paiN^(MU - 2)*(MU - 1)] [(p*paiN)/epsi] [0] [0] [0] [0] [0] [0] [0] [0] [0] [paiN] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-pvmc] [0] [-pvmc] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [pvmr] [-pvmr] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-(p*paiN)/epsi] [0] [0] [0] [0] [0] [epsi] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-pai] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-output] [(cT - yT)/output] [0] [0] [0] [0] [-d/(4*output)] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-1] [0] [0] [0] [0] [0] [-d/(r + 1)^2] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-1] [la/(r + 1)^2] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [-1]];
nfyp=[[0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [BETTA*THETA*pvmr*(1/paiN)^(1 - MU)] [(BETTA*THETA*pvmc)/(1/paiN)^(MU/ALFA)] [0] [BETTA*la] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [(BETTA*THETA*pvmr*(MU - 1))/(paiN*(1/paiN)^MU)] [(BETTA*MU*THETA*pvmc)/(ALFA*paiN*(1/paiN)^(MU/ALFA + 1))] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [(BETTA*THETA*pvmr*(MU - 1))/(paiN*(1/paiN)^MU)] [(BETTA*MU*THETA*pvmc)/(ALFA*paiN*(1/paiN)^(MU/ALFA + 1))] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [(BETTA*THETA*pvmc)/(1/paiN)^(MU/ALFA)] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [BETTA*THETA*pvmr*(1/paiN)^(1 - MU)] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0] [0]];
nf=[[yT - d - cT + d/(r + 1)] [A/(cT^(1/XI)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1) + 1)*(1/(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1)))^SIGG) - la] [- p - (cT^(1/XI)*(A - 1))/(A*yN^(1/XI))] [la*w - VARPHI/(HBAR - h)] [(h/s)^ALFA - yN] [THETA*paiN^(MU/ALFA)*s - (THETA - 1)/pNtilde^(MU/ALFA) - s] [THETA*paiN^(MU - 1) - pNtilde^(1 - MU)*(THETA - 1) - 1] [(p*paiN)/epsi - p] [pvmr - pvmc] [BETTA*THETA*pvmr*(1/paiN)^(1 - MU) - pvmr + (p*pNtilde^(1 - MU)*yN*(MU - 1))/MU] [(BETTA*THETA*pvmc)/(1/paiN)^(MU/ALFA) - pvmc - (w*yN^(1/ALFA)*(TAU - 1))/(ALFA*pNtilde^(MU/ALFA))] [rstar - r + PSSI*(exp(d - DBAR) - 1)] [BETTA*la - la/(r + 1)] [epsi - 1] [yT - output + p*yN] [- tby - (cT - yT)/output] [- a2 - (A - 1)/(yN^(1/XI)*(A*cT^(1 - 1/XI) - yN^(1 - 1/XI)*(A - 1))^(1/(1/XI - 1) + 1))] [paiN - pai] [a12*log((rstar + 1)/(RSTAR + 1)) - log(yT) + a11*log(yT)] [a22*log((rstar + 1)/(RSTAR + 1)) - log((rstar + 1)/(RSTAR + 1)) + a21*log(yT)] [d/(4*output) - dy]];
nETASHOCK=[[0] [0] [0] [0] [eta11] [eta21] [0] [0] [0] [0] [eta12] [eta22]];
nfx=reshape(nfx,[21   6]);
nfxp=reshape(nfxp,[21   6]);
nfy=reshape(nfy,[21  15]);
nfyp=reshape(nfyp,[21  15]);
nf=reshape(nf,[21   1]);

nETASHOCK=reshape(nETASHOCK,[6  2]);
