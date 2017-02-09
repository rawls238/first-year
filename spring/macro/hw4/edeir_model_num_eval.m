%edeir_model_num_eval.m was created by running  /Users/garidor/Desktop/first-year/spring/macro/hw4/edeir_model.m
nfx = [
[ r + 1, d,                                         0,                        0]
[     0, 0,        ALFA*a*h^(1 - ALFA)*k*k^(ALFA - 1),    a*h^(1 - ALFA)*k^ALFA]
[     0, 0,                                         0,                        0]
[     0, 0, -(ALFA*a*k*(k/h)^(ALFA - 1)*(ALFA - 1))/h, -a*(k/h)^ALFA*(ALFA - 1)]
[     0, 0,                                         0,                        0]
[     0, 0,                                  PHI*k*la,                        0]
[     0, 0,                                         0,                        0]
[     0, 0,                             k*(DELTA - 1),                        0]
[     0, 0,                                         0,                        0]
[     0, 0,                                         0,                        0]
[     1, 0,                                         0,                        0]
[     0, 0,                                         0,                        0]
[     0, 0,                                         0,                      RHO]
[     0, 0,                                         0,                        0]
];
nfxp = [
[                 -1,        0,                                                                        0,                                0]
[                  0,        0,                                                                        0,                                0]
[                  0,        0,                                                                        0,                                0]
[                  0,        0,                                                                        0,                                0]
[                  0, BETTA*la,                                                                        0,                                0]
[                  0,        0, - PHI*k*la - BETTA*la*(PHI*k - (ALFA*a*k*(k/h)^(ALFA - 2)*(ALFA - 1))/h), ALFA*BETTA*a*la*(k/h)^(ALFA - 1)]
[ PSSI*exp(d - DBAR),       -1,                                                                        0,                                0]
[                  0,        0,                                                                        k,                                0]
[                  0,        0,                                                                        0,                                0]
[                  0,        0,                                                                        0,                                0]
[                 -1,        0,                                                                        0,                                0]
[                  0,        0,                                                                        0,                                0]
[                  0,        0,                                                                        0,                               -1]
[                  0,        0,                                                                        k,                                0]
];
nfy = [
[                                        c,  ivv,    -output,                                                                      0,   0,    0,        0,  0,        0,  0]
[                                        0,    0,    -output,                                        -(a*h*k^ALFA*(ALFA - 1))/h^ALFA,   0,    0,        0,  0,        0,  0]
[ -(SIGG*c)/(c - h^OMEGA/OMEGA)^(SIGG + 1),    0,          0,                  (SIGG*h*h^(OMEGA - 1))/(c - h^OMEGA/OMEGA)^(SIGG + 1), -la,    0,        0,  0,        0,  0]
[                                        0,    0,          0, (ALFA*a*k*(k/h)^(ALFA - 1)*(ALFA - 1))/h - h*h^(OMEGA - 2)*(OMEGA - 1),   0,    0,        0,  0,        0,  0]
[                                        0,    0,          0,                                                                      0, -la,    0,        0,  0,        0,  0]
[                                        0,    0,          0,                                                                      0, -la,    0,        0,  0,        0,  0]
[                                        0,    0,          0,                                                                      0,   0,    0,        0,  0,        0,  0]
[                                        0, -ivv,          0,                                                                      0,   0,    0,        0,  0,        0,  0]
[                                       -c, -ivv,     output,                                                                      0,   0,    0,       -1,  0,        0,  0]
[                                        0,    0, -tb/output,                                                                      0,   0,    0, 1/output, -1,        0,  0]
[                                        0,    0,          0,                                                                      0,   0,    0,        0,  0,       -1,  0]
[                                        0,    0, -ca/output,                                                                      0,   0,    0,        0,  0, 1/output, -1]
[                                        0,    0,          0,                                                                      0,   0,    0,        0,  0,        0,  0]
[                                        0,    0,          0,                                                                      0,   0, -kfu,        0,  0,        0,  0]
];
nfyp = [
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                BETTA*la*(r + 1),                0, 0, 0, 0, 0]
[ 0, 0, 0, -(ALFA*BETTA*a*k*la*(k/h)^(ALFA - 2)*(ALFA - 1))/h, -BETTA*la*(DELTA + PHI*(k - kfu) - ALFA*a*(k/h)^(ALFA - 1) - 1), BETTA*PHI*kfu*la, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
[ 0, 0, 0,                                                  0,                                                               0,                0, 0, 0, 0, 0]
];
nf = [
[ c - d + ivv - output + d*(r + 1), a*h^(1 - ALFA)*k^ALFA - output, 1/(c - h^OMEGA/OMEGA)^SIGG - la, - h^(OMEGA - 1) - a*(k/h)^ALFA*(ALFA - 1), BETTA*la*(r + 1) - la, - la - BETTA*la*(DELTA + PHI*(k - kfu) - ALFA*a*(k/h)^(ALFA - 1) - 1), RSTAR - r + PSSI*(exp(d - DBAR) - 1), k - ivv + k*(DELTA - 1), output - ivv - c - tb, tb/output - tby, -ca, ca/output - cay, RHO*log(a) - log(a), k - kfu]
];
nf=transpose(nf);
nETASHOCK = [
        0
        0
        0
 ETATILDE
];
nvarshock = [
[ 0, 0, 0,                       0]
[ 0, 0, 0,                       0]
[ 0, 0, 0,                       0]
[ 0, 0, 0, ETATILDE*conj(ETATILDE)]
];
nvarme = [
];
nd = 1;
nr = 2;
nk = 3;
na = 4;
nstate = 4;
nc = 1;
nivv = 2;
noutput = 3;
nh = 4;
nla = 5;
nkfu = 6;
ntb = 7;
ntby = 8;
nca = 9;
ncay = 10;
ncontrol = 10;
