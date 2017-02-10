clear all

% Create the symbolic expressions for the model. 
%You need to run this only once. 
edeir_model %This program generates the file edeir_model_num_eval.m
OMEGA = 1.3975; %Frisch ela st. from Mendoza 1991
OLD_OMEGA = OMEGA;
DBAR = 0.7439; %debt
OLD_DBAR = DBAR;
PSSI = 0.0001; %debt elasticity of interest rate
OLD_PSSI = PSSI;
PHI = 0.0271; %capital adjustment cost
OLD_PHI = PHI;
RHO = 0.5300; %persistence of TFP shock
OLD_RHO = RHO;
ETATILDE = 0.0095; %standard deviation of innovation to TFP shock
OLD_ETATILDE = ETATILDE;
edeir_ss_4



%target sigma_y,rho_y,sigma_i,rho_i,sigma_h,tb/y
target = [2.81;0.62;9.82;0.31;2.02;0.02];
TARGET_TB_Y = 6;
TARGET_SIGMA_H = 5;
TARGET_SERIAL_Y = 2;
TARGET_SIGMA_I = 3;
TARGET_SERIAL_I = 4;
TARGET_SIGMA_Y = 1;
threshold = 0.1;
[stds, scorr, actual] = iterate(OMEGA,DBAR,PSSI,PHI,RHO,ETATILDE);
current_opt = compute_distance(actual, target);
prev = [];
best = [];
opt_vals = [];

%implementation of simulated annealing
iter = 1;
num_change = 1;
while current_opt > threshold
    current_opt
    OLD_OMEGA = OMEGA;
    OLD_DBAR = DBAR;
    OLD_PSSI = PSSI;
    OLD_PHI = PHI;
    OLD_RHO = RHO;
    OLD_ETATILDE = ETATILDE;
    if ~isempty(prev) && prev(TARGET_SIGMA_Y) > target(TARGET_SIGMA_Y)
        ETATILDE = perturb(ETATILDE, -1);
    elseif ~isempty(prev) ~= 0 && prev(TARGET_SIGMA_Y) < target(TARGET_SIGMA_Y)
        ETATILDE = perturb(ETATILDE, 1);
    else
        ETATILDE = perturb(ETATILDE, 0);
    end
    
    if ~isempty(prev) && prev(TARGET_SERIAL_Y) > target(TARGET_SERIAL_Y)
        RHO = perturb(RHO, -1);
    elseif ~isempty(prev) && prev(TARGET_SERIAL_Y) < target(TARGET_SERIAL_Y)
        RHO = perturb(RHO, 1);
    else
        RHO = perturb(RHO, 0);
    end
    
    if ~isempty(prev) && prev(TARGET_SERIAL_I) < target(TARGET_SERIAL_I)
        PSSI = perturb(PSSI, -1);
    elseif ~isempty(prev) && prev(TARGET_SERIAL_I) > target(TARGET_SERIAL_I)
        PSSI = perturb(PSSI, 1);
    else
        PSSI = perturb(PSSI, 0);
    end
    
    if ~isempty(prev) && prev(TARGET_SIGMA_H) > target(TARGET_SIGMA_H)
        OMEGA = perturb(OMEGA, 1);
    elseif ~isempty(prev) && prev(TARGET_SIGMA_H) < target(TARGET_SIGMA_H)
        OMEGA = perturb(OMEGA, -1);
    else
        OMEGA = perturb(OMEGA, 0);
    end
    
    if ~isempty(prev) && (prev(TARGET_SIGMA_I) < target(TARGET_SIGMA_I))
        PHI = perturb(PHI, 1);
    elseif ~isempty(prev) && (prev(TARGET_SIGMA_I) > target(TARGET_SIGMA_I))
        PHI = perturb(PHI, -1);
    else
        PHI = perturb(PHI,0);
    end
        
    y = ((1-ALFA)*KAPA^(ALFA * OMEGA))^(1/(OMEGA-1));
    DBAR = y*tby/RSTAR;
    [stds, scorr, actual] = iterate(OMEGA,DBAR,PSSI,PHI,RHO,ETATILDE);
    prev = actual;
    distance = compute_distance(actual, target);
    if distance > current_opt
        num_change = num_change + 1;
        OMEGA = OLD_OMEGA;
        DBAR = OLD_DBAR;
        PSSI = OLD_PSSI;
        PHI = OLD_PHI;
        RHO = OLD_RHO;
        ETATILDE = OLD_ETATILDE;
    else
       if actual(1) ~= 1000 && ~isnan(actual(1))
            current_opt = distance;
            best = actual;
            opt_vals = [OMEGA; DBAR; PSSI;PHI;RHO;ETATILDE];
       end
    end
    iter = iter + 1;
end

function new_value = perturb(old_value, sign)
    ran = 50;
    if sign == 1
        perturb_percent = (-.1 + rand(1)) / ran;
    elseif sign == -1
        perturb_percent = (-.9 + rand(1)) / ran;
    else
        perturb_percent = (-.5 + rand(1)) / ran;
    end
    new_value = old_value + perturb_percent * old_value;
end


function distance = compute_distance(actual, target)
    z = actual - target;
    distance = sqrt(z' * z);
end

function [stds, scorr, targets] = iterate(OMEGA1,DBAR1,PSSI1,PHI1,RHO1,ETATILDE1)
    OMEGA = OMEGA1;
    DBAR = DBAR1;
    PSSI = PSSI1;
    RHO = RHO1;
    ETATILDE = ETATILDE1;
    PHI=PHI1;
    edeir_ss_4
    edeir_model_num_eval;
    [gx,hx,exitflag]=gx_hx(nfy,nfx,nfyp,nfxp);
    varshock = nETASHOCK*nETASHOCK';

    %standard deviations
    try
        [sigy0,sigx0]=mom(gx,hx,varshock);
        stds = sqrt(diag(sigy0));
    
    %serial correlations
        [sigy1,sigx1]=mom(gx,hx,varshock,1);
        scorr = diag(sigy1)./diag(sigy0);
        targets = [stds(noutput)*100; scorr(noutput); stds(nivv)*100; scorr(nivv);stds(nh)*100;tby];
    catch
        infinity = 1000;
        stds = [infinity;infinity;infinity;infinity;infinity;infinity];
        scorr = stds;
        targets = stds;
    end
 end
