clear all

% Create the symbolic expressions for the model. 
%You need to run this only once. 
edeir_model %This program generates the file edeir_model_num_eval.m

PERTURBATION_PERCENT = 0.01;

OMEGA = 1.455; %Frisch ela st. from Mendoza 1991
OLD_OMEGA = OMEGA;
DBAR = 0.74421765717098; %debt
OLD_DBAR = DBAR;
PSSI = 0.11135/150; %debt elasticity of interest rate
OLD_PSSI = PSSI;
PHI = 0.028; %capital adjustment cost
OLD_PHI = PHI;
RHO = 0.42; %persistence of TFP shock
OLD_RHO = RHO;
ETATILDE = 0.0129; %standard deviation of innovation to TFP shock
OLD_ETATILDE = ETATILDE;

%target tb/y,sigma_y,rho_y,sigma_i,rho_i,sigma_h
target = [0.02;2.81;0.62;9.82;0.31;2.02];
threshold = 0.1;
current_error = 1;

while current_error > threshold
    param = randi(6);
    if param == 1
        OLD_OMEGA = OMEGA;
        OMEGA = perturb(OMEGA);
    elseif param == 2
        OLD_DBAR = DBAR;
        DBAR = perturb(DBAR);
    elseif param == 3
        OLD_PSSI = PSSI;
        PSSI = perturb(PSSI);
    elseif param == 4
        OLD_PHI = PHI;
        PHI = perturb(PHI);
    elseif param == 5
        OLD_RHO = RHO;
        RHO = perturb(RHO);
    elseif param == 6
        OLD_ETATILDE = ETATILDE;
        ETATILDE = perturb(ETATILDE);
    end
    [stds, scorrm actual] = iterate(OMEGA,DBAR,PSSI,PHI,RHO,ETATILDE);
    distance = compute_distance(actual, target);
    distance
    if distance > current_error
        if param == 1
            OMEGA = OLD_OMEGA;
        elseif param == 2
            DBAR = OLD_DBAR;
        elseif param == 3
            PSSI = OLD_PSSI;
        elseif param == 4
            PHI = OLD_PHI;
        elseif param == 5
            RHO = OLD_RHO;
        elseif param == 6
            ETATILDE = OLD_ETATILDE;
        end
    else
        current_error = distance;
    end
end

function new_value = perturb(old_value) 
    val = randi(2);
    PERTURBATION_PERCENT = 0.10;
    if val == 1
        new_value = old_value + PERTURBATION_PERCENT * old_value;
    else
        new_value = old_value - PERTURBATION_PERCENT * old_value;
    end
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
    [sigy0,sigx0]=mom(gx,hx,varshock);
    stds = sqrt(diag(sigy0));
    
    %serial correlations
    [sigy1,sigx1]=mom(gx,hx,varshock,1);
    scorr = diag(sigy1)./diag(sigy0);
    targets = [0.02;stds(noutput); scorr(noutput); stds(nivv); scorr(nivv);stds(nh)];
end
