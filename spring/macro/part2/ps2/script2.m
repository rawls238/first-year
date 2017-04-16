clear all;
close all;

rng(6216);

rho      = 0.66;
sigma_epsi  = 0.08;
mu       = 0.0021;
sigma_eta  = 0.0032;

beta = .96^(1/12);
K = 0.0245;
num_A    = 31;
sig_A  = sigma_epsi/sqrt(1-rho^2);
a = linspace(-3*sig_A,3*sig_A,num_A)';
num_P    = 301;
rp = linspace(-.3,.3,num_P)';

naa = length(a);
P_a = NaN(naa, naa);
for i = 1:naa
    for j = 1:naa
        if j==1
            P_a(i,j) = normcdf((a(j)+(a(j+1)-a(j))/2-rho*a(i))/sigma_epsi,0,1);
        elseif j==naa
            P_a(i,j) = 1-normcdf((a(j)-(a(j)-a(j-1))/2-rho*a(i))/sigma_epsi,0,1);
        else
            first = normcdf((a(j)+(a(j)-a(j-1))/2-rho*a(i))/sigma_epsi,0,1);
            second = normcdf((a(j)-(a(j)-a(j-1))/2-rho*a(i))/sigma_epsi,0,1);
            P_a(i,j) = first - second;
        end
    end
end


nrp = length(rp);
P_rp = NaN(nrp, nrp);
for i = 1: nrp
    for j = 1: nrp
        if j==1
            P_rp(i,j) = normcdf((rp(j)+(rp(j+1)-rp(j))/2+mu-rp(i))/sigma_eta,0,1);
        elseif j==nrp
            P_rp(i,j) = 1-normcdf((rp(j)-(rp(j)-rp(j-1))/2+mu-rp(i))/sigma_eta,0,1);
        else
            first = normcdf((rp(j)+(rp(j)-rp(j-1))/2+mu-rp(i))/sigma_eta,0,1);
            second = normcdf((rp(j)-(rp(j)-rp(j-1))/2+mu-rp(i))/sigma_eta,0,1);
            P_rp(i,j) = first - second;
        end
    end
end

theta = 4;
C = 1;
first = exp(rp).^(-theta)*ones(1, num_A);
second = exp(rp)*ones(1,num_A)- ones(num_P,1)*((theta-1)/theta)*(exp(a)'.^(-1));
profits = C * first .* second;

prev_v = zeros(num_P,num_A);

tolerance = 1e-6;
distance = 1;
while distance>tolerance
    prev_val = profits+beta*P_rp*prev_v*(P_a');
    potential_cur_val= ones(num_P,1)*max(prev_val-(theta-1)/theta*K,[],1);
    cur_v   = max(prev_val,potential_cur_val);
    distance = max(max(abs(cur_v-prev_v)));
    prev_v = cur_v;
end

[~, delta_pol] = max(prev_val-(theta-1)/theta*K,[],1);
delta_pol = ones(num_P,1)*delta_pol;
no_delta_pol = kron((1:num_P)',ones(1,num_A));
diff = (potential_cur_val>prev_val);
policy_fn = diff.*delta_pol+(ones(size(diff))-diff).*no_delta_pol;

figure;
surf(rp,a,rp(policy_fn)');
xlabel('log(p/P)');
ylabel('log(A)');
zlabel('Policy value');
figure;
surf(rp,a,cur_v');
xlabel('log(p/P)');
ylabel('log(A)');
zlabel('Value function');


nrp = length(rp);
grid = rp;
P_inf  = NaN(1, nrp);
for k = 1: nrp
    if k==1
        P_inf(1,k) = normcdf((rp(k)+(rp(k+1)-rp(k))/2-mu)/sigma_eta,0,1);
    elseif k==nrp
        P_inf(1,k) = 1-normcdf((grid(k)-(rp(k)-rp(k-1))/2-mu)/sigma_eta,0,1);
    else
        first = normcdf((rp(k)+(rp(k)-rp(k-1))/2-mu)/sigma_eta,0,1);
        second = normcdf((rp(k)-(rp(k)-rp(k-1))/2-mu)/sigma_eta,0,1);
        P_inf(1,k) = first - second;
    end
end
P_inf = ones(nrp,1)*P_inf;


num_periods = 200;
init_a = ceil(num_A/2);
init_p = ceil(num_P/2);

chain = markovchain(P_a,num_periods,init_a);
infl = markovchain(P_inf,num_periods,init_p);
log_pt = NaN(num_periods,1);
log_Pt = NaN(num_periods,1);
log_mc = NaN(num_periods,1);

log_Pt(1) = 0;
pol_sim = NaN(num_periods,1);
prob = NaN(num_periods,1);
prob(1) = init_p;


for i=1:num_periods
    if i==1
        log_Pt(i)= rp(infl(i));
    else
        log_Pt(i) = log_Pt(i-1) + rp(infl(i));
        p_st  = log_pt(i-1)-log_Pt(i);
        state1 = (rp<=p_st+eps);
        state2 = (rp>=p_st-eps);
        q = ((state1+state2)==2);
       if sum(q)>0
           prob(i) = find((state1+state2)==2);
       else
           if sum(state1)==num_P
               prob(i) = num_P;
           else
               prob(i) = 1;
           end
       end
    end
    pol_sim(i) = policy_fn(prob(i),chain(i));
    log_pt(i)  = rp(pol_sim(i))+log_Pt(i);
    log_mc(i)  = log((theta-1)/theta)+log_Pt(i)-a(chain(i));
end

h = figure(); 
plot(1:num_periods,log_pt,'--');
hold on;
plot(1:num_periods,log_Pt,'-.');
hold on;
plot(1:num_periods,log_mc,'r-');
hold on;
legend('log(p)','log(P)','log(mc)');

function S = markovchain(P,T,S0)
    [n,nc] = size(P);
    S      = zeros(T,1);
    if nargin < 3
        S(1,1) = ceil(rand*n);
    else
        S(1,1) = S0;
    end;
    C = cumsum(P');
    for t = 2:T
        j  = find(rand < C(:,S(t-1)));
    S(t,1) = j(1);
    end;
end