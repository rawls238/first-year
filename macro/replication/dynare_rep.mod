var y z k n c lambda;
varexo e;
parameters mu beta alpha gamma theta v sigma delta rho; 
 
mu    = 0.34;
beta = 0.99;
alpha = 1.0;
gamma = -1.0;
theta = 0.36;
v = 3.0;
sigma = 0.0099;
delta = 0.025;
A = [.906 .088; .088 .906];

model;
y_h = ((lambda*(k_h)^theta*(n_h)^(1-theta))^(-v)+sigma*(z_h)^(-v))^(-1/v);
y_f = ((lambda*(k_f)^theta*(n_f)^(1-theta))^(-v)+sigma*(z_f)^(-v))^(-1/v);
mu*c_h^(gamma*mu-1)*(1-n_h)^(gamma*(1-mu)) * (1-theta)*y_h^(v+1)*(y_h^(-v)-sigma*z_h^(-v))/n_h=-(mu-1)*c_h^(gamma*mu)*(1-n_h)^(gamma*(1-mu)-1);
mu*c_f^(gamma*mu-1)*(1-n_f)^(gamma*(1-mu)) * (1-theta)*y_f^(v+1)*(y_f^(-v)-sigma*z_f^(-v))/n_f=-(mu-1)*c_f^(gamma*mu)*(1-n_h)^(gamma*(1-mu)-1);
mu*c_h^(gamma*mu-1)*(1-n_h)^(gamma*(1-mu))=beta*((1-delta) + theta*y_h(+1)^(v+1)*(y_h(+1)^(-v)-sigma*z_h(+1)^(-v))/k_h(+1)) * mu*c_h(+1)^(gamma*mu-1)*(1-n_h(+1))^(gamma*(1-mu));
mu*c_f^(gamma*mu-1)*(1-n_f)^(gamma*(1-mu))=beta*((1-delta) + theta*y_f(+1)^(v+1)*(y_f(+1)^(-v)-sigma*z_f(+1)^(-v))/k_f(+1)) * mu*c_f(+1)^(gamma*mu-1)*(1-n_f(+1))^(gamma*(1-mu));
mu*c_h^(gamma*mu-1)*(1-n_h)^(gamma*(1-mu))=beta * mu*c_h(+1)^(gamma*mu-1)*(1-n_h(+1))^(gamma*(1-mu)) * (1+sigma*(y_h(+1)^(v+1))*z_h(+1)^(-v-1));
mu*c_f^(gamma*mu-1)*(1-n_f)^(gamma*(1-mu))=beta * mu*c_f(+1)^(gamma*mu-1)*(1-n_f(+1))^(gamma*(1-mu)) * (1+sigma*(y_h(+1)^(v+1))*z_h(+1)^(-v-1));
y_h = c_h + k_h(+1)-(1-delta)*k_h+z_h(+1)-z_h;
y_f = c_f + k_f(+1)-(1-delta)*k_f+z_f(+1)-z_f;
lambda_h = A(1,2)*lambda_h(-1)+ A(1,2)*lambda_h + e;
lambda_f = A(2,1) * lambda_h(-1) + A(2,2)*lambda_f + e;
end;

initval;
k = 1;
n = 0.5;
c = 1;
y = 1;
z = 1;
lambda = 1;
e = 0;
end;


shocks;
var E_H; stderr 0.00852;
var E_F; stderr 0.00852;
corr E_H, E_F = 0.258;
end;


steady;
check;

stoch_simul(order=1, hp_filter=1600);
