function [f,g] = refresh_f_g(q,user_Gamma)
%This function refreshes the cellular f,g vectors given an input vector q
f = zeros(4,1);
g = f;

rho = q(1);
u = q(2)/q(1);
v = q(3)/q(1);
E = q(4)/q(1);
p = (user_Gamma-1)*rho*(E-0.5*(u^2+v^2));

f = [rho*u;rho*u^2+p;rho*u*v;rho*(E+p/rho)*u];
g = [rho*v;rho*u*v;rho*v^2+p;rho*(E+p/rho)*v];
end

