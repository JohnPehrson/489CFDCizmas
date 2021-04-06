function [cells_q,cells_f,cells_g] = setInitialConditions(user_Mach,user_Gamma,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax)
%This function defines the initial conditions of the flow
%Assume 2 layers of ghost nodes
%Using info from page 29 of my handwritten notes

rho = 1;
alpha = 0;
u = user_Mach*cos(alpha);
v = user_Mach*sin(alpha);
p = P_static; %can forcefully set p = pstatic because we are also forcefully setting mach number
E = p/((user_Gamma-1))+0.5*user_Mach^2;
H = E+p;

%assume that all cells have the same initial conditions
q = [rho;rho*u;rho*v;rho*E];
f = [rho*u;rho*u^2+p;rho*u*v;rho*H*u];
g = [rho*v;rho*u*v;rho*v^2+p;rho*H*v];

    %assign every interior cell the initial condition
    for i = 3:cells_Imax-2
        for j = 3:cells_Jmax-2
    cells_q(i,j,:) = q;
    cells_f(i,j,:) = f;
    cells_g(i,j,:) = g;
        end
    end

end

