function [q_b,f_b,g_b] = singleCellInletBC(user_Gamma,user_Mach,P_static,alpha,q_in,f_in,g_in)
%This sub-function sets the ghost cells outside the inlet of the grid.
%The program uses Riemann invariants, pressure, and flow AoA to find the
%ghost cells.
%alpha is in degrees

%q_in, f_in, and g_in all reference the first cell inside the wall
%q_b, f_b, and g_b reference the TWO ghost cells that act as the boundary

%Find P_tot,rho_tot in the inlet far field, use to find R1
P_inf = P_static; % pressurefinder(user_Mach,user_Gamma,P_static);
rho_inf = 1; %1/((1+(user_Gamma-1)/2*user_Mach^2)^(-1/(user_Gamma-1))); %assume far field rho is 1
R1 = (2/(user_Gamma-1))*sqrt(user_Gamma*P_inf/rho_inf)+user_Mach; %Riemann invariant 1 in the farfield where there is no flow

%squeeze flow parameters
q = squeeze(q_in);
f = squeeze(f_in);
g = squeeze(g_in);
%Find R2, invariant from the flow field
uflow = q(2)/q(1);
vflow = q(3)/q(1);
Vflowmag = sqrt(uflow^2+vflow^2);
p = f(2)-(q(2)^2)/q(1);
c = sqrt(user_Gamma*p/q(1));
R2 = Vflowmag-2*c/(user_Gamma-1);

%find v at inlet using riemans
Vinletmag = 0.5*(R1+R2);
uinlet = Vinletmag*cos(alpha);
vinlet = Vinletmag*sin(alpha);
cinlet = 0.25*(user_Gamma-1)*(R1-R2);
Minlet = Vinletmag/c;
Pinlet = P_static;
rhoinlet = user_Gamma*Pinlet/(cinlet^2);
Einlet = Pinlet/((user_Gamma-1)*rhoinlet)+0.5*(uinlet^2+vinlet^2);

%Put variables into the q,f,g form
q_b_cell = [rhoinlet;rhoinlet*uinlet;rhoinlet*vinlet;rhoinlet*Einlet];
f_b_cell = [rhoinlet*uinlet;rhoinlet*uinlet^2+Pinlet;rhoinlet*uinlet*vinlet;rhoinlet*(Einlet+Pinlet)*uinlet];
g_b_cell = [rhoinlet*vinlet;rhoinlet*uinlet*vinlet;rhoinlet*vinlet^2+Pinlet;rhoinlet*(Einlet+Pinlet)*vinlet];

%create empty output matrixes in the right format
q_b = NaN(2,1,4);
f_b = q_b;
g_b = q_b;

%convert the gc into the right format for output
q_b(1,:,:) = q_b_cell;
q_b(2,:,:) = q_b_cell;
f_b(1,:,:) = f_b_cell;
f_b(2,:,:) = f_b_cell;
g_b(1,:,:) = g_b_cell;
g_b(2,:,:) = g_b_cell;


end

