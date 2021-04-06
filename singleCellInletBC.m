function [q_b,f_b,g_b] = singleCellInletBC(user_Gamma,user_Mach,P_resevoir,alpha,q_in,f_in,g_in)
%This sub-function sets the ghost cells outside the inlet of the grid.
%The program uses Riemann invariants, pressure, and flow AoA to find the
%ghost cells.
%alpha is in degrees

%q_in, f_in, and g_in all reference the first cell inside the wall
%q_b, f_b, and g_b reference the TWO ghost cells that act as the boundary


normalinlet = [-1,0]; %outward pointing normal

%Calculate R-, referred to as R_gc (ghost cell) -data comes from ref/inlet
P_gc = pressurefinder(P_resevoir,user_Mach,user_Gamma);
rho_gc = 1; %by nondimensionalization
c_gc = sqrt(user_Gamma*P_gc/rho_gc);
V_gc = user_Mach*[cosd(alpha),sind(alpha)];
R_gc = dot(V_gc,normalinlet)-(2*c_gc)/(user_Gamma-1);

%Calculate R+, referred to as R_in (interior of the computational grid)
%-data comes from the interior cells in q_in,f_in,g_in
V_in = [q_in(2)/q_in(1),q_in(3)/q_in(1)];
rho_in = q_in(1);
P_in = pressurefinder(P_resevoir,norm(V_in),user_Gamma);
c_in = sqrt(user_Gamma*P_in/rho_in);
R_in = dot(V_in,normalinlet)+(2*c_in)/(user_Gamma-1);

%Calculate Vb and cb
Vb = (R_in+R_gc)/2;
c_b = (user_Gamma-1)*(R_in-R_gc)/4;

%Vb vector calculated using
Vbvec = V_gc + (Vb-dot(V_gc,normalinlet))*normalinlet;
u_b = Vbvec(1);
v_b = Vbvec(2);

%Finding variables to put into q,f,g boundary cells
s_b = (c_gc^2)/(user_Gamma*rho_gc^(user_Gamma-1));
rho_b = (c_b^2)/(user_Gamma*s_b);
p_b = rho_b*(c_b^2)/user_Gamma;
E_b = p_b/(rho_b*(user_Gamma-1))+0.5*Vb^2;

%Put variables into the q,f,g form
q_b_cell = [rho_b;rho_b*u_b;rho_b*v_b;rho_b*E_b];
f_b_cell = [rho_b*u_b;rho_b*u_b^2+p_b;rho_b*u_b*v_b;rho_b*(E_b+p_b)*u_b];
g_b_cell = [rho_b*v_b;rho_b*u_b*v_b;rho_b*v_b^2+p_b;rho_b*(E_b+p_b)*v_b];

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

