function [f_gc,g_gc] = RK_outletBC(user_Gamma,user_Mach,P_static,q_in,f_in,g_in,f_in_adj,g_in_adj)
%This is a function that modifies f and g in the RK loop, to be called when
%the RK cell is bordering the outlet. The f and g for the ghost cell are
%sent back to the RK loop.


%create vectors for rho,m,n for the cells, index 3 is ghost cells
m = [f_in_adj(1),q_in(2),NaN];
n = [g_in_adj(1),q_in(3),NaN];
uadj = g_in_adj(2)/g_in_adj(1);
rho = [f_in_adj(1)/uadj,q_in(1),NaN];

%calculate new ghost cells rho,m,n, shared between both ghost cells
rho(3) = 2*rho(2)-rho(1);
m(3) = 2*m(2)-m(1);
n(3) = 2*n(2)-n(1);

%calculate new epsilon for ghost cells
u = m(3)/rho(3);
v = n(3)/rho(3);


% epsilon = (P_static)/(user_Gamma-1)+0.5*(u^2+v^2)*(rho(3));
E = (P_static)/((user_Gamma-1)*rho(3))+0.5*(u^2+v^2);

%make single vectors
% q_out1cell = [rho(3),m(3),n(3),rho(3)*E];
f_gc = [m(3),(m(3)^2)/rho(3)+P_static,m(3)*n(3)/rho(3),m(3)*(E+P_static/rho(3))];
g_gc = [n(3),m(3)*n(3)/rho(3),(n(3)^2)/rho(3)+P_static,n(3)*(E+P_static/rho(3))];




end

