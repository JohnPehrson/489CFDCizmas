function [f_cell_gc,g_cell_gc] = RK_wallBC(user_Gamma,q_cell,f_cell,g_cell,walltan)
%This function applies a wall boundary condition to a nearby cell if that
%cell has been identified as a ghost cell. This function is only called in
%the RK loop.

%data from flow cell
rho = q_cell(1);
u = q_cell(2)/rho;
v = q_cell(3)/rho;
p = f_cell(2)-rho*u^2;
e = p/((user_Gamma-1)*rho);

%find alpha for velocity switch
flowcellvelocity = [u;v];
% angle = acosd(dot(flowcellvelocity/norm(flowcellvelocity),walltan/norm(walltan)));

%new method
beta = atand(walltan(2)/walltan(1));
gamma = atand(flowcellvelocity(2)/flowcellvelocity(1));
alpha = gamma-beta;
theta = beta-alpha;

%find mirrored velocities
velocity_gc = [cosd(theta);sind(theta)]*norm(flowcellvelocity);
u_gc = velocity_gc(1);
v_gc = velocity_gc(2);
E_gc = e+0.5*(u_gc^2+v_gc^2);
H_gc = E_gc+p/rho;

%fill out ghost cell q,f,g
% q_cell_gc = [rho;rho*u_gc;rho*v_gc;rho*E_gc];
f_cell_gc = [rho*u_gc;rho*u_gc^2+p;rho*u_gc*v_gc;rho*H_gc*u_gc];
g_cell_gc = [rho*v_gc;rho*u_gc*v_gc;rho*v_gc^2+p;rho*H_gc*v_gc];

end

