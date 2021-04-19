function [q_cell_gc,f_cell_gc,g_cell_gc] = singleCellWallBC(q_cell,f_cell,g_cell,walltan)
%This function applies a wall boundary condition on a single ghost cell
%when given inputs of the flow cell
%_cell is flow cell, _cell_gc is a ghost cell

%data from flow cell
rho = q_cell(1);
u = q_cell(2)/q_cell(1);
v = q_cell(3)/q_cell(1);
p = f_cell(2)-rho*u^2;
e = q_cell(4)/q_cell(1)-0.5*(u^2+v^2);

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
H_gc = E_gc+p;

%fill out ghost cell q,f,g
q_cell_gc = [rho;rho*u_gc;rho*v_gc;rho*E_gc];
f_cell_gc = [rho*u_gc;rho*u_gc^2+p;rho*u_gc*v_gc;rho*H_gc*u_gc];
g_cell_gc = [rho*v_gc;rho*u_gc*v_gc;rho*v_gc^2+p;rho*H_gc*v_gc];


% %new method
% beta = atand(walltan(2)/walltan(1));
% gamma = atand(flowcellvelocity(2)/flowcellvelocity(1));
% alpha = gamma-beta;
% theta = beta-alpha;
% angle = theta; %to get into the old format

end

