function [q_out,f_out,g_out] = singleCellOutletBC(gammag,q_in,f_in,g_in,Pstatic_inf)
%This function sets the outlet boundary conditons in a single cell. Used in
%"applyOutletBC".

%q_in,f_in, and g_in are all 2x1x4 vectors, and contain the information for
%both the right-most and second-right-most cells in the flow

%q_out,f_out,g_out are also 2x1x4 and represent the ghost cells in the same
%row as the inputs

%create vectors for rho,m,n for the cells, index 3 is ghost cells
rho = [q_in(1,1,1),q_in(2,1,1),NaN];
m = [q_in(1,1,2),q_in(2,1,2),NaN];
n = [q_in(1,1,3),q_in(2,1,3),NaN];

%calculate new ghost cells rho,m,n, shared between both ghost cells
rho(3) = 2*rho(2)-rho(1);
m(3) = 2*m(2)-m(1);
n(3) = 2*n(2)-n(1);

%calculate new epsilon for ghost cells
epsilon = (Pstatic_inf)/(gammag-1)+0.5*(m(3)^2+n(3)^2)/(rho(3));

%make single vectors
q_out1cell = [rho(3),m(3),n(3),epsilon];
f_out1cell = [m(3),(m(3)^2)/rho(3)+Pstatic_inf,m(3)*n(3)/rho(3),m(3)*((epsilon/rho(3))+Pstatic_inf)];
g_out1cell = [n(3),m(3)*n(3)/rho(3),(n(3)^2)/rho(3)+Pstatic_inf,n(3)*((epsilon/rho(3))+Pstatic_inf)];

%create empty output matrixes in the right format
q_out = NaN(2,1,4);
f_out = q_out;
g_out = q_out;

%convert the gc into the right format for output
q_out(1,:,:) = q_out1cell;
q_out(2,:,:) = q_out1cell;
f_out(1,:,:) = f_out1cell;
f_out(2,:,:) = f_out1cell;
g_out(1,:,:) = g_out1cell;
g_out(2,:,:) = g_out1cell;


end

