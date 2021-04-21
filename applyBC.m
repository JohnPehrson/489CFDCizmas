function [cells_q,cells_f,cells_g] = applyBC(nodes_x,nodes_y,user_alpha,user_Gamma,user_Mach,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax)
%This function contains the 4 boundary condition functions. Call this
%function to refresh all boundary conditions. 
    [cells_q,cells_f,cells_g] = applyUpperWallBC(nodes_x,nodes_y,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    [cells_q,cells_f,cells_g] = applyBottomWallBC(nodes_x,nodes_y,cells_q,cells_f,cells_g,cells_Imax);
    [cells_q,cells_f,cells_g] = applyOutletBC(user_Gamma,user_Mach,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax,P_static);
    [cells_q,cells_f,cells_g] = applyInletBC(user_alpha,user_Gamma,user_Mach,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
end

