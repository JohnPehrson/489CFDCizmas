function [cells_q,cells_f,cells_g] = applyOutletBC(user_Gamma,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax)
%This function applies the outlet boundary conditions by modifying the
%ghost nodes outside the right edge of the grid. 

%rho,m=rho*u,and n=rho*v are found from a second order approximation from
%inside the grid.
%Pstatic_inf is defined (right now as 1/gamma)
%epsilon = rho*E is defined using the static pressure, m,n

%reference indexes for grid cells and ghost cells
i_g = (cells_Imax-3):(cells_Imax-2);
i_gc = (cells_Imax-1):cells_Imax;

%loop through vertical cells at the outlet
    for j = 3:(cells_Jmax-2) 
        [cells_q(i_gc,j,:),cells_f(i_gc,j,:),cells_g(i_gc,j,:)] = singleCellOutletBC(user_Gamma,cells_q(i_g,j,:),cells_f(i_g,j,:),cells_g(i_g,j,:),P_static);
        %[q_out,f_out,g_out] = singleCellOutletBC(gammag,q_in,f_in,g_in,Pstatic_inf)
    end


end

