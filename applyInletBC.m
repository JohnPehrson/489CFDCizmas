function [cells_q,cells_f,cells_g] = applyInletBC(alpha,user_Gamma,user_Mach,P_resevoir,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax)
%This function applies the inlet boundary conditions by modifying the
%ghost nodes outside the left edge of the grid. 

%This function just loops through the various cells on the boundary and
%passes necessary data down the line.
%P_resevoir is used, along with Mach, to find static pressure

%reference indexes for grid cells and ghost cells
i_g = 3;
i_gc = 1:2;

%loop through vertical cells at the inlet
    for j = 3:(cells_Jmax-2) 
        [cells_q(i_gc,j,:),cells_f(i_gc,j,:),cells_g(i_gc,j,:)] = singleCellInletBC(user_Gamma,user_Mach,P_resevoir,alpha,cells_q(i_g,j,:),cells_f(i_g,j,:),cells_g(i_g,j,:));
        %[q_b,f_b,g_b] = singleCellInletBC(user_Gamma,user_Mach,P_resevoir,alpha,q_in,f_in,g_in)
    end


end

