function [cells_p,cells_c] = update_p_c(user_Gamma,cells_p,cells_c,cells_q,cells_f,cells_Imax,cells_Jmax)
%This function updates the p,c,eigenvalue matraxies for the grid. This
%function should be called after the boundary conditions are reapplied

    for i = 1:cells_Imax
        for j = 1:cells_Jmax
            f = cells_f(i,j,:);
            q = cells_q(i,j,:);
            cells_p(i,j) = f(2)-(q(2)^2)/q(1);
            cells_c(i,j) = sqrt(user_Gamma*cells_p(i,j)/q(1));
        end
    end

end

