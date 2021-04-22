function [cells_p,cells_c,cells_eig] = update_p_c(user_Gamma,cells_p,cells_c,cells_eig,cells_q,cells_f,nodes_x,nodes_y,cells_Imax,cells_Jmax)
%This function updates the p,c,eigenvalue matraxies for the grid. This
%function should be called after the boundary conditions are reapplied

    for i = 1:cells_Imax
        for j = 1:cells_Jmax
            f = cells_f(i,j,:);
            q = cells_q(i,j,:);
            cells_p(i,j) = f(2)-(q(2)^2)/q(1);
            cells_c(i,j) = sqrt(user_Gamma*cells_p(i,j)/q(1));
            [x_abcd,y_abcd] = nodes_touch_cell(i,j,nodes_x,nodes_y);
            [cells_eig(i,j,1),cells_eig(i,j,2)] = eigenvaluefinder(x_abcd,y_abcd,q,cells_c(i,j));

        end
    end

end

