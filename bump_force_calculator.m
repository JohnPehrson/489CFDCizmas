function [F] = bump_force_calculator(cells_pressure,nodes_x,nodes_y,cells_range_i)
%This function calculates the force exerted on the bump using the cell
%pressures. The node locations are used to find the cell edge length on the
%bump and to find the normal vector to the wall. The cells_range_i is the
%range of i values of cells that have their bottom face touching the bump.

%This function outputs the force in the x and y direction acting on the
%bump. 
j = 3; %bottom most row of cells in the grid
F = [0;0];

    for i = cells_range_i(1):cells_range_i(end) %loop over bump cells

        %find wall normal
        [X_abcd,Y_abcd] = nodes_touch_cell(i,j,nodes_x,nodes_y);
        [nx,ny] = cellnormal(X_abcd(1),X_abcd(2),Y_abcd(1),Y_abcd(2));
        n = [nx;ny];
        %find wall edge length
        [length] = edgelength(X_abcd,Y_abcd,'S');
        %find force in the wall direction
        Fcell = -cells_pressure(i,j)*n*length;
        %sum force
        F = F+Fcell;
    end

end

