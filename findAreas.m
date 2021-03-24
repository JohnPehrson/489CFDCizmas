function [A] = findAreas(nodes_x,nodes_y,cells_Imax,cells_Jmax,A)
%Finds the area of each cell based on the node positions that define the
%cell. This is done only once becasue this is a static grid

    %double loop to sweep through all cells
    for i = 1:cells_Imax
        for j = 1:cells_Jmax
            [X_abcd,Y_abcd] = nodes_touch_cell(i,j,nodes_x,nodes_y);
            A(i,j) = 0.5*((X_abcd(3)-X_abcd(1))*(Y_abcd(4)-Y_abcd(2))-(X_abcd(4)-X_abcd(2))*(Y_abcd(3)-Y_abcd(1)));
        end
    end
end

