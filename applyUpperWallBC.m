function [q,f,g] = applyUpperWallBC(nodes_x,nodes_y,q,f,g,cells_Imax,cells_Jmax)
%This function applies boundary conditions on the upper wall through the
%use of ghost cells. The ghost cells have their pressure, energy, and
%density matched to interior cells, and their velocities mirrored over the
%wall. 

%This function uses the x and y locations of the nodes to find the wall
%direction. The function then uses the q,f,g vectors to identify the
%necessary ghost cells, and then reports the new
%q,f,g values with modified ghost nodes. 

%q,f,g will be Imax,Jmax,4 matraxies each, representing the whole grid
%cells_Imax and cells_Jmax represent the number of cells in each direction

%loop through ghost cells on the bottom of the wall, make changes
j_gc1 = cells_Jmax; %The "last" row of ghost cells
j_gc2 = cells_Jmax-1; %The "second to last" row of ghost cells
    for i = 3:(cells_Imax-2) %aka for main cells in the i direction
        %find wall tangent vector
            [X_wall_abcd,Y_wall_abcd] = nodes_touch_cell(i,j_gc2-1,nodes_x,nodes_y);
            walltan = [X_wall_abcd(3)-X_wall_abcd(4);Y_wall_abcd(3)-Y_wall_abcd(4)];
        %do boundary calculations to change the second row, by the wall
            [q(i,j_gc2,:),f(i,j_gc2,:),g(i,j_gc2,:)] = singleCellWallBC(q(i,j_gc2-1,:),f(i,j_gc2-1,:),g(i,j_gc2-1,:),walltan);            
            
        %do boundary calculations to change the first row
            [q(i,j_gc1,:),f(i,j_gc1,:),g(i,j_gc1,:)] = singleCellWallBC(q(i,j_gc1-3,:),f(i,j_gc1-3,:),g(i,j_gc1-3,:),walltan);
                       
    end


end

