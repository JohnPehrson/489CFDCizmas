function [q,f,g] = applyBottomWallBC(nodes_x,nodes_y,q,f,g,cells_Imax)
%This function applies boundary conditions on the bottom wall through the
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
j_gc1 = 1; %The "first" row of ghost nodes is at j = 1
j_gc2 = 2; %The "second" row of ghost nodes is at j = 2
    for i = 3:(cells_Imax-2) %aka for main cells in the i direction
        %find wall tangent vector
            [X_wall_abcd,Y_wall_abcd] = nodes_touch_cell(i,j_gc2+1,nodes_x,nodes_y);
            walltan = [X_wall_abcd(2)-X_wall_abcd(1);Y_wall_abcd(2)-Y_wall_abcd(1)];
        %do boundary calculations to change the second row, by the wall
            [q(i,j_gc2,:),f(i,j_gc2,:),g(i,j_gc2,:)] = singleCellWallBC(q(i,j_gc2+1,:),f(i,j_gc2+1,:),g(i,j_gc2+1,:),walltan);            
            
        %do boundary calculations to change the first row
            [q(i,j_gc1,:),f(i,j_gc1,:),g(i,j_gc1,:)] = singleCellWallBC(q(i,j_gc1+3,:),f(i,j_gc1+3,:),g(i,j_gc1+3,:),walltan);
                       
    end


end

