function [eigenvalueout] = eigenvalueDis(c,i,j,x,y,q,face)
%This function calculates the absolute value of the largest eigenvalue of
%the euler equation aligned with the xi/n directions. 
%c is speed of sound, i and j are the cell indicies, x and y are the whole
%x and y vector of nodes
%face determines what face we are looking at: 'N','E','S','W'
[X_abcd,Y_abcd] = nodes_touch_cell(i,j,x,y);
velomagcell = sqrt((q(i,j,2)^2+q(i,j,3)^2)/(q(i,j,1)^2)); %abs(u) at the cell of interest
flowdirdegcell = atand(q(i,j,3)/q(i,j,2)); %finds the direction of flow in the main cell

    switch face
        case 'N'
            velomagadj = sqrt((q(i,j+1,2)^2+q(i,j+1,3)^2)/(q(i,j+1,1)^2)); %abs(u) for the other cell bordering the edge
            flowdirdegadj = atand(q(i,j+1,3)/q(i,j+1,2)); %finds the direction of flow in the adjacent cell
             Startpoint = [X_abcd(4), Y_abcd(4)];
             Endpoint = [X_abcd(3), Y_abcd(3)];
             direction = Endpoint - Startpoint;
             direction_scaled = direction/norm(direction); % vector with length 1
             normal = [-direction_scaled(2),direction_scaled(1)]; %finds outward pointing normal vector
             
        case 'E'
            velomagadj = sqrt((q(i+1,j,2)^2+q(i+1,j,3)^2)/(q(i+1,j,1)^2)); %abs(u) for the other cell bordering the edge
            flowdirdegadj = atand(q(i+1,j,3)/q(i+1,j,2));%finds the direction of flow in the adjacent cell
             Startpoint = [X_abcd(2), Y_abcd(2)];
             Endpoint = [X_abcd(3), Y_abcd(3)];
             direction = Endpoint - Startpoint;
             direction_scaled = direction/norm(direction); % vector with length 1
             normal = [direction_scaled(2),-direction_scaled(1)]; %finds outward pointing normal vector
        case 'S'
            velomagadj = sqrt((q(i,j-1,2)^2+q(i,j-1,3)^2)/(q(i,j-1,1)^2)); %abs(u) for the other cell bordering the edge
            flowdirdegadj = atand(q(i,j-1,3)/q(i,j-1,2));%finds the direction of flow in the adjacent cell
             Startpoint = [X_abcd(1), Y_abcd(1)];
             Endpoint = [X_abcd(2), Y_abcd(2)];
             direction = Endpoint - Startpoint;
             direction_scaled = direction/norm(direction); % vector with length 1
             normal = [direction_scaled(2),-direction_scaled(1)]; %finds outward pointing normal vector
        case 'W'
            velomagadj = sqrt((q(i-1,j,2)^2+q(i-1,j,3)^2)/(q(i-1,j,1)^2)); %abs(u) for the other cell bordering the edge
            flowdirdegadj = atand(q(i-1,j,3)/q(i-1,j,2));%finds the direction of flow in the adjacent cell
             Startpoint = [X_abcd(1), Y_abcd(1)];
             Endpoint = [X_abcd(4), Y_abcd(4)];
             direction = Endpoint - Startpoint;
             direction_scaled = direction/norm(direction); % vector with length 1
             normal = [-direction_scaled(2),direction_scaled(1)]; %finds outward pointing normal vector
    end
            velomag = 1/2*(velomagcell+velomagadj); %average of cell velocity magnitudes
            flowdirdeg = 1/2*(flowdirdegadj+flowdirdegcell); %average of cell velocity directions
            velo = velomag*[cosd(flowdirdeg),sind(flowdirdeg)];  %average velocity in components
            
            eigenvalueout = dot(velo,normal)+c; 
end

