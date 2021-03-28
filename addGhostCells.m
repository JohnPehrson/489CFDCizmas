function [node_x_output,node_y_output] = addGhostCells(nodes_x_input,nodes_y_input)
%This function adds 2 rows/columns of ghost cells to every boundary,
%required for the current implementation of boundary conditions

%% Measure data and port to the expanded matrix
%identify input grid size Imax and Jmax
Imax = length(nodes_x_input(:,1))+4; %+4 because of the new nodes
Jmax = length(nodes_x_input(1,:))+4; %+4 because of the new nodes

%initialize output matraxies
nodes_x_output = NaN(Imax,Jmax);
nodes_y_output = NaN(Imax,Jmax);

%port input data into the center of the output matraxies
nodes_x_output(3:Imax-2,3:Jmax-2) = nodes_x_input;
nodes_y_output(3:Imax-2,3:Jmax-2) = nodes_y_input;

%% Add ghost nodes
    %assume corners are square, ghost nodes are a 'continuation', 16
    %special cases that we do more-or less manually because they don't have
    %attached edges
    [nodes_x_output,nodes_y_output] = GhostNodesCorners(nodes_x_output,nodes_y_output,Imax,Jmax);
    
    %for writing convenience
    x = nodes_x_output;
    y = nodes_y_output;
    
    %do the top face, two rows
    for celloop = 2:-1:1
        jrow = Jmax-celloop; %so jrow+1 is the new nodes
        for i = 4:Imax-3  %only consider nodes that have two attached faces 
            [x(i,jrow+1),y(i,jrow+1)] = addGhostNode(x(i-1,jrow),y(i-1,jrow),x(i,jrow),y(i,jrow),x(i+1,jrow),y(i+1,jrow),x(i,jrow-1),y(i,jrow-1));
            %[gx,gy] = addGhostNode(x_ex_left,y_ex_left,x_ex_cent,y_ex_cent,x_ex_right,y_ex_right,x_int,y_int))
        end
    end
    
    %left face, two rows
    for celloop = 3:-1:2
        irow = celloop; %so irow-1 is the row of new nodes
        for j = 4:Jmax-3  %only consider nodes that have two attached faces 
            [x(irow-1,j),y(irow-1,j)] = addGhostNode(x(irow,j-1),y(irow,j-1),x(irow,j),y(irow,j),x(irow,j+1),y(irow,j+1),x(irow+1,j),y(irow+1,j));
            %[gx,gy] = addGhostNode(x_ex_left,y_ex_left,x_ex_cent,y_ex_cent,x_ex_right,y_ex_right,x_int,y_int)
        end
    end
    
    %right face, two rows
    for celloop = 2:-1:1
        irow = Imax-celloop; %so irow+1 is the row of new nodes
        for j = 4:Jmax-3  %only consider nodes that have two attached faces 
            [x(irow+1,j),y(irow+1,j)] = addGhostNode(x(irow,j+1),y(irow,j+1),x(irow,j),y(irow,j),x(irow,j-1),y(irow,j-1),x(irow-1,j),y(irow-1,j));
            %[gx,gy] = addGhostNode(x_ex_left,y_ex_left,x_ex_cent,y_ex_cent,x_ex_right,y_ex_right,x_int,y_int)
        end
    end
    
    %bottom face, two rows
    for celloop = 3:-1:2
        jrow = celloop; %so irow-1 is the row of new nodes
        for i = 4:Imax-3  %only consider nodes that have two attached faces 
            [x(i,jrow-1),y(i,jrow-1)] = addGhostNode(x(i-1,jrow),y(i-1,jrow),x(i,jrow),y(i,jrow),x(i+1,jrow),y(i+1,jrow),x(i,jrow+1),y(i,jrow+1));
            %[gx,gy] = addGhostNode(x_ex_left,y_ex_left,x_ex_cent,y_ex_cent,x_ex_right,y_ex_right,x_int,y_int)
        end
    end   
    
    nodes_x_output = x;
    nodes_y_output = y;
    
    
    %% Plot input and output
%     figure;
%         for k = 1:Imax-4
%         scatter(nodes_x_input(k,:),nodes_y_input(k,:));
%         hold on;
%         end
%     title('Input nodes');
%     grid on;
    
    figure;
        for k = 1:Imax
        scatter(nodes_x_output(k,:),nodes_y_output(k,:));
        hold on;
        end
    title('Output nodes');
    grid on;
end

