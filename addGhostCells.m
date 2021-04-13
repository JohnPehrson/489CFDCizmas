function [nodes_x_output,nodes_y_output] = addGhostCells(nodes_x_input,nodes_y_input)
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
   
    %find dx and dy, used everywhere
    dx = nodes_x_input(2,Jmax-4)-nodes_x_input(1,Jmax-4);
    dy = nodes_y_input(1,Jmax-4)-nodes_y_input(1,Jmax-5);
    
    %go wide first, add ghost nodes in dx direction
        %left side
    nodes_x_output(1:2,(2):(Jmax-1)) = [nodes_x_output(3,(2):(Jmax-1))-2*dx;nodes_x_output(3,(2):(Jmax-1))-dx];
    nodes_y_output(1:2,(2):(Jmax-1)) = [nodes_y_output(3,(2):(Jmax-1));nodes_y_output(3,(2):(Jmax-1))];
        %right side
    nodes_x_output((Imax-1):Imax,(2):(Jmax-1)) = [nodes_x_output(Imax-2,(2):(Jmax-1))+dx;nodes_x_output(Imax-2,(2):(Jmax-1))+2*dx];
    nodes_y_output((Imax-1):Imax,(2):(Jmax-1)) = [nodes_y_output(Imax-2,(2):(Jmax-1));nodes_y_output(Imax-2,(2):(Jmax-1))];
    
    %top and bottom, all the way along
        %top
    nodes_x_output(1:Imax,(Jmax-1):Jmax) = [nodes_x_output(1:Imax,Jmax-2),nodes_x_output(1:Imax,Jmax-2);];
    nodes_y_output(1:Imax,(Jmax-1):Jmax) = [nodes_y_output(1:Imax,Jmax-2)+dy,nodes_y_output(1:Imax,Jmax-2)+2*dy];
        %bottom
    nodes_x_output(1:Imax,1:2) = [nodes_x_output(1:Imax,3),nodes_x_output(1:Imax,3);];
    nodes_y_output(1:Imax,1:2) = [nodes_y_output(1:Imax,3)-2*dy,nodes_y_output(1:Imax,3)-dy];
        
        
%% Plot input and output
%     figure;
%         for k = 1:Imax-4
%         scatter(nodes_x_input(k,:),nodes_y_input(k,:));
%         hold on;
%         end
%     title('Input nodes');
%     grid on;
%     
%     figure;
%         for k = 1:Imax
%         scatter(nodes_x_output(k,:),nodes_y_output(k,:));
%         hold on;
%         end
%     title('Output nodes');
%     grid on;
%     xlim([-1,6]);
%     ylim([-.5,1.5]);
end

