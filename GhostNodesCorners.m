function [nodes_x_output,nodes_y_output] = GhostNodesCorners(nodes_x_output,nodes_y_output,Imax,Jmax)
%A subroutine that will add the 16 special corner cases of the ghost nodes 

%top left-i(left-right)
dx = nodes_x_output(4,Jmax-2)-nodes_x_output(3,Jmax-2);
dy = nodes_y_output(4,Jmax-2)-nodes_y_output(3,Jmax-2);
nodes_x_output(1:2,Jmax-2) = [nodes_x_output(3,Jmax-2)-2*dx,nodes_x_output(3,Jmax-2)-dx];
nodes_y_output(1:2,Jmax-2) = [nodes_y_output(3,Jmax-2)-2*dy,nodes_y_output(3,Jmax-2)-dy];

%top left-j(up-down)
dx = nodes_x_output(3,Jmax-2)-nodes_x_output(3,Jmax-3);
dy = nodes_y_output(3,Jmax-2)-nodes_y_output(3,Jmax-3);
nodes_x_output(3,(Jmax-1):Jmax) = [nodes_x_output(3,Jmax-2)+dx,nodes_x_output(3,Jmax-2)+2*dx];
nodes_y_output(3,(Jmax-1):Jmax) = [nodes_y_output(3,Jmax-2)+dy,nodes_y_output(3,Jmax-2)+2*dy];

%top right-i
dx = nodes_x_output(Imax-2,Jmax-2)-nodes_x_output(Imax-3,Jmax-2);
dy = nodes_y_output(Imax-2,Jmax-2)-nodes_y_output(Imax-3,Jmax-2);
nodes_x_output((Imax-1):Imax,Jmax-2) = [nodes_x_output(Imax-2,Jmax-2)+dx,nodes_x_output(Imax-2,Jmax-2)+2*dx];
nodes_y_output((Imax-1):Imax,Jmax-2) = [nodes_y_output(Imax-2,Jmax-2)+dy,nodes_y_output(Imax-2,Jmax-2)+2*dy];

%top right-j
dx = nodes_x_output(Imax-2,Jmax-2)-nodes_x_output(Imax-2,Jmax-3);
dy = nodes_y_output(Imax-2,Jmax-2)-nodes_y_output(Imax-2,Jmax-3);
nodes_x_output(Imax-2,(Jmax-1):Jmax) = [nodes_x_output(Imax-2,Jmax-2)+dx,nodes_x_output(Imax-2,Jmax-2)+2*dx];
nodes_y_output(Imax-2,(Jmax-1):Jmax) = [nodes_y_output(Imax-2,Jmax-2)+dy,nodes_y_output(Imax-2,Jmax-2)+2*dy];

%bottom left-i
dx = nodes_x_output(4,3)-nodes_x_output(3,3);
dy = nodes_y_output(4,3)-nodes_y_output(3,3);
nodes_x_output(1:2,3) = [nodes_x_output(3,3)-2*dx,nodes_x_output(3,3)-dx];
nodes_y_output(1:2,3) = [nodes_y_output(3,3)-2*dy,nodes_y_output(3,3)-dy];

%bottom left -j
dx = nodes_x_output(3,4)-nodes_x_output(3,3);
dy = nodes_y_output(3,4)-nodes_y_output(3,3);
nodes_x_output(3,1:2) = [nodes_x_output(3,3)-2*dx,nodes_x_output(3,3)-dx];
nodes_y_output(3,1:2) = [nodes_y_output(3,3)-2*dy,nodes_y_output(3,3)-dy];

%bottom right-i
dx = nodes_x_output(Imax-2,3)-nodes_x_output(Imax-3,3);
dy = nodes_y_output(Imax-2,3)-nodes_y_output(Imax-3,3);
nodes_x_output((Imax-1:Imax),3) = [nodes_x_output(Imax-2,3)+dx,nodes_x_output(Imax-2,3)+2*dx];
nodes_y_output((Imax-1:Imax),3) = [nodes_y_output(Imax-2,3)+dy,nodes_y_output(Imax-2,3)+2*dy];

%bottom right-j
dx = nodes_x_output(Imax-2,4)-nodes_x_output(Imax-2,3);
dy = nodes_y_output(Imax-2,4)-nodes_y_output(Imax-2,3);
nodes_x_output(Imax-2,1:2) = [nodes_x_output(Imax-2,3)-2*dx,nodes_x_output(Imax-2,3)-dx];
nodes_y_output(Imax-2,1:2) = [nodes_y_output(Imax-2,3)-2*dy,nodes_y_output(Imax-2,3)-dy];

end

