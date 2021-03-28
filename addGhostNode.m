function [gx,gy] = addGhostNode(x_ex_left,y_ex_left,x_ex_cent,y_ex_cent,x_ex_right,y_ex_right,x_int,y_int)
%This function gives the x and y location of a ghost node
%Requires an input Px,Py (location of interior node guiding the ghost node)
    %called x and y int in the program
%Also requires the exterior node location to find distance
    %called x and y int central
%Needs adjacent nodes for boundary normal calculations
    %x and y int left and right
%Swapsigns would change the sign/direction of the normal

    [n1x,n1y] = cellnormal(x_ex_left,x_ex_cent,y_ex_left,y_ex_cent); %left and cent
    [n2x,n2y] = cellnormal(x_ex_cent,x_ex_right,y_ex_cent,y_ex_right); %right and cent
    %[nx,ny] = cellnormal(x1,x2,y1,y2)    
    
    %average the normal vectors
    nx = (n1x+n2x)/2;
    ny = (n1y+n2y)/2;
    
    %distance
    dx = abs(-x_int+x_ex_cent);
    dy = -y_int+y_ex_cent;
    
    %ghost node location
    gx = x_int+2*dx*nx;
    gy = y_int+2*dy*ny;    
end

