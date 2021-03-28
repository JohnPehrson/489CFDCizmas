function [nx,ny] = cellnormal(x1,x2,y1,y2)
%Finds the normal vector to a line given the endpoints of the line
%Used to find the normal vector to a cell boundary face

%**flip the signs of the outputs nx and ny to suit your needs**%

Tx = x2-x1;
Ty = y2-y1;
Tmag = sqrt(Tx^2+Ty^2);
nx = -Ty/Tmag;
ny = Tx/Tmag;
end

