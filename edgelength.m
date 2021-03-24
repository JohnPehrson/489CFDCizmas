function [length] = edgelength(i,j,x,y,face)
%This function finds the length of the edge of a cell. There are 4 possible
%edges to each cell, denoted by 'face' = 'N','E','S','W'
%i and j refer to cells, not nodes
%x and y are currently the whole matrix of node values

[X_abcd,Y_abcd] = nodes_touch_cell(i,j,x,y);
    switch face
        case 'N'  %use c and d
            length = sqrt((X_abcd(3)-X_abcd(4))^2+(Y_abcd(3)-Y_abcd(4))^2);
        case 'E' %use c and b
            length = sqrt((X_abcd(3)-X_abcd(2))^2+(Y_abcd(3)-Y_abcd(2))^2);
        case 'S' %use a and b
            length = sqrt((X_abcd(2)-X_abcd(1))^2+(Y_abcd(2)-Y_abcd(1))^2);
        case 'W' %use d and a
            length = sqrt((X_abcd(4)-X_abcd(1))^2+(Y_abcd(4)-Y_abcd(1))^2);
    end
end

