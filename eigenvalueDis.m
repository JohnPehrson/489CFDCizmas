function [eigenvalueout] = eigenvalueDis(c,x_abcd,y_abcd,q3by3,face)
%This function calculates the absolute value of the largest eigenvalue of
%the euler equation aligned with the xi/n directions. 
%c is speed of sound
%x_abcd,y_abcd are the node locations of the 4 adjacent nodes to the cell
%q is a 3x3x4 matrix that gives the q values of the central cell and the 4
    %bordering cells, used for velocity calculations
%face determines what face we are looking at: 'N','E','S','W'

q = q3by3;
%calculate information about the central cell
ucent = q(2,2,2)/q(2,2,1);
vcent = q(2,2,3)/q(2,2,1);
velo_central = [ucent;vcent];

%split into different cases based on which edge is being used
    switch face
        case 'N'
            [nx,ny] = cellnormal(x_abcd(4),x_abcd(3),y_abcd(4),y_abcd(3));
            normal = [nx;ny];
            uadj = q(1,2,2)/q(1,2,1);
            vadj = q(1,2,3)/q(1,2,1);
            velo_adj = [uadj;vadj];
        case 'E'
            [nx,ny] = cellnormal(x_abcd(3),x_abcd(2),y_abcd(3),y_abcd(2));
            normal = [nx;ny];
            uadj = q(2,3,2)/q(2,3,1);
            vadj = q(2,3,3)/q(2,3,1);
            velo_adj = [uadj;vadj];
        case 'S'
            [nx,ny] = cellnormal(x_abcd(2),x_abcd(1),y_abcd(2),y_abcd(1));
            normal = [nx;ny];
            uadj = q(3,2,2)/q(3,2,1);
            vadj = q(3,2,3)/q(3,2,1);
            velo_adj = [uadj;vadj];
        case 'W'
            [nx,ny] = cellnormal(x_abcd(1),x_abcd(4),y_abcd(1),y_abcd(4));
            normal = [nx;ny];
            uadj = q(2,1,2)/q(2,1,1);
            vadj = q(2,1,3)/q(2,1,1);
            velo_adj = [uadj;vadj];
    end

eigenvalueout = c+0.5*(abs(dot(velo_central,normal))+abs(dot(velo_adj,normal)));
end

