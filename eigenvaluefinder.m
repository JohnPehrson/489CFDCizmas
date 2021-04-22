function [eigout] = eigenvaluefinder(x_abcd,y_abcd,q,c)
%This function calculates the absolute value of the largest eigenvalue of
%the euler equation aligned with the xi/n directions on the N,E,S,W faces
%c is speed of sound
%x_abcd,y_abcd are the node locations of the 4 adjacent nodes to the cell
%q is a 1x4? matrix 
%I need to report the eigenvalue on each of the 4 faces


%calculate information about the cell cell
q = squeeze(q);
u = q(2)/q(1);
v = q(3)/q(1);
V = [u;v];

%% Split by face

    %N
    [nx,ny] = cellnormal(x_abcd(4),x_abcd(3),y_abcd(4),y_abcd(3));
    n_N = [nx;ny];
    un_N = dot(V,n_N);
    eig_N = abs(un_N)+c;

    %E
    [nx,ny] = cellnormal(x_abcd(3),x_abcd(2),y_abcd(3),y_abcd(2));
    n_E = [nx;ny];
    un_E = dot(V,n_E);
    eig_E = abs(un_E)+c;
    
    %S
    [nx,ny] = cellnormal(x_abcd(2),x_abcd(1),y_abcd(2),y_abcd(1));
    n_S = [nx;ny];
    un_S = dot(V,n_S);
    eig_S = abs(un_S)+c;

    %W
    [nx,ny] = cellnormal(x_abcd(1),x_abcd(4),y_abcd(1),y_abcd(4));
    n_W = [nx;ny];
    un_W = dot(V,n_W);
    eig_W = abs(un_W)+c;
    
%Put it into the right output format for the overall matrix
eigout = zeros(1,1,4);
eigout(1,1,:) = [eig_N,eig_E,eig_S,eig_W];

end

