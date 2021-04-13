function [eigenvalueout_z,eigenvalueout_n] = eigenvaluefinder(x_abcd,y_abcd,q,c)
%This function calculates the absolute value of the largest eigenvalue of
%the euler equation aligned with the xi/n directions. 
%c is speed of sound
%x_abcd,y_abcd are the node locations of the 4 adjacent nodes to the cell
%q is a 1x4? matrix 
%p is a single pressure
%neighboring cells, used for speed of sound calcs
%face determines what face we are looking at: 'N','E','S','W'


%calculate information about the cell cell
u = q(2)/q(1);
v = q(3)/q(1);
V = [u;v];
%average the normal of the two faces in each direciton
[normal_zrx,normal_zry] = cellnormal(x_abcd(2),x_abcd(3),y_abcd(2),y_abcd(3));
[normal_zlx,normal_zly] = cellnormal(x_abcd(1),x_abcd(4),y_abcd(1),y_abcd(4));
normal_z = [0.5.*(normal_zrx+normal_zlx),0.5.*(normal_zry+normal_zly)];

[normal_nux,normal_nuy] = cellnormal(x_abcd(3),x_abcd(4),y_abcd(3),y_abcd(4));
[normal_nbx,normal_nby] = cellnormal(x_abcd(2),x_abcd(1),y_abcd(2),y_abcd(1));
normal_n = [0.5.*(normal_nux+normal_nbx),0.5.*(normal_nuy+normal_nby)];


eigenvalueout_z = c+abs(dot(V,normal_z));
eigenvalueout_n = c+abs(dot(V,normal_n));
end

