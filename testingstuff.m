clear all;close all;clc;

% load('testgridX.mat');
% load('testgridY.mat');
% nodes_Imax = length(xout(:,1));
% nodes_Jmax = length(xout(1,:));
% cells_Imax = nodes_Imax-1;
% cells_Jmax = nodes_Jmax-1;
% cells_p = zeros(cells_Imax,cells_Jmax,5);
% cells_f = cells_p;
% cells_g = cells_p;
% 
% for k = 1:nodes_Imax
% scatter(xout(k,:),yout(k,:))
% hold on;
% end
% 
%  Startpoint = [0 0];
%  Endpoint = [-0.1 1];
%  direction = Endpoint - Startpoint;
%  direction_scaled = direction/norm(direction) % vector with length 1
% %  up =  [-direction_scaled(2),direction_scaled(1)]
% % right = [direction_scaled(2),-direction_scaled(1)]
% % down = [direction_scaled(2),-direction_scaled(1)]
% left = [-direction_scaled(2),direction_scaled(1)]

% vector = [1;2;3;4];
% testmatrixstuff = 5*11*0.1*vector

velocityvector = [1;1];
walltan = [1;0];
angle = acosd(dot(velocityvector/norm(velocityvector),walltan/norm(walltan)));
if (velocityvector(2)-walltan(2))<0
    angle = angle*-1;
end
vnew = [cosd(2*angle),sind(2*angle);-sind(2*angle),cosd(2*angle)]*velocityvector

