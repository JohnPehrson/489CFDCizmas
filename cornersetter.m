function [cells_q,cells_f,cells_g] = cornersetter(cells_q,cells_f,cells_g,cells_Imax,cells_Jmax)
%The goal of this function is to change the q,f,g in the corners of the
%ghost boundary. This is a superficial change that is meant to change the
%VISUALS of the grid, which will make later data visualization better

%q
cells_q(1:2,1:2,:) = cells_q(1:2,3:4,:);
cells_q((cells_Imax-1):cells_Imax,1:2,:) = cells_q((cells_Imax-1):cells_Imax,3:4,:);
cells_q((cells_Imax-1):cells_Imax,(cells_Jmax-1):cells_Jmax,:) = cells_q((cells_Imax-1):cells_Imax,(cells_Jmax-3):(cells_Jmax-2),:);
cells_q(1:2,(cells_Jmax-1):cells_Jmax,:) = cells_q(1:2,(cells_Jmax-3):(cells_Jmax-2),:);

%f
cells_f(1:2,1:2,:) = cells_f(1:2,3:4,:);
cells_f((cells_Imax-1):cells_Imax,1:2,:) = cells_f((cells_Imax-1):cells_Imax,3:4,:);
cells_f((cells_Imax-1):cells_Imax,(cells_Jmax-1):cells_Jmax,:) = cells_f((cells_Imax-1):cells_Imax,(cells_Jmax-3):(cells_Jmax-2),:);
cells_f(1:2,(cells_Jmax-1):cells_Jmax,:) = cells_f(1:2,(cells_Jmax-3):(cells_Jmax-2),:);

%g
cells_g(1:2,1:2,:) = cells_g(1:2,3:4,:);
cells_g((cells_Imax-1):cells_Imax,1:2,:) = cells_g((cells_Imax-1):cells_Imax,3:4,:);
cells_g((cells_Imax-1):cells_Imax,(cells_Jmax-1):cells_Jmax,:) = cells_g((cells_Imax-1):cells_Imax,(cells_Jmax-3):(cells_Jmax-2),:);
cells_g(1:2,(cells_Jmax-1):cells_Jmax,:) = cells_g(1:2,(cells_Jmax-3):(cells_Jmax-2),:);


end

