function [xout,yout] = LoadGrid(gridtype)
%Loads the laplace grid based on the input type defined by the user

    if strcmp(gridtype,'coarse') %coarsegrid
        load('coarsegridX.mat');
        load('coarsegridY.mat');
    elseif strcmp(gridtype,'medium') %medium grid
        load('mediumgridX.mat');
        load('mediumgridY.mat');
    elseif strcmp(gridtype,'fine')
        %fine grid
        
    elseif strcmp(gridtype,'testalg') %test grid
        load('51_10_alg_gridX.mat');
        load('51_10_alg_gridY.mat');
        xout = x;
        yout = y;
   elseif strcmp(gridtype,'testlaplace') %test grid
        load('51_10_laplace_gridX.mat');
        load('51_10_laplace_gridY.mat');
    else
        fprintf('Please input the correct grid descriptor');
    end
end

