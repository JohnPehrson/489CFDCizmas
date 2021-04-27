function [xout,yout] = LoadGrid(gridtype)
%Loads the laplace grid based on the input type defined by the user

    if strcmp(gridtype,'coarse') %coarsegrid
        load('coarselaplace_gridX.mat');
        load('coarselaplace_gridY.mat');
    elseif strcmp(gridtype,'medium') %medium grid
        load('mediumlaplace_gridX.mat');
        load('mediumlaplace_gridY.mat');
    elseif strcmp(gridtype,'fine')  %fine grid
        load('finelaplace_gridX.mat');
        load('finelaplace_gridY.mat');
    elseif strcmp(gridtype,'testalg') %test grid
        load('testgridX.mat');
        load('testgridY.mat');
        xout = x;
        yout = y;
    else
        fprintf('Please input the correct grid descriptor');
    end
end

