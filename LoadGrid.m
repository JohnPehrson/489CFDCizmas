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
        
    else
        fprintf('Please input the correct grid descriptor');
    end
end

