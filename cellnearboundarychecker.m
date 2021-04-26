function [ind1,ind2] = cellnearboundarychecker(cells_i,cells_j,cells_Imax,cells_Jmax)
%Identifies if the cell in question is bordering any ghost cells. Used by
%the RK time loop to identify what boundary conditions to apply to nearby
%cells
%Assume 2 layers of ghost cells

imin = 3;
jmin = 3;
imax = cells_Imax-2;
jmax = cells_Jmax-2;

%check E/W using i
    if cells_i == imin
        ind1 = 'W';
    elseif cells_i == imax
        ind1 = 'E';
    else
        ind1 = 'neither';
    end

%check N/S using j
    if cells_j == jmin
        ind2 = 'S';
    elseif cells_j == jmax
        ind2 = 'N';
    else
        ind2 = 'neither';
    end
end

