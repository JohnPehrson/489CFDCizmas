function [cells_range_i] = wallbumpcellfinder(nodes_y,nodes_Imax)
%This function identifies which cells are along the bottom wall bump. this
%is used to identify which cells should be used for wall-force calculations

%assume 2 rounds of ghost cells, cells have 1 fewer index than nodes do
j = 3; %bottom wall interior
i_checker = NaN(1,nodes_Imax-4);
    for i = 3:nodes_Imax-2
        if nodes_y(i,j)>0
            i_checker(i) = i;
        else
            %nothing
        end
    end
[null,Imin] = min(i_checker);
[null,Imax] = max(i_checker);

cells_range_i =linspace(Imin-1,Imax,Imax-Imin+2);
end

