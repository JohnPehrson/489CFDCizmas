function [R] = findResiduals(nodes_x,nodes_y,cells_f,cells_g,cells_Imax,cells_Jmax,celltype)
%Calculates the residuals (flux integral) for the whole grid. Requires node
%locations and cell boundary fluxes. Note that the residual is the output,
%and is recalculated throughout the solving process


R = NaN(cells_Imax,cells_Jmax);
%loop through each cell in the grid
for i = 1:cells_Imax
    for j = 1:cells_Jmax
        [x,y] = nodes_touch_cell(i,j,nodes_x,nodes_y); %abcd
        [f,g] = cell_fluxes(i,j,cells_f,cells_g,celltype(i,j)); %NESW
        
        R(i,j) = f(1)*(y(4)-y(3))-g(1)*(x(4)-x(3))+f(4)*(y(1)-y(4))-g(4)*(x(1)-x(4))+...
                 f(3)*(y(2)-y(1))-g(3)*(x(2)-x(1))+f(2)*(y(3)-y(2))-g(2)*(x(3)-x(2));
    
    end
end



end

