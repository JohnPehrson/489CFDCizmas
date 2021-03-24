function [nearbynodesX_abcd,nearbynodesY_abcd] = nodes_touch_cell(cells_i,cells_j,nodes_x,nodes_y)
%Reports the x and y location of nodes around a cell given a cell location
%and the x and y matraxies of all node locations

nearbynodesX_abcd = [nodes_x(cells_i,cells_j),nodes_x(cells_i+1,cells_j),nodes_x(cells_i+1,cells_j+1),nodes_x(cells_i,cells_j+1)];
nearbynodesY_abcd = [nodes_y(cells_i,cells_j),nodes_y(cells_i+1,cells_j),nodes_y(cells_i+1,cells_j+1),nodes_y(cells_i,cells_j+1)];
end
