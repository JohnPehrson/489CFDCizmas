function [celltype] = cellidentifier(cells_Imax,cells_Jmax)
%Identifies the type of faces that the cell has. Ex: The cell has the
%bottom face on the wall-boundary.

    %There are 9 types of cells in this grid. This program will identify them
    %as such. The label will be used for later computations to ensure that no
    %requests are made for cells that don't exist

    %1: Interior cell bounded on all sides by other cells
    %2: Cell on the bottom wall but otherwise connected to other cells
    %3: Cell on the top wall but otherwise connected to other cells
    %4: Cell on the left face/inlet that is otherwise connected to other cells
    %5: Cell on the right face/outlet that is otherwise connected to other
    %cells
    %6: Bottom left corner cell, connected to the bottom wall and the inlet
    %face
    %7: Bottom right corner cell, connected to the bottom wall and the outlet
    %face
    %8: Top right corner cell, connected to the top wall and the outlet face
    %9: Top left corner cell, connected to the top wall and the inlet face

    %initialize the returned matrix
    celltype = zeros(cells_Imax,cells_Jmax);

    %loop through the cells to identify types
    for i = 1:cells_Imax
        for j = 1:cells_Jmax
            
            if i>=2 && i<=(cells_Imax-1) && j>=2 && j<=(cells_Jmax-1) %central cell
                celltype(i,j) = 1;
            elseif (i ==1 && j==1) %bottom left corner
                celltype(i,j) = 6;
            elseif (i == 1 && j == cells_Jmax) %bottom right corner
                celltype(i,j) = 7;
            elseif (i == cells_Imax && j == 1) %top left corner
                celltype(i,j) = 9;
            elseif (i == cells_Imax && j == cells_Jmax) %top right corner
                celltype(i,j) = 8;
            elseif ((i>=2 && i<=(cells_Imax-1)) && j == 1) %bottom face not corners
                celltype(i,j) = 2;
            elseif ((i>=2 && i<=(cells_Imax-1)) && j == cells_Jmax) %top face not corners
                celltype(i,j) = 3;                
            elseif (i == 1 && (j>=2 && j<=(cells_Jmax-1))) %left face not corners
                celltype(i,j) = 4;
            elseif (i == cells_Imax && (j>=2 && j<=(cells_Jmax-1))) %right face not corners
                celltype(i,j) = 5; 
        end
    end
    





end

