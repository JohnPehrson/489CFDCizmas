function [cellfluxesf_NESW,cellfluxesg_NESW] = cell_fluxes(i,j,f,g,celltype)
%This function calculates the fluxes that pass through the edges of a given
%cell. This function outputs vectors that correspond to different faces of
%the defined cell. Used in the calculation of the Residuals.

%The calculation of fluxes is dependent on the type of cell (location of
%the cell in the grid). Celltype is a constant that specifies location. 

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

    switch celltype %what celltype?
        case 1 % interior cell
            fN = 0.5*(f(i,j)+f(i,j+1));
            fE = 0.5*(f(i,j)+f(i+1,j));
            fS = 0.5*(f(i,j)+f(i,j-1));
            fW = 0.5*(f(i,j)+f(i-1,j));
            cellfluxesf_NESW = 0.5*[fN,fE,fS,fW];
            gN = 0.5*(g(i,j)+g(i,j+1));
            gE = 0.5*(g(i,j)+g(i+1,j));
            gS = 0.5*(g(i,j)+g(i,j-1));
            gW = 0.5*(g(i,j)+g(i-1,j));
            cellfluxesg_NESW = 0.5*[gN,gE,gS,gW];
        case 2 %bottom wall cell, no flux through the south
            fN = 0.5*(f(i,j)+f(i,j+1));
            fE = 0.5*(f(i,j)+f(i+1,j));
            fS = 0;
            fW = 0.5*(f(i,j)+f(i-1,j));
            cellfluxesf_NESW = 0.5*[fN,fE,fS,fW];
            gN = 0.5*(g(i,j)+g(i,j+1));
            gE = 0.5*(g(i,j)+g(i+1,j));
            gS = 0;
            gW = 0.5*(g(i,j)+g(i-1,j));
            cellfluxesg_NESW = 0.5*[gN,gE,gS,gW];
        case 3 %top wall cell, no flux through the north
            fN = 0;
            fE = 0.5*(f(i,j)+f(i+1,j));
            fS = 0.5*(f(i,j)+f(i,j-1));
            fW = 0.5*(f(i,j)+f(i-1,j));
            cellfluxesf_NESW = 0.5*[fN,fE,fS,fW];
            gN = 0;
            gE = 0.5*(g(i,j)+g(i+1,j));
            gS = 0.5*(g(i,j)+g(i,j-1));
            gW = 0.5*(g(i,j)+g(i-1,j));
            cellfluxesg_NESW = 0.5*[gN,gE,gS,gW];
        case 4 %inlet cell, 
            fN = 0.5*(f(i,j)+f(i,j+1));
            fE = 0.5*(f(i,j)+f(i+1,j));
            fS = 0.5*(f(i,j)+f(i,j-1));
            fW = 0.5*(f(i,j)+f(i-1,j));
            cellfluxesf_NESW = 0.5*[fN,fE,fS,fW];
            gN = 0.5*(g(i,j)+g(i,j+1));
            gE = 0.5*(g(i,j)+g(i+1,j));
            gS = 0.5*(g(i,j)+g(i,j-1));
            gW = 0.5*(g(i,j)+g(i-1,j));
            cellfluxesg_NESW = 0.5*[gN,gE,gS,gW];
        case 5 %outlet cell, 
            
        case 6 %inlet & bottom wall cell, 
            
        case 7 %outlet & bottom wall cell, 
            
        case 8 %outlet & top wall cell, 
            
        case 9 %inlet & top wall cell, 
            
    end


end

