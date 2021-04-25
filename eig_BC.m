function [cells_eig] = eig_BC(cells_Imax,cells_Jmax,cells_eig)
%This function defines the eigenvalues associated with the boundary
%conditions in order to specify dissipation on specific boundaries

%bottom wall
for i = 3:(cells_Imax-2)
    cells_eig(i,2,1) = cells_eig(i,3,3);
    cells_eig(i,2,3) = cells_eig(i,3,1);
end
%top wall
jtop = cells_Jmax-2;
for i = 3:(cells_Imax-2)
    cells_eig(i,jtop+1,3) = cells_eig(i,jtop,1);
    cells_eig(i,jtop+1,1) = cells_eig(i,jtop,3);
end
%inlet
for j = 3:(cells_Jmax-2)
    cells_eig(2,j,2) = cells_eig(3,j,4);
    cells_eig(2,j,4) = cells_eig(3,j,2);
end
%outlet
iout = cells_Imax-2;
for j = 3:(cells_Jmax-2)
    cells_eig(iout+1,j,4) = cells_eig(iout,j,2);
    cells_eig(iout+1,j,2) = cells_eig(iout,j,4);
end

end

