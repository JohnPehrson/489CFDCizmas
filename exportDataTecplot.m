function exportDataTecplot(nodes_x,nodes_y,cells_q,nodes_Imax,nodes_Jmax,cells_Imax,cells_Jmax)
%This function takes the final results of the program and writes it to a
%text file that can be imported into tecPlot to visualize the data

%write x,y,q1,q2,q3,q4 vectors that give the data in 1d
x = NaN(1,nodes_Imax*nodes_Jmax);
y = x;
q1 = NaN(1,cells_Imax*cells_Jmax);
q2 = q1;
q3 = q1;
q4 = q1;

nodeit = 1;
for j = 1:nodes_Jmax
    for i = 1:nodes_Imax
        x(nodeit) = nodes_x(i,j);
        y(nodeit) = nodes_y(i,j);
        nodeit = nodeit+1;
    end
end

cellit = 1;
for j = 1:cells_Jmax
    for i = 1:cells_Imax
        q1(cellit) = cells_q(i,j,1);
        q2(cellit) = cells_q(i,j,2);
        q3(cellit) = cells_q(i,j,3);
        q4(cellit) = cells_q(i,j,4); 
        cellit = cellit+1;
    end
end




%find and open the data file
fileID = fopen('testTecPlot.txt','w');
fprintf(fileID,' VARIABLES = "X", "Y", "q1", "q2", "q3", "q4"\n');
fprintf(fileID,'ZONE I=%2g J=%2g DATAPACKING=BLOCK VARLOCATION=([3,4,5,6]=CELLCENTERED)',cells_Imax,cells_Jmax); %[3,4,5,6] correspond to the number of cell-centered variables, in this case qvec




fclose(fileID);
end

