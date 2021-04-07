function exportDataTecplot(user_Mach,iterations,nodes_x,nodes_y,plot_cells_q,plot_cells_f,plot_cells_g,nodes_Imax,nodes_Jmax,cells_Imax,cells_Jmax,plot_it_reduced,plot_it)
%This function takes the final results of the program and writes it to a
%text file that can be imported into tecPlot to visualize the data


%find and open the data file
newmach = erase(num2str(user_Mach),'.');
filetitle = ['Pehrson_P2_' newmach '_' num2str(iterations-1) '.txt'];
fileID = fopen(filetitle,'w');


for timeloop = 1:(1+plot_it_reduced) %for loop for individual sets of data in time/iterations
  
  %For each loop, only consider one q,f,g set of data
  cells_q = plot_cells_q(timeloop,:,:,:);
  cells_f = plot_cells_g(timeloop,:,:,:);
  cells_g = plot_cells_g(timeloop,:,:,:);
    
    %write x,y,q1,q2,q3,q4 vectors that give the data in 1d
    x = NaN(1,nodes_Imax*nodes_Jmax);
    y = x;
    q1 = NaN(1,cells_Imax*cells_Jmax);
    q2 = q1;
    q3 = q1;
    q4 = q1;
    pstatic = q1;
    Mach = q1;

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
                q1(cellit) =  plot_cells_q(timeloop,i,j,1);
                q2(cellit) = plot_cells_q(timeloop,i,j,2);
                q3(cellit) = plot_cells_q(timeloop,i,j,3);
                q4(cellit) = plot_cells_q(timeloop,i,j,4);
                pstatic(cellit) = plot_cells_f(timeloop,i,j,2)-((plot_cells_g(timeloop,i,j,2))^2)/plot_cells_q(timeloop,i,j,1);   
                Mach(cellit) = sqrt((q2(cellit)^2+q3(cellit)^2)/(q1(cellit)^2));
                cellit = cellit+1;
            end
        end

        %Set time variable
        IterationsCount = plot_it*(timeloop-1);
        IterationsName = num2str(IterationsCount);


        %write to the data file
        fprintf(fileID,' VARIABLES = "X", "Y", "q1", "q2", "q3", "q4","Pressure","Mach"\n');
        fprintf(fileID,'ZONE T="%2s Iterations", I=%2g, J=%2g, DATAPACKING=BLOCK VARLOCATION=([3,4,5,6,7,8]=CELLCENTERED)\n',IterationsName,nodes_Imax,nodes_Jmax); %[3,4,5,6] correspond to the number of cell-centered variables, in this case qvec
        for i = 1:nodes_Imax*nodes_Jmax
            fprintf(fileID,'%.5g ', x(i));
        end
        fprintf(fileID,'\n');
        for i = 1:nodes_Imax*nodes_Jmax
            fprintf(fileID,'%.5g ', y(i));
        end
        fprintf(fileID,'\n');
        for i = 1:cells_Imax*cells_Jmax
            fprintf(fileID,'%.5g ', q1(i));
        end
        fprintf(fileID,'\n');
        for i = 1:cells_Imax*cells_Jmax
            fprintf(fileID,'%.5g ', q2(i));
        end
        fprintf(fileID,'\n');
        for i = 1:cells_Imax*cells_Jmax
            fprintf(fileID,'%.5g ', q3(i));
        end
        fprintf(fileID,'\n');
        for i = 1:cells_Imax*cells_Jmax
            fprintf(fileID,'%.5g ', q4(i));
        end
        fprintf(fileID,'\n');
        for i = 1:cells_Imax*cells_Jmax
            fprintf(fileID,'%.5g ', pstatic(i));
        end
        fprintf(fileID,'\n');
        for i = 1:cells_Imax*cells_Jmax
            fprintf(fileID,'%.5g ', Mach(i));
        end
        fprintf(fileID,'\n');

end











    fclose(fileID);
end

