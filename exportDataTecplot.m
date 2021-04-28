function exportDataTecplot(user_Mach,user_MeshQual,iterations,nodes_x,nodes_y,plot_cells_q,plot_cells_f,plot_cells_g,plot_Residual,plot_cells_pressure,plot_cells_c,plot_cells_dissipation,nodes_Imax,nodes_Jmax,cells_Imax,cells_Jmax,user_itmax,report_freq,plot_full)
%This function takes the final results of the program and writes it to a
%text file that can be imported into tecPlot to visualize the data


%find and open the data file
newmach = erase(num2str(user_Mach),'.');
filetitle = ['Pehrson_P2_' newmach '_' user_MeshQual '_' num2str(iterations-1) '.txt'];
fileID = fopen(filetitle,'w');

if plot_full == 1
    for timeloop = 1:report_freq:(1+user_itmax) %for loop for individual sets of data in time/iterations


        %write x,y,q1,q2,q3,q4 vectors that give the data in 1d
        x = NaN(1,nodes_Imax*nodes_Jmax);
        y = x;
        q1 = NaN(1,cells_Imax*cells_Jmax);
        q2 = q1;
        q3 = q1;
        q4 = q1;
        pstatic = q1;
        Mach = q1;
        c = q1;
        f1 = q1;
        f2 = q1;
        f3 = q1;
        f4 = q1;
        g1 = q1;
        g2 = q1;
        g3 = q1;
        g4 = q1;
        r1 = q1;
        r2 = q1;
        r3 = q1;
        r4 = q1;
        d1 = q1;
        d2 = q2;
        d3 = q3;
        d4 = q4;

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
                    f1(cellit) =  plot_cells_f(timeloop,i,j,1);
                    f2(cellit) = plot_cells_f(timeloop,i,j,2);
                    f3(cellit) = plot_cells_f(timeloop,i,j,3);
                    f4(cellit) = plot_cells_f(timeloop,i,j,4);
                    g1(cellit) =  plot_cells_g(timeloop,i,j,1);
                    g2(cellit) = plot_cells_g(timeloop,i,j,2);
                    g3(cellit) = plot_cells_g(timeloop,i,j,3);
                    g4(cellit) = plot_cells_g(timeloop,i,j,4);
                    r1(cellit) =  plot_Residual(timeloop,i,j,1);
                    r2(cellit) = plot_Residual(timeloop,i,j,2);
                    r3(cellit) = plot_Residual(timeloop,i,j,3);
                    r4(cellit) = plot_Residual(timeloop,i,j,4);
                    d1(cellit) = plot_cells_dissipation(timeloop,i,j,1);
                    d2(cellit) = plot_cells_dissipation(timeloop,i,j,2);
                    d3(cellit) = plot_cells_dissipation(timeloop,i,j,3);
                    d4(cellit) = plot_cells_dissipation(timeloop,i,j,4);
                    pstatic(cellit) =  plot_cells_pressure(timeloop,i,j); 
                    Mach(cellit) = sqrt((q2(cellit)^2+q3(cellit)^2)/(q1(cellit)^2));
                    c(cellit) = plot_cells_c(timeloop,i,j);
                    cellit = cellit+1;
                end
            end

            %Set time variable
            IterationsCount = (timeloop-1);
            IterationsName = num2str(IterationsCount);


            %write to the data file
            fprintf(fileID,' VARIABLES = "X", "Y", "q1", "q2", "q3", "q4","f1", "f2", "f3", "f4","g1","g2","g3","g4","r1","r2","r3","r4","d1","d2","d3","d4","Pressure","Mach","c" \n');
            fprintf(fileID,'ZONE T="%2s Iterations", I=%2g, J=%2g, DATAPACKING=BLOCK VARLOCATION=([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]=CELLCENTERED)\n',IterationsName,nodes_Imax,nodes_Jmax); %[3,4,5,6] correspond to the number of cell-centered variables, in this case qvec
            for i = 1:nodes_Imax*nodes_Jmax
                fprintf(fileID,'%.4g ', x(i));
            end
            fprintf(fileID,'\n');
            for i = 1:nodes_Imax*nodes_Jmax
                fprintf(fileID,'%.4g ', y(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', q1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', q2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', q3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', q4(i));
            end
                    fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', f1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', f2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', f3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', f4(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', g1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', g2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', g3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', g4(i));
            end
                    fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', r1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', r2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', r3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', r4(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', d1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', d2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', d3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', d4(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', pstatic(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', Mach(i));
            end
            fprintf(fileID,'\n');
            for i = 1:cells_Imax*cells_Jmax
                fprintf(fileID,'%.4g ', c(i));
            end
            fprintf(fileID,'\n');

    end
else %only plot the interior cells
    
    for timeloop = 1:report_freq:(1+user_itmax) %for loop for individual sets of data in time/iterations

        %write x,y,q1,q2,q3,q4 vectors that give the data in 1d
        x = NaN(1,(nodes_Imax-4)*(nodes_Jmax-4));
        y = x;
        q1 = NaN(1,(cells_Imax-4)*(cells_Jmax-4));
        q2 = q1;
        q3 = q1;
        q4 = q1;
        pstatic = q1;
        Mach = q1;
        c = q1;
        f1 = q1;
        f2 = q1;
        f3 = q1;
        f4 = q1;
        g1 = q1;
        g2 = q1;
        g3 = q1;
        g4 = q1;
        r1 = q1;
        r2 = q1;
        r3 = q1;
        r4 = q1;
        d1 = q1;
        d2 = q2;
        d3 = q3;
        d4 = q4;

        nodeit = 1;
            for j = 3:nodes_Jmax-2
                for i = 3:nodes_Imax-2
                    x(nodeit) = nodes_x(i,j);
                    y(nodeit) = nodes_y(i,j);
                    nodeit = nodeit+1;
                end
            end

            cellit = 1;
            for j = 3:cells_Jmax-2
                for i = 3:cells_Imax-2
                    q1(cellit) =  plot_cells_q(timeloop,i,j,1);
                    q2(cellit) = plot_cells_q(timeloop,i,j,2);
                    q3(cellit) = plot_cells_q(timeloop,i,j,3);
                    q4(cellit) = plot_cells_q(timeloop,i,j,4);
                    f1(cellit) =  plot_cells_f(timeloop,i,j,1);
                    f2(cellit) = plot_cells_f(timeloop,i,j,2);
                    f3(cellit) = plot_cells_f(timeloop,i,j,3);
                    f4(cellit) = plot_cells_f(timeloop,i,j,4);
                    g1(cellit) =  plot_cells_g(timeloop,i,j,1);
                    g2(cellit) = plot_cells_g(timeloop,i,j,2);
                    g3(cellit) = plot_cells_g(timeloop,i,j,3);
                    g4(cellit) = plot_cells_g(timeloop,i,j,4);
                    r1(cellit) =  plot_Residual(timeloop,i,j,1);
                    r2(cellit) = plot_Residual(timeloop,i,j,2);
                    r3(cellit) = plot_Residual(timeloop,i,j,3);
                    r4(cellit) = plot_Residual(timeloop,i,j,4);
                    d1(cellit) = plot_cells_dissipation(timeloop,i,j,1);
                    d2(cellit) = plot_cells_dissipation(timeloop,i,j,2);
                    d3(cellit) = plot_cells_dissipation(timeloop,i,j,3);
                    d4(cellit) = plot_cells_dissipation(timeloop,i,j,4);
                    pstatic(cellit) =  plot_cells_pressure(timeloop,i,j); 
                    Mach(cellit) = sqrt((q2(cellit)^2+q3(cellit)^2)/(q1(cellit)^2));
                    c(cellit) = plot_cells_c(timeloop,i,j);
                    cellit = cellit+1;
                end
            end

            %Set time variable
            IterationsCount = (timeloop-1);
            IterationsName = num2str(IterationsCount);
                
            %write to the data file
            fprintf(fileID,' VARIABLES = "X", "Y", "q1", "q2", "q3", "q4","f1", "f2", "f3", "f4","g1","g2","g3","g4","r1","r2","r3","r4","d1","d2","d3","d4","Pressure","Mach","c" \n');
            fprintf(fileID,'ZONE T="%2s Iterations", I=%2g, J=%2g, DATAPACKING=BLOCK VARLOCATION=([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]=CELLCENTERED)\n',IterationsName,(nodes_Imax-4),(nodes_Jmax-4)); %[3,4,5,6] correspond to the number of cell-centered variables, in this case qvec
          
            for i = 1:(nodes_Imax-4)*(nodes_Jmax-4)
                fprintf(fileID,'%.4g ', x(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(nodes_Imax-4)*(nodes_Jmax-4)
                fprintf(fileID,'%.4g ', y(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', q1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', q2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', q3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', q4(i));
            end
                    fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', f1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', f2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', f3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', f4(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', g1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', g2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', g3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', g4(i));
            end
                    fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', r1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', r2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', r3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', r4(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', d1(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', d2(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', d3(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', d4(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', pstatic(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', Mach(i));
            end
            fprintf(fileID,'\n');
            for i = 1:(cells_Imax-4)*(cells_Jmax-4)
                fprintf(fileID,'%.4g ', c(i));
            end
            fprintf(fileID,'\n');

    end
end

    fclose(fileID);
end

