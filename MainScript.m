clear all;close all;clc;

%% Aero 489 Project 2: Euler Solver
%John Clark Pehrson
%Aero 489 Introduction to CFD with Prof. Cizmas
%March 21, 2021

%This program solves the euler equations for a geometry specified using an
%input laplace grid. The program allows the user to specify an input mach
%number and mesh quality, and reports the residuals, iteration number, bump
%forces, and mach contours for the final solution. Data visualization is
%handed off to TecPlot via a .txt file. 

%% User-Defined Variables
user_Mach = 0.3;            %choose either 0.3, 0.6, or 0.9
user_alpha = 0;             %direction of incoming flow into the inlet. Recommended to keep at 0. [deg]
user_Gamma = 1.4;           %the ratio of specific heats of the gas
user_MeshQual = 'coarse';   %choose either coarse, medium, or fine (or test for the algebraic test grid)
user_itmax = 3;            %maximum number of iterations made when solving
user_tol = 0.0000005;        %acceptable nondimensional error/tolerance of the residual when solving
v2 = 0.25;                  %[0,0.5] dissipation switch second order
v4 = 0.004;                 %[0.0001,0.01] dissipation switch fourth order
CFL = 1.25;                  %0.5 recommended from Cizmas
report_freq = 1;            %the frequency with which data is exported and plotted. Only for output
plot_full = 0;              %1 if visualize the ghost cells, 0 to not visualize ghost cells


%% Input and modify the grid
fwait = waitbar(0,'Loading and configuring data');
%read node locations in from the specified grid, put into matraxies
[nodes_x_input,nodes_y_input] = LoadGrid(user_MeshQual); %no ghost nodes/cells

%Modify Grid by adding Ghost cells on the boundaries
[nodes_x,nodes_y] = addGhostCells(nodes_x_input,nodes_y_input);

%% Initialize Matraxies of various sizes for data tracking and plotting
%Identify the Imax and Jmax for cells and nodes
nodes_Imax = length(nodes_x(:,1));     
nodes_Jmax = length(nodes_y(1,:));     
cells_Imax = nodes_Imax-1;
cells_Jmax = nodes_Jmax-1;

%Setup Residual plotting
meanresidual = NaN(4,user_itmax); %4 quantities to match residual (q vec change), one value for each iteration
maxresidual = NaN(4,user_itmax);

%Set up 2d matraxies for cells
cells_eig = NaN(cells_Imax,cells_Jmax,4);           %eigenvalues for individual cells. 1:N,2:E,3:S,4:W. Face directions for last index
cells_pressure = zeros(cells_Imax,cells_Jmax);      %static pressures for individual cells
cells_c = cells_pressure;                                %speed of sound for individual cells
A = zeros(cells_Imax,cells_Jmax);                   %initialize empty Area matrix
A = findAreas(nodes_x,nodes_y,cells_Imax,cells_Jmax,A); %fill out the unchanging area matrix

%Set up 3d matraxies for cells
cells_q = zeros(cells_Imax,cells_Jmax,4);
cells_f = cells_q;
cells_g = cells_q;
Residual = zeros(cells_Imax,cells_Jmax,4); %used to track residuals in a single iteration

%Set up 4d matraxies for cells over time (3d matraxies tracked for
%iterations)
plot_cells_q = NaN(1+user_itmax,cells_Imax,cells_Jmax,4); %The 1+plot_it_reduced is for bc, then the iterations
plot_cells_f = plot_cells_q;
plot_cells_g = plot_cells_q;
plot_Residual = plot_cells_q;
plot_cells_dissipation = zeros(user_itmax+1,cells_Imax,cells_Jmax,4);


%Set up 3d matraxies for cells over time (2d matraxies tracked for
%iterations)
plot_cells_pressure = NaN(1+user_itmax,cells_Imax,cells_Jmax);
plot_cells_c = plot_cells_pressure;

%Set up 2d matrix for bump force over time
plot_bump_force = NaN(1+user_itmax,2); % (:,1) is x, (:,2) is y

% Set up frozen variables for the loop iterations prior to RK
qfreeze_Mat = cells_q;
Dfreeze_Mat = cells_q;
dt_Mat = zeros(cells_Imax,cells_Jmax);
R_rk = zeros(cells_Imax,cells_Jmax,4);

%% Grid Initialization

%identify which bottom wall cells are on the bump
[cells_range_i] = wallbumpcellfinder(nodes_y,nodes_Imax);

%Set static pressure to be 1/gamma at inlet/outlet far field
P_static = 1/user_Gamma;
%set initial conditions by filling out the q vector for every cell
[cells_q,cells_f,cells_g] = setInitialConditions(user_Mach,user_Gamma,user_alpha,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax,nodes_x,nodes_y);

%set boundary conditions after interior initial conditions
    [cells_q,cells_f,cells_g] = applyBC(nodes_x,nodes_y,user_alpha,user_Gamma,user_Mach,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    [cells_q,cells_f,cells_g] = cornersetter(cells_q,cells_f,cells_g,cells_Imax,cells_Jmax); %visual change only, used once
    [cells_pressure,cells_c] = update_p_c(user_Gamma,cells_pressure,cells_c,cells_q,cells_f,cells_Imax,cells_Jmax);
    
        %update eigenvalues for interior cells
        for i = 3:(cells_Imax-2)
            for j = 3:(cells_Jmax-2)
                [x_abcd,y_abcd] = nodes_touch_cell(i,j,nodes_x,nodes_y);
                cells_eig(i,j,:) = eigenvaluefinder(x_abcd,y_abcd,cells_q(i,j,:),cells_c(i,j));
            end
        end
        %update boundary eigenvalues for ghost cells
        cells_eig = eig_BC(cells_Imax,cells_Jmax,cells_eig);
        
%calculate the force on the bump
[plot_bump_force(1,:)] = bump_force_calculator(cells_pressure,nodes_x,nodes_y,cells_range_i);
    
%save the initialized grid for visualization    
plot_cells_q(1,:,:,:) = cells_q; %save data for visualization
plot_cells_f(1,:,:,:) = cells_f;
plot_cells_g(1,:,:,:) = cells_g;
plot_Residual(1,:,:,:) = zeros(cells_Imax,cells_Jmax,4);
plot_cells_pressure(1,:,:) = cells_pressure;
plot_cells_c(1,:,:) = cells_c;
plot_i = 2;

%% Iteration Loop for solving

 %Define some preliminary variables/constants for the while loop
 iterations = 1;
 a_rk = [1/4,1/3,1/2,1]; %for rk loop
 
 while (iterations<(user_itmax+1)) %while loop that iterates the solution in time
    waitbar(iterations/user_itmax,fwait,'Iterating');
     
    %loops through cells for the frozen data
    for j = 3:(cells_Jmax-2)
        for i = 3:(cells_Imax-2)
            %freeze a q
                qfreeze_Mat(i,j,:) = squeeze(cells_q(i,j,:));
            %calculate and freeze a D
                [x_abcd,y_abcd] = nodes_touch_cell(i,j,nodes_x,nodes_y);
                q5by5 = cells_q((i-2):(i+2),(j-2):(j+2),:);
                p5by5 = cells_pressure((i-2):(i+2),(j-2):(j+2));
                eig3by3 = cells_eig((i-1):(i+1),(j-1):(j+1),:);
                Dfreeze_Mat(i,j,:) = Dis(v2,v4,user_Gamma,x_abcd,y_abcd,q5by5,p5by5,eig3by3);           
                    %save dissipation for visualization
                    plot_cells_dissipation(plot_i,i,j,:) = Dfreeze_Mat(i,j,:);
            %grab the cell area
                A_cell = A(i,j);
            %get the delta_t
                dt_Mat(i,j) = time_step(CFL,A_cell,x_abcd,y_abcd,eig3by3);
        end
    end

    %replacement RK loop
    for k = 1:4 %loops through 4 rk steps
        
        for j = 3:(cells_Jmax-2) %loop through cells
            for i = 3:(cells_Imax-2)
                %pseudo step a single cell (find residual then find the
                %pseudo step q
                f_cNESW = squeeze([cells_f(i,j,:);cells_f(i,j+1,:);cells_f(i+1,j,:);cells_f(i,j-1,:);cells_f(i-1,j,:)]);
                g_cNESW = squeeze([cells_g(i,j,:);cells_g(i,j+1,:);cells_g(i+1,j,:);cells_g(i,j-1,:);cells_g(i-1,j,:)]);
                R_rk(i,j,:) = findResidual(x_abcd,y_abcd,f_cNESW,g_cNESW);
                cells_q(i,j,:) = qfreeze_Mat(i,j,:)-a_rk(k)*(dt_Mat(i,j)/A(i,j))*(R_rk(i,j,:)-Dfreeze_Mat(i,j,:));
                
                %update the f and g matrix?
                [cells_f(i,j,:),cells_g(i,j,:)] = refresh_f_g(cells_q(i,j,:),user_Gamma);
            end
        end
        
        %apply BC after q is updated
        [cells_q,cells_f,cells_g] = applyBC(nodes_x,nodes_y,user_alpha,user_Gamma,user_Mach,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    end
    Residual = R_rk;
    
    
    %loops through cells to update p, c, eigenvalues (which aren't used
    %again until the next iteration
    for j = 3:(cells_Jmax-2)
        for i = 3:(cells_Imax-2)
            %update the eigenvalues, pressure, and c for the cell
            cells_pressure(i,j) = cells_f(i,j,2)-(cells_q(i,j,2)^2)/cells_q(i,j,1);
            cells_c(i,j) = sqrt(user_Gamma*cells_pressure(i,j)/cells_q(i,j,1));
            [x_abcd,y_abcd] = nodes_touch_cell(i,j,nodes_x,nodes_y);
            cells_eig(i,j,:) = eigenvaluefinder(x_abcd,y_abcd,cells_q(i,j,:),cells_c(i,j));
        end
    end
    
    %calculate the force on the bump
    [plot_bump_force(plot_i,:)] = bump_force_calculator(cells_pressure,nodes_x,nodes_y,cells_range_i);
    
    %residual tracking for a single iteration
    spani = 3:(cells_Imax-2);
    spanj = 3:(cells_Jmax-2);
    meanresidual(:,iterations) = [abs(mean(mean(Residual(spani,spanj,1))));abs(mean(mean(Residual(spani,spanj,2)))); abs(mean(mean(Residual(spani,spanj,3))));abs(mean(mean(Residual(spani,spanj,4))))];
    maxresidual(:,iterations) = [max(max(abs(Residual(spani,spanj,1))));max(max(abs(Residual(spani,spanj,2))));max(max(abs(Residual(spani,spanj,3))));max(max(abs(Residual(spani,spanj,4))));];

    %reapply BC
    [cells_q,cells_f,cells_g] = applyBC(nodes_x,nodes_y,user_alpha,user_Gamma,user_Mach,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    
    %update p,c,eigenvalues
    [cells_pressure,cells_c] = update_p_c(user_Gamma,cells_pressure,cells_c,cells_q,cells_f,cells_Imax,cells_Jmax);
    cells_eig = eig_BC(cells_Imax,cells_Jmax,cells_eig);
    
    %Save data for plotting
        plot_cells_q(plot_i,:,:,:) = cells_q; %save data for visualization
        plot_cells_f(plot_i,:,:,:) = cells_f;
        plot_cells_g(plot_i,:,:,:) = cells_g;
        plot_Residual(plot_i,:,:,:) = Residual;
        plot_cells_pressure(plot_i,:,:) = cells_pressure;
        plot_cells_c(plot_i,:,:) = cells_c;
        plot_i = plot_i+1;

    
    %increase iteration count
    iterations = iterations+1;
   
end


%% Report Data

waitbar(1,fwait,'exporting data');
%Format and export data to visualize in TecPlot
exportDataTecplot(user_Mach,user_MeshQual,iterations,nodes_x,nodes_y,plot_cells_q,plot_cells_f,plot_cells_g,plot_Residual,plot_cells_pressure,plot_cells_c,plot_cells_dissipation,nodes_Imax,nodes_Jmax,cells_Imax,cells_Jmax,user_itmax,report_freq,plot_full);

%% Plot Residuals and bump force
figure;
for i = 1:4
plot(1:report_freq:user_itmax,log10(meanresidual(i,1:report_freq:user_itmax)),'Linewidth',2);
hold on;
end
title('Mean Residuals');
legend('R rho','R rho*u','R rho*v','R rho*E');
xlabel('Iteration #');
ylabel('Log Base 10 of Mean Residual');

figure;
for i = 1:4
plot(1:report_freq:user_itmax,log10(maxresidual(i,1:report_freq:user_itmax)),'Linewidth',2);
hold on;
end
title('Max Residuals');
legend('R rho','R rho*u','R rho*v','R rho*E');
xlabel('Iteration #');
ylabel('Log Base 10 of Max Residual');

figure;
yyaxis left;
plot(1:report_freq:user_itmax+1,plot_bump_force(1:report_freq:user_itmax+1,1),'Linewidth',2);
hold on;
yyaxis right;
plot(1:report_freq:user_itmax+1,plot_bump_force(1:report_freq:user_itmax+1,2),'Linewidth',2);
title('Forces at the bump');
xlabel('Iteration number');

yyaxis left;
ylabel('Nondimensional Force in the X direction');
yyaxis right;
ylabel('Nondimensional Force in the Y direction');
legend('Force in X','Force in Y');

close(fwait);