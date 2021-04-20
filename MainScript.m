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
user_MeshQual = 'coarse';   %choose either coarse, medium, or fine
user_itmax = 300;            %maximum number of iterations made when solving
user_tol = 0.00005;         %acceptable nondimensional error/tolerance of the residual when solving
v2 = 0.25;                  %[0,0.5] dissipation switch second order
v4 = 0.004;                 %[0.0001,0.01] dissipation switch fourth order
CFL = 2;                   %0.5 recommended from Cizmas


%% Input and modify the grid
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
cells_eig = zeros(cells_Imax,cells_Jmax,2);           %eigenvalues for individual cells
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

%Set up 3d matraxies for cells over time (2d matraxies tracked for
%iterations)
plot_cells_pressure = NaN(1+user_itmax,cells_Imax,cells_Jmax);
plot_cells_c = plot_cells_pressure;


%% Grid Initialization

%Set static pressure to be 1/gamma at inlet/outlet far field
P_static = 1/user_Gamma;
%set initial conditions by filling out the q vector for every cell
[cells_q,cells_f,cells_g] = setInitialConditions(user_Mach,user_Gamma,user_alpha,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax,nodes_x,nodes_y);

%set boundary conditions after interior initial conditions
    [cells_q,cells_f,cells_g] = applyBC(nodes_x,nodes_y,user_alpha,user_Gamma,user_Mach,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    [cells_q,cells_f,cells_g] = cornersetter(cells_q,cells_f,cells_g,cells_Imax,cells_Jmax); %visual change only, used once
    [cells_pressure,cells_c,cells_eig] = update_p_c_eig(user_Gamma,cells_pressure,cells_c,cells_eig,cells_q,cells_f,nodes_x,nodes_y,cells_Imax,cells_Jmax);
    
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
 residual_it = 1; %just to get the loop going
 a_rk = [1/4,1/3,1/2,1]; %for rk loop
 
 while (iterations<(user_itmax+1)) && (residual_it>user_tol) %while loop that iterates the solution in time

    residualmax = 0;
    for j = 3:(cells_Jmax-2) %loop through the interior cells in the grid
        for i = 3:(cells_Imax-2)
            %freeze a q
                q_freeze = squeeze(cells_q(i,j,:));
            %calculate and freeze a D
                [x_abcd,y_abcd] = nodes_touch_cell(i,j,nodes_x,nodes_y);
                q5by5 = cells_q((i-2):(i+2),(j-2):(j+2),:);
                p5by5 = cells_f((i-2):(i+2),(j-2):(j+2),2)-((cells_q((i-2):(i+2),(j-2):(j+2),2)).^2)./cells_q((i-2):(i+2),(j-2):(j+2),1);
                eig3by3 = cells_eig((i-1):(i+1),(j-1):(j+1),:);
                D_freeze = Dis(v2,v4,user_Gamma,x_abcd,y_abcd,q5by5,p5by5,eig3by3);
            %grab the cell area
                A_cell = A(i,j);
            %get the delta_t
                dt = time_step(CFL,A_cell,x_abcd,y_abcd,eig3by3);
            %calculate f and g into the rk timestep thing
                 f_cNESW = squeeze([cells_f(i,j,:);cells_f(i,j+1,:);cells_f(i+1,j,:);cells_f(i,j-1,:);cells_f(i-1,j,:)]);
                 g_cNESW = squeeze([cells_g(i,j,:);cells_g(i,j+1,:);cells_g(i+1,j,:);cells_g(i,j-1,:);cells_g(i-1,j,:)]);                     
            %Runge-Kutta Temporal Incrementing
                [cells_q(i,j,:),cells_f(i,j,:),cells_g(i,j,:),Residual(i,j,:)] = RK_time_step(user_Gamma,a_rk,q_freeze,D_freeze,A_cell,dt,x_abcd,y_abcd,f_cNESW,g_cNESW);
        
            %update the eigenvalues, pressure, and c for the cell
            cells_pressure(i,j) = cells_f(i,j,2)-(cells_q(i,j,2)^2)/cells_q(i,j,1);
            cells_c(i,j) = sqrt(user_Gamma*cells_pressure(i,j)/cells_q(i,j,1));
            [x_abcd,y_abcd] = nodes_touch_cell(i,j,nodes_x,nodes_y);
            [cells_eig(i,j,1),cells_eig(i,j,2)] = eigenvaluefinder(x_abcd,y_abcd,cells_q(i,j,:),cells_c(i,j));
            
        end
    end
    
    %residual tracking for a single iteration
    spani = 3:(cells_Imax-2);
    spanj = 3:(cells_Jmax-2);
    meanresidual(:,iterations) = [abs(mean(mean(Residual(spani,spanj,1))));abs(mean(mean(Residual(spani,spanj,2)))); abs(mean(mean(Residual(spani,spanj,3))));abs(mean(mean(Residual(spani,spanj,4))))];
    maxresidual(:,iterations) = [max(max(Residual(spani,spanj,1)));max(max(Residual(spani,spanj,2)));max(max(Residual(spani,spanj,3)));max(max(Residual(spani,spanj,4)));];

    %reapply BC
    [cells_q,cells_f,cells_g] = applyBC(nodes_x,nodes_y,user_alpha,user_Gamma,user_Mach,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    
    %update p,c,eigenvalues
    [cells_pressure,cells_c,cells_eig] = update_p_c_eig(user_Gamma,cells_pressure,cells_c,cells_eig,cells_q,cells_f,nodes_x,nodes_y,cells_Imax,cells_Jmax);
    
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

%Format and export data to visualize in TecPlot
exportDataTecplot(user_Mach,iterations,nodes_x,nodes_y,plot_cells_q,plot_cells_f,plot_cells_g,plot_Residual,plot_cells_pressure,plot_cells_c,nodes_Imax,nodes_Jmax,cells_Imax,cells_Jmax,user_itmax);

%% Plot Residuals
figure;
for i = 1:4
plot(1:user_itmax,meanresidual(i,:),'Linewidth',2);
hold on;
end
title('Mean Residuals');
legend('Res_rho','Res_rho*u','Res_rho*v','Res_rho*E');
xlabel('Iteration #');
ylabel('Mean Residual');

figure;
for i = 1:4
plot(1:user_itmax,maxresidual(i,:),'Linewidth',2);
hold on;
end
title('Max Residuals');
legend('Res_rho','Res_rho*u','Res_rho*v','Res_rho*E');
xlabel('Iteration #');
ylabel('Max Residual');

