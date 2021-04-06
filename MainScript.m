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
user_itmax = 10;           %maximum number of iterations made when solving
user_tol = 0.00005;         %acceptable nondimensional error/tolerance of the residual when solving
v2 = 0.25;                  %[0,0.5] dissipation switch second order
v4 = 0.005;                %[0.0001,0.01] dissipation switch fourth order
CFL = 2*sqrt(2);
c = 1;                      %speed of sound is reference??

%% Input and modify the grid
%read node locations in from the specified grid, put into matraxies
[nodes_x_input,nodes_y_input] = LoadGrid(user_MeshQual); %no ghost nodes/cells

%Modify Grid by adding Ghost cells on the boundaries
[nodes_x,nodes_y] = addGhostCells(nodes_x_input,nodes_y_input);

%% Set up Matraxies
%Identify the Imax and Jmax, use to setup node and cell matraxies
nodes_Imax = length(nodes_x(:,1));     
nodes_Jmax = length(nodes_y(1,:));     
cells_Imax = nodes_Imax-1;
cells_Jmax = nodes_Jmax-1;
cells_q = zeros(cells_Imax,cells_Jmax,4);
cells_f = cells_q;
cells_g = cells_q;

%Setup A,R,D matraxies for cells
A = zeros(cells_Imax,cells_Jmax);
A = findAreas(nodes_x,nodes_y,cells_Imax,cells_Jmax,A); %fill out the unchanging area matrix
R = zeros(cells_Imax,cells_Jmax,4);
D = zeros(cells_Imax,cells_Jmax,4);

%% Grid Initialization

%find resulting static pressure from the desired resevoir pressure and mach
%number
P_resevoir = 1/user_Gamma;
[P_static] = pressurefinder(P_resevoir,user_Mach,user_Gamma);

%set initial conditions by filling out the q vector for every cell
[cells_q,cells_f,cells_g] = setInitialConditions(user_Mach,user_Gamma,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
%set boundary conditions after interior initial conditions
[cells_q,cells_f,cells_g] = applyBC(nodes_x,nodes_y,user_alpha,user_Gamma,user_Mach,P_resevoir,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
[cells_q,cells_f,cells_g] = cornersetter(cells_q,cells_f,cells_g,cells_Imax,cells_Jmax); %visual change only, used once
    
%% Iteration Loop for solving

%While loop that limits runtime based on tolerance and max iterations
iterations = 0;
residual_it = 1;   
residual_vec = [];
residualsum = 0;
residualavg_vec = [];
a_rk = [1/4,1/3,1/2,1]; %for rk loop
while (iterations<user_itmax) && (residual_it>user_tol) 
    residualmax = 0;
    for j = 3:(cells_Jmax-2) %loop through the interior cells in the grid
        for i = 3:(cells_Imax-2)
            %freeze a q
                q_freeze = squeeze(cells_q(i,j,:));
            %calculate and freeze a D
                [x_abcd,y_abcd] = nodes_touch_cell(i,j,nodes_x,nodes_y);
                q5by5 = cells_q((i-2):(i+2),(j-2):(j+2),:);
                p5by5 = cells_f((i-2):(i+2),(j-2):(j+2),2)-((cells_q((i-2):(i+2),(j-2):(j+2),2)).^2)./cells_q((i-2):(i+2),(j-2):(j+2),1);
                D_freeze = Dis(v2,v4,c,x_abcd,y_abcd,q5by5,p5by5);
                D(i,j,:) = D_freeze(:); %put the dissipation in a matrix for visualization
            %grab the cell area
                A_cell = A(i,j);
            %get the delta_t
                [dt] = time_step(c,CFL,A_cell,x_abcd,y_abcd,q5by5(2:4,2:4,:));
            %calculate f and g into the rk timestep thing
                 f_cNESW = squeeze([cells_f(i,j,:);cells_f(i,j+1,:);cells_f(i+1,j,:);cells_f(i,j-1,:);cells_f(i-1,j,:)]);
                 g_cNESW = squeeze([cells_g(i,j,:);cells_g(i,j+1,:);cells_g(i+1,j,:);cells_g(i,j-1,:);cells_g(i-1,j,:)]);                     
            %Runge-Kutta Temporal Incrementing
                [cells_q(i,j,:),cells_f(i,j,:),cells_g(i,j,:),cell_Res] = RK_time_step(user_Gamma,a_rk,q_freeze,D_freeze,A_cell,dt,x_abcd,y_abcd,f_cNESW,g_cNESW);
            %check if residual is larger than maximum grid residual
                if cell_Res>residualmax
                residualmax = cell_Res;
                end
            %add the cell residual to the total residual, used for mean
            %calculations
            residualsum = cell_Res+residualsum;
        end
    end
    
    %store residual for that iteration
    residual_vec = [residual_vec,residualmax];
    %average the summed residual
    residualavg = residualsum/(cells_Imax*cells_Jmax);
    residualavg_vec = [residualavg_vec,residualavg];
    %reapply BC
    [cells_q,cells_f,cells_g] = applyBC(nodes_x,nodes_y,user_alpha,user_Gamma,user_Mach,P_resevoir,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    %increase iteration count
    iterations = iterations+1;
end


%% Report Data

%plot residuals vs iteration number for the grid

%Format and export data to visualize in TecPlot
exportDataTecplot(user_Mach,iterations,nodes_x,nodes_y,cells_q,cells_f,cells_g,nodes_Imax,nodes_Jmax,cells_Imax,cells_Jmax);


