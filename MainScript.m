clear all;close all;clc;

%% Aero 489 Project 2: Euler Solver
%John Clark Pehrson
%Aero 489 Introduction to CFD with Prof. Cizmas
%March 21, 2021

%This program solves the euler equations for a geometry specified using an
%input laplace grid. The program allows the user to specify an input mach
%number and mesh quality, and reports the residuals, iteration number, bump
%forces, and mach contours for the final solution.

%% User-Defined Variables
user_Mach = 0.3;            %choose either 0.3, 0.6, or 0.9
user_alpha = 0;             %direction of incoming flow into the inlet. Recommended to keep at 0. [deg]
user_Gamma = 1.4;           %the ratio of specific heats of the gas
user_MeshQual = 'coarse';   %choose either coarse, medium, or fine
user_itmax = 100;           %maximum number of iterations made when solving
user_tol = 0.00005;         %acceptable nondimensional error/tolerance of the residual when solving
v2 = 0.25;                  %[0,0.5] dissipation switch second order
v4 = 0.005;                %[0.0001,0.01] dissipation switch fourth order

%% Input and modify the grid
%read node locations in from the specified grid, put into matraxies
[nodes_x_input,nodes_y_input] = LoadGrid(user_MeshQual); %no ghost nodes/cells

%Modify Grid by adding Ghost cells on the boundaries
[nodes_x,nodes_y] = addGhostCells(nodes_x_input,nodes_y_input);

% %% Set up Matraxies
%Identify the Imax and Jmax, use to setup node and cell matraxies
nodes_Imax = length(nodes_x(:,1));     
nodes_Jmax = length(nodes_y(1,:));     
cells_Imax = nodes_Imax-1;
cells_Jmax = nodes_Jmax-1;
cells_q = NaN(cells_Imax,cells_Jmax,4);
cells_f = cells_q;
cells_g = cells_q;

% %Setup A,R,D matraxies for cells
% A = zeros(cells_Imax,cells_Jmax);
% A = findAreas(nodes_x,nodes_y,cells_Imax,cells_Jmax,A); %fill out the unchanging area matrix
% R = zeros(cells_Imax,cells_Jmax);
% %R = findResiduals(nodes_x,nodes_y,cells_f,cells_g,cells_Imax,cells_Jmax);
% D = zeros(cells_Imax,cells_Jmax);
% 

%% Grid Initialization

%find resulting static pressure from the desired resevoir pressure and mach
%number
P_resevoir = 1/user_Gamma;
[P_static] = pressurefinder(P_resevoir,user_Mach,user_Gamma);

%set initial conditions by filling out the q vector for every cell
[cells_q,cells_f,cells_g] = setInitialConditions(user_Mach,user_Gamma,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    %set boundary conditions after interior initial conditions
    [cells_q,cells_f,cells_g] = applyBottomWallBC(nodes_x,nodes_y,cells_q,cells_f,cells_g,cells_Imax);
    [cells_q,cells_f,cells_g] = applyUpperWallBC(nodes_x,nodes_y,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    [cells_q,cells_f,cells_g] = applyOutletBC(user_Gamma,P_resevoir,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    [cells_q,cells_f,cells_g] = applyInletBC(user_alpha,user_Gamma,user_Mach,P_resevoir,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax);
    
% %% Iteration Loop for solving
% 
% %While loop that limits runtime based on tolerance and max iterations
% 
% %While loop iterates through time, if statements to loop through individual
% %cells...?
% 
% %call function to do node logic to identify the location of the node in the
% %grid
% 
% %timestep limitations?
% 
% 
% %% Report Data
% 
% %plot residuals vs iteration number for the grid
% 
% %plot the forces at the bump in the x and y direction as a function of the
% %number of points in the grid
% 
% %plot mach number contours
% 
% 
