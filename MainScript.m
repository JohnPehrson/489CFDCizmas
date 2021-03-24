clear all;close all;clc;

%% Aero 489 Project 2: Euler Solver
%John Clark Pehrson
%Aero 489 Introduction to CFD with Prof. Cizmas
%March 21, 2021

%This program solves the euler equations for a geometry specified using an
%input laplace grid. The program allows the user to specify an input mach
%number and mesh quality, and reports the residuals, iteration number, bump
%forces, and mach contours for the final solution.

%% To Do:
%make ghost nodes near the input and output
%look into LaTex
%Write subfunctions for each value in the Dij calculation (switches,
%lengths, eivenvalues, f(p), etc.)



%% User-Defined Variables
user_Mach = 0.3;            %choose either 0.3, 0.6, or 0.9
user_MeshQual = 'coarse';   %choose either coarse, medium, or fine
user_itmax = 100;           %maximum number of iterations made when solving
user_tol = 0.00005;         %acceptable nondimensional error/tolerance of the residual when solving
v2 = 0.25;                  %[0,0.5] dissipation switch second order
v4 = 0.005;                %[0.0001,0.01] dissipation switch fourth order


%% Set up Matraxies

%read node locations in from the specified grid, put into matraxies
[nodes_x,nodes_y] = LoadGrid(user_MeshQual);

%Identify the Imax and Jmax, use to setup node and cell matraxies
nodes_Imax = length(nodes_x(:,1));
nodes_Jmax = length(nodes_y(1,:));
cells_Imax = nodes_Imax-1;
cells_Jmax = nodes_Jmax-1;
cells_q = zeros(cells_Imax,cells_Jmax,5);
cells_f = cells_q;
cells_g = cells_q;

%Setup A,R,D matraxies for cells
A = zeros(cells_Imax,cells_Jmax);
A = findAreas(nodes_x,nodes_y,cells_Imax,cells_Jmax,A); %fill out the unchanging area matrix
R = zeros(cells_Imax,cells_Jmax);
%R = findResiduals(nodes_x,nodes_y,cells_f,cells_g,cells_Imax,cells_Jmax);
D = zeros(cells_Imax,cells_Jmax);

%identify the type of cells for eventual flux stuff
[celltype] = cellidentifier(cells_Imax,cells_Jmax);

%% Grid Initialization

%set initial conditions by filling out the q vector for every cell

%% Iteration Loop for solving

%While loop that limits runtime based on tolerance and max iterations

%While loop iterates through time, if statements to loop through individual
%cells...?

%call function to do node logic to identify the location of the node in the
%grid

%timestep limitations?


%% Report Data

%plot residuals vs iteration number for the grid

%plot the forces at the bump in the x and y direction as a function of the
%number of points in the grid

%plot mach number contours


