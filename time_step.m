function [delta_t] = time_step(CFL,cell_area,x_abcd,y_abcd,ei3by3)
%This function calculates the time step delta_t at a cell for a step in
%runge kutta. Because the simulation is steady, time steps for individual
%cells doesn't need to be the same

%Inputs: Area of the cell
    %CFL number: CFL, should be less than 2sqrt(2) for stability
    %nearby node locations: x_abcd and y_abcd
    %q3by3 for the face velocies in the calculation of the eigenvalues
    %ei3by3 is the matrix of local eigenvalues 
%Outputs: The time step delta_t
    
i = 2;
j = 2;
lambda_xi_E = 0.5*(ei3by3(i,j,1)+ei3by3(i+1,j,1));
lambda_xi_W = 0.5*(ei3by3(i,j,1)+ei3by3(i-1,j,1));
lambda_n_N = 0.5*(ei3by3(i,j,2)+ei3by3(i,j+1,2));
lambda_n_S = 0.5*(ei3by3(i,j,2)+ei3by3(i,j-1,2));

sumeigenvaluefaces = lambda_n_N*edgelength(x_abcd,y_abcd,'N')+...
                     lambda_xi_E*edgelength(x_abcd,y_abcd,'E')+...
                     lambda_n_S*edgelength(x_abcd,y_abcd,'S')+...
                     lambda_xi_W*edgelength(x_abcd,y_abcd,'W');

                 
delta_t_max = 2*cell_area/sumeigenvaluefaces;
delta_t = delta_t_max/CFL;
end

