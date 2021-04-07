function [delta_t] = time_step(user_Gamma,CFL,cell_area,x_abcd,y_abcd,q3by3,p3by3)
%This function calculates the time step delta_t at a cell for a step in
%runge kutta. Because the simulation is steady, time steps for individual
%cells doesn't need to be the same

%Inputs: Area of the cell
    %CFL number: CFL, should be less than 2sqrt(2) for stability
    %nearby node locations: x_abcd and y_abcd
    %speed of sound: c
    %q3by3 for the face velocies in the calculation of the eigenvalues
%Outputs: The time step delta_t


sumeigenvaluefaces = eigenvalueDis(user_Gamma,x_abcd,y_abcd,q3by3,p3by3,'N')*edgelength(x_abcd,y_abcd,'N')+...
                     eigenvalueDis(user_Gamma,x_abcd,y_abcd,q3by3,p3by3,'E')*edgelength(x_abcd,y_abcd,'E')+...
                     eigenvalueDis(user_Gamma,x_abcd,y_abcd,q3by3,p3by3,'S')*edgelength(x_abcd,y_abcd,'S')+...
                     eigenvalueDis(user_Gamma,x_abcd,y_abcd,q3by3,p3by3,'W')*edgelength(x_abcd,y_abcd,'W');

                 
delta_t_max = 2*cell_area/sumeigenvaluefaces;
delta_t = delta_t_max/CFL;
end

