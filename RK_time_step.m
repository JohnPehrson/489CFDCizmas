function [cell_q_out,cell_f_out,cell_g_out,Resmax] = RK_time_step(user_Gamma,a_rk,q_freeze,D_freeze,A_cell,dt,x_abcd,y_abcd,f_cNESW,g_cNESW)
%This is a modified Runge-Kutta algorithm that steps from n to n+1 (iterations)
%This algorithm is called at a single cell at a time n, meaning that all
%solutions should be reported only for that cell

%D_freeze is a 1x1x4
%q_freeze is a 1x1x4
%f_cNESW and g_cNESW are both 5x4

q_rk = zeros(4,4);
R_rk = zeros(4,4);
for k = 1:4 %loops through the 4 steps in the RK algorithm
    R_rk(:,k) = findResidual(x_abcd,y_abcd,f_cNESW,g_cNESW);
    q_rk(:,k) = q_freeze(:)-a_rk(k)*(dt/A_cell)*(R_rk(:,k)-D_freeze(:));
    
        %update f,g using q and pressure
        [f,g] = refresh_f_g(q_rk(:,k),user_Gamma);
        f_cNESW(1,:) = f';
        g_cNESW(1,:) = g';
end
cell_q_out = q_rk(:,4);
cell_f_out = f_cNESW(1,:)';
cell_g_out = g_cNESW(1,:)';
Resmax = max(max(R_rk));
end

