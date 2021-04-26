function [f_cNESW_out,g_cNESW_out] = RK_BC(user_Gamma,user_Mach,user_alpha,P_static,ind1,ind2,q_cell,f_cNESW,g_cNESW,x_abcd,y_abcd)
%This function will modify the fluxes near the cell in the RK loop (during
%the pseudo-time-stepping) to obey the boundary conditions it may have.
   
    %define the output matrix size
    f_cNESW_out = NaN(5,4);
    g_cNESW_out = NaN(5,4);
    
    %set the central cell, it never changes because of BC
    f_cNESW_out(1,:) = f_cNESW(1,:);
    g_cNESW_out(1,:) = g_cNESW(1,:);

    %% Logic for BC
    switch ind1
        case 'E' %outlet BC
            [f_cNESW_out(3,:) ,g_cNESW_out(3,:)] = RK_outletBC(user_Gamma,user_Mach,P_static,q_cell,f_cNESW(1,:),g_cNESW(1,:));
            f_cNESW_out(5,:) = f_cNESW(5,:);
            g_cNESW_out(5,:) = g_cNESW(5,:);
        case 'W' %inlet BC
            f_cNESW_out(3,:) = f_cNESW(3,:);
            g_cNESW_out(3,:) = g_cNESW(3,:);
            [f_cNESW_out(5,:),g_cNESW_out(5,:)] = RK_inletBC(user_Gamma,user_alpha,user_Mach,P_static,q_cell,f_cNESW(1,:),g_cNESW(1,:));
        case 'neither' %no BC in E or W
            f_cNESW_out(3,:) = f_cNESW(3,:);
            g_cNESW_out(3,:) = g_cNESW(3,:);
            f_cNESW_out(5,:) = f_cNESW(5,:);
            g_cNESW_out(5,:) = g_cNESW(5,:);
    end
    
    switch ind2
        case 'N' %top wall BC
            walltan = [x_abcd(3)-x_abcd(4);y_abcd(3)-y_abcd(4)];
            [f_cNESW_out(2,:),g_cNESW_out(2,:)] = RK_wallBC(user_Gamma,q_cell,f_cNESW(1,:),g_cNESW(1,:),walltan);
            f_cNESW_out(4,:) = f_cNESW(4,:);
            g_cNESW_out(4,:) = g_cNESW(4,:);
        case 'S' %bottom wall BC
            f_cNESW_out(2,:) = f_cNESW(2,:);
            g_cNESW_out(2,:) = g_cNESW(2,:);
            walltan = [x_abcd(2)-x_abcd(1);y_abcd(2)-y_abcd(1)];
            [f_cNESW_out(4,:),g_cNESW_out(4,:)] = RK_wallBC(user_Gamma,q_cell,f_cNESW(1,:),g_cNESW(1,:),walltan);
        case 'neither' %no BC in N or S
            f_cNESW_out(2,:) = f_cNESW(2,:);
            g_cNESW_out(2,:) = g_cNESW(2,:);
            f_cNESW_out(4,:) = f_cNESW(4,:);
            g_cNESW_out(4,:) = g_cNESW(4,:);
    end

end

