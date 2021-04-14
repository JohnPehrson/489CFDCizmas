function [cells_q,cells_f,cells_g,cells_pressure,cells_eig,cells_c] = setInitialConditions(user_Mach,user_Gamma,user_alpha,P_static,cells_q,cells_f,cells_g,cells_Imax,cells_Jmax,nodes_x,nodes_y)
%This function defines the initial conditions of the flow
%Assume 2 layers of ghost nodes
%Using info from page 29 of my handwritten notes

%constant parameters for the flow (assume Mach is constant, but that the
%alpha is variable depending on the location in the geometry)
rho = 1;
p = P_static; 
E = p/((user_Gamma-1)*rho)+0.5*user_Mach^2;
H = E+p/rho;
c = sqrt(user_Gamma*p/rho);

%loop for individual cells (loop j outside, for i in that loop)
for j = 3:(cells_Jmax-2)
    for i = 3:(cells_Imax-2)
        %assign alpha,u,v depending on the cell
        [AoAtangent] = bottomwalltangent(nodes_x(i,j)); %find AoA in rad of the flow at that point
        alpha = deg2rad(user_alpha)+AoAtangent*(1-nodes_y(i,j)).^2; %AoA component to account for wall, nodes_y to damp for the closeness to the top wall
        u = user_Mach*cos(alpha);
        v = user_Mach*sin(alpha);
        
            q = [rho;rho*u;rho*v;rho*E];
            f = [rho*u;rho*u^2+p;rho*u*v;rho*H*u];
            g = [rho*v;rho*u*v;rho*v^2+p;rho*H*v];
        
            cells_q(i,j,:) = q;
            cells_f(i,j,:) = f;
            cells_g(i,j,:) = g;
            
    end
end


end

