function [q_b,f_b,g_b] = singleCellInletBC(user_Gamma,user_Mach,P_static,alpha,q_in,f_in,g_in)
%This sub-function sets the ghost cells outside the inlet of the grid.
%The program uses Riemann invariants, pressure, and flow AoA to find the
%ghost cells.
%alpha is in degrees
%P_static is freestream static pressure

%q_in, f_in, and g_in all reference the first cell inside the wall
%q_b, f_b, and g_b reference the TWO ghost cells that act as the boundary

%Implicitly assume that the boundary normal is perfectly aligned with 'x'
%direction

%find R-, which is free stream riemann invariant
uinf = -user_Mach*cosd(alpha);
rhoinf = 1;
cinf = sqrt(user_Gamma*P_static/rhoinf);
R_neg = uinf-2*cinf/(user_Gamma-1);

%find R+, which is the interior-flow riemann invariant
rhoc = q_in(1);
uc = -q_in(2)/rhoc;
pc = f_in(2)-rhoc*uc^2;
cc = sqrt(user_Gamma*pc/rhoc);
R_pos = uc+2*cc/(user_Gamma-1);

%find velocity at inlet using riemans
Vbmag = 0.5*(R_neg+R_pos);
ub = -Vbmag*cos(alpha);
vb = -Vbmag*sin(alpha);
cb = 0.25*(user_Gamma-1)*(R_pos-R_neg);
sb = cinf^2/(user_Gamma*(rhoinf^(user_Gamma-1)));
rhob = ((cb^2)/(user_Gamma*sb))^(1/(user_Gamma-1));
pb = (rhob*cb^2)/user_Gamma;
Eb = pb/((user_Gamma-1)*rhob)+0.5*(ub^2+vb^2);


%Put variables into the q,f,g form
q_b_cell = [rhob;rhob*ub;rhob*vb;rhob*Eb];
f_b_cell = [rhob*ub;rhob*ub^2+pb;rhob*ub*vb;rhob*(Eb+pb)*ub];
g_b_cell = [rhob*vb;rhob*ub*vb;rhob*vb^2+pb;rhob*(Eb+pb)*vb];

%create empty output matrixes in the right format
q_b = NaN(2,1,4);
f_b = q_b;
g_b = q_b;

%convert the gc into the right format for output
q_b(1,:,:) = q_b_cell;
q_b(2,:,:) = q_b_cell;
f_b(1,:,:) = f_b_cell;
f_b(2,:,:) = f_b_cell;
g_b(1,:,:) = g_b_cell;
g_b(2,:,:) = g_b_cell;


if pb<0.5
    fprintf('whats up');
end

end

