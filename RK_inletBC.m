function [f_gc,g_gc] = RK_inletBC(user_Gamma,alpha,user_Mach,P_static,q_in,f_in,g_in)
%This function implements the inlet BC on the adjacent ghost cell when
%doing the RK time loop. Called from the Rk loop only.

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
%q_b_cell = [rhob;rhob*ub;rhob*vb;rhob*Eb];
f_gc = [rhob*ub;rhob*ub^2+pb;rhob*ub*vb;rhob*(Eb+pb/rhob)*ub];
g_gc = [rhob*vb;rhob*ub*vb;rhob*vb^2+pb;rhob*(Eb+pb/rhob)*vb];

end

