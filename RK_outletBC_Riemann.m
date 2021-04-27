function [f_gc,g_gc] = RK_outletBC_Riemann(user_Gamma,user_Mach,P_static,q_in,f_in,g_in)
%This is a function that modifies f and g in the RK loop, to be called when
%the RK cell is bordering the outlet. The f and g for the ghost cell are
%sent back to the RK loop.

%% R- is cells, R+ is infinity


%R-
uinf = user_Mach;
rhoinf = 1;
pinf = P_static;
cinf = sqrt(user_Gamma*pinf/rhoinf);
R_neg = uinf-2*cinf/(user_Gamma-1);

%R+
rhoc = q_in(1);
uc = q_in(2)/rhoc;
pc = f_in(2)-rhoc*uc^2;
cc = sqrt(user_Gamma*pc/rhoc);
R_pos = uc+2*cc/(user_Gamma-1);

%Calculate ub and cb
ub = 0.5*(R_pos+R_neg);
cb = 0.25*(user_Gamma-1)*(R_pos-R_neg);

%calculate other stuff
sb = (cc^2)/(user_Gamma*(rhoc)^(user_Gamma-1)); %for outflow
rhob = ((cb^2)/(user_Gamma*sb))^((1)/(user_Gamma-1));
pb = rhob*cb^2/user_Gamma;

%getting everything in the right format
%rhob
%ub
vb = 0;
Eb = pb/((user_Gamma-1)*rhob)+0.5*(ub^2+vb^2);
Hb = Eb+pb/rhob;

%q_out1cell = [rhob,rhob*ub,rhob*vb,rhob*Eb];
f_gc = [rhob*ub,rhob*ub^2+pb,rhob*ub*vb,rhob*Hb*ub];
g_gc = [rhob*vb,rhob*ub*vb,rhob*vb^2+pb,rhob*Hb*vb];



end

