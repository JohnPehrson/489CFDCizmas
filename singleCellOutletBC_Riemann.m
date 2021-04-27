function [q_out,f_out,g_out] = singleCellOutletBC_Riemann(user_Gamma,user_Mach,q_in,f_in,g_in,P_static)
%This function sets the outlet boundary conditons in a single cell. Used in
%"applyOutletBC".

%q_in,f_in, and g_in are all 2x1x4 vectors, and contain the information for
%both the right-most and second-right-most cells in the flow

%q_out,f_out,g_out are also 2x1x4 and represent the ghost cells in the same
%row as the inputs

% P_static is freestream pressure of the flow at infinity

%% R- is cells, R+ is infinity

%old
%R-
uinf = user_Mach;
rhoinf = 1;
pinf = P_static;
cinf = sqrt(user_Gamma*pinf/rhoinf);
R_neg = uinf-2*cinf/(user_Gamma-1);

%R+
rhoc = q_in(1,1,1);
uc = q_in(1,1,2)/rhoc;
pc = f_in(1,1,2)-rhoc*uc^2;
cc = sqrt(user_Gamma*pc/rhoc);
R_pos = uc+2*cc/(user_Gamma-1);

%Calculate ub and cb
ub = 0.5*(R_pos+R_neg);
cb = 0.25*(user_Gamma-1)*(R_pos-R_neg);

%calculate other stuff
sb = (cc^2)/(user_Gamma*(rhoc)^(user_Gamma-1)); %for outflow
rhob = ((cb^2)/(user_Gamma*sb))^((1)/(user_Gamma-1));
pb = (rhob*cb^2)/user_Gamma;

%getting everything in the right format
%rhob
%ub
vb = 0;
Eb = pb/((user_Gamma-1)*rhob)+0.5*(ub^2+vb^2);
Hb = Eb+pb/rhob;

q_out1cell = [rhob,rhob*ub,rhob*vb,rhob*Eb];
f_out1cell = [rhob*ub,rhob*ub^2+pb,rhob*ub*vb,rhob*Hb*ub];
g_out1cell = [rhob*vb,rhob*ub*vb,rhob*vb^2+pb,rhob*Hb*vb];


%% Output formatting

%create empty output matrixes in the right format
q_out = NaN(2,1,4);
f_out = q_out;
g_out = q_out;

%convert the gc into the right format for output
q_out(1,:,:) = q_out1cell;
q_out(2,:,:) = q_out1cell;
f_out(1,:,:) = f_out1cell;
f_out(2,:,:) = f_out1cell;
g_out(1,:,:) = g_out1cell;
g_out(2,:,:) = g_out1cell;


end

