function [Ptotal] = pressurefinder(user_Mach,user_Gamma,Pstatic)
%This function calculates the static pressure in a subsonic compressible
%flow given a resevoir (stagnation) pressure, local mach number, and gamma

%This function returns the static pressure at the mach number

Ptotal = Pstatic*(1+((user_Gamma-1)/2)*user_Mach^2)^((user_Gamma)/(user_Gamma-1));

end

