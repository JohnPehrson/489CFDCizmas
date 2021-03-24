function [Disout] = Dis(i,j,x,y,q,p,c,v2,v4)
%The dissipation function finds the dissipation at a single cell using a
%bunch of sub-functions
%calls switches, edge lengths, eigenvalues, and q

%uses: 
    %i and j, cel locations
    %x and y, node locaitons (which implicitly node-cell assumptions)
    %q, 3d matrix
    %p, a 2d matrix of static pressures *** should change ***
    %c, gas speed of sound
    %v2 and v4, dissipation constants (like viscosity things?)
    

[length] = edgelength(i,j,x,y,face)
[switch2out] = switch2(i,j,v2,p,type)
[switch4out] = switch4(i,j,v4,p,type,pm)
[eigenvalueout] = eigenvalueDis(c,i,j,x,y,q,face)

s2east = 1/2*(switch2(i,j,v2,p,'xi')+switch2(i+1,j,v2,p,'xi'))*edgelength(i,j,x,y,'E')*eigenvalueDis(c,i,j,x,y,q,'E')*(q(i+1,j,:)-q(i,j,:));
s2west = 1/2*(switch2(i,j,v2,p,'xi')+switch2(i-1,j,v2,p,'xi'))*edgelength(i,j,x,y,'E')*eigenvalueDis(c,i,j,x,y,q,'E')*(q(i,j,:)-q(i-1,j,:));
s2north = 
s2south = 
s4east = 
s4west = 
s4north = 
s4south = 

Disout = (s2east+ s2west+ s2north+ s2south+ s4east+ s4west+ s4north+ s4south);




end

