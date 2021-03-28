function [Residual] = findResidual(x_abcd,y_abcd,f_cNESW,g_cNESW)
%Calculates the residual R(i,j) of a single cell
    %Required inputs: f_cNESW where f(1) = f(i,j), f(2) = f(i,j+1), f(3) = f(i+1,j), vector of length 5 f(i+-1,j+-1)
                    % g_cNESW where g(1) = g(i,j), g(2) = g(i,j+1), g(3) = g(i+1,j), vector of length 5 g(i+-1,j+-1)
                    % x_abcd where x(1) is x_a, x(2) is x_b, etc.
                    % y_abcd where y(1) is y_a, y(2) is y_b, etc.
%Calculate fluxes on the boundaries
fN = 0.5*(f(1)+f(2));
fE = 0.5*(f(1)+f(3));
fS = 0.5*(f(1)+f(4));
fW = 0.5*(f(1)+f(5));
gN = 0.5*(g(1)+g(2));
gE = 0.5*(g(1)+g(3));
gS = 0.5*(g(1)+g(4));
gW = 0.5*(g(1)+g(5));
%Calculate the Residual        
Residual = fN*(y(4)-y(3))-gN*(x(4)-x(3))+fW*(y(1)-y(4))-gW*(x(1)-x(4))+...
           fS*(y(2)-y(1))-gS*(x(2)-x(1))+fE*(y(3)-y(2))-gE*(x(3)-x(2));
end

