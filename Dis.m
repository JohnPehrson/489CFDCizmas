function [Disout] = Dis(v2,v4,user_Gamma,x_abcd,y_abcd,q5by5,p5by5,ei3by3)
%The dissipation function finds the dissipation at a single cell using a
%bunch of sub-functions
%calls switches, edge lengths, eigenvalues, and q

%Inputs:    %v2 and v4 are dissipation constants
            %c is speed of the gas
            %x_abcd and y_abcd are the locations of the nodes bounding the
                %cell
            %q5by5 is a 5x5x4 matrix that is centered on the cell 
            %p5by5 is a 5x5x4 matrix that is centered on the cell
            %eig3by3 is a 3x3x4 matrix of the eigenvalues centered at cell
                %eig(2,2,1) is lambda_n N, eig(2,2,2) is lambda_xi E, 3 is
                %S, 4 is W
            

%uses:   [switch2out] = switch2(v2,p5by5,type,face)
            %p5by5 should be a 5x5 matrix of pressures near the cell
        %[switch4out] = switch4(v2,v4,p5by5,type,face)
            %p5by5 should be a 5x5 vector of pressures near the cell
        %[length] = edgelength(x_abcd,y_abcd,face)
        %[q_op_vec_out] = q_opDis(q5by5,type,order,face)
            %a 5x5x4 q matrix with only + terms being nonzero, centered at the cell being
            %calculated in the dissipation function

%define eigenvalues on the faces
i = 2; 
j = 2;
lambda_xi_E = 0.5*(ei3by3(i,j,2)+ei3by3(i+1,j,4));
lambda_xi_W = 0.5*(ei3by3(i,j,4)+ei3by3(i-1,j,2));
lambda_n_N = 0.5*(ei3by3(i,j,1)+ei3by3(i,j+1,3));
lambda_n_S = 0.5*(ei3by3(i,j,3)+ei3by3(i,j-1,1));


%define individual terms in the dissipation
s2east = switch2(v2,p5by5,'xi','E')*edgelength(x_abcd,y_abcd,'E')*lambda_xi_E*q_opDis(q5by5,'xi',1,'E');
s2west = switch2(v2,p5by5,'xi','W')*edgelength(x_abcd,y_abcd,'W')*lambda_xi_W*q_opDis(q5by5,'xi',1,'W');
s2north = switch2(v2,p5by5,'n','N')*edgelength(x_abcd,y_abcd,'N')*lambda_n_N*q_opDis(q5by5,'n',1,'N');
s2south = switch2(v2,p5by5,'n','S')*edgelength(x_abcd,y_abcd,'S')*lambda_n_S*q_opDis(q5by5,'n',1,'S');
s4east = switch4(v2,v4,p5by5,'xi','E')*edgelength(x_abcd,y_abcd,'E')*lambda_xi_E*q_opDis(q5by5,'xi',3,'E');
s4west = switch4(v2,v4,p5by5,'xi','W')*edgelength(x_abcd,y_abcd,'W')*lambda_xi_W*q_opDis(q5by5,'xi',3,'W');
s4north = switch4(v2,v4,p5by5,'n','N')*edgelength(x_abcd,y_abcd,'N')*lambda_n_N*q_opDis(q5by5,'n',3,'N');
s4south = switch4(v2,v4,p5by5,'n','S')*edgelength(x_abcd,y_abcd,'S')*lambda_n_S*q_opDis(q5by5,'n',3,'S');
                 
Disout = (s2east-s2west+s2north-s2south)-(s4east-s4west+s4north-s4south);
Disout = squeeze(Disout);
end

