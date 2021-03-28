function [Disout] = Dis(v2,v4,c,x_abcd,y_abcd,q5by5,p5by5)
%The dissipation function finds the dissipation at a single cell using a
%bunch of sub-functions
%calls switches, edge lengths, eigenvalues, and q

%Inputs:    %v2 and v4 are dissipation constants
            %c is speed of the gas
            %x_abcd and y_abcd are the locations of the nodes bounding the
                %cell
            %q5by5 is a 5x5x4 matrix that is centered on the cell 
            %p5by5 is a 5x5x4 matrix that is centered on the cell

%uses:   [switch2out] = switch2(v2,p5by5,type,face)
            %p5by5 should be a 5x5 matrix of pressures near the cell
        %[switch4out] = switch4(v2,v4,p5by5,type,face)
            %p5by5 should be a 5x5 vector of pressures near the cell
        %[length] = edgelength(x_abcd,y_abcd,face)
        %[eigenvalueout] = eigenvalueDis(c,x_abcd,y_abcd,q3by3,face)
            %q is a 3x3x4 matrix that gives the q values of the central cell and the 4
            %bordering cells,
        %[q_op_vec_out] = q_opDis(q5by5,type,order,face)
            %a 5x5x4 q matrix with only + terms being nonzero, centered at the cell being
            %calculated in the dissipation function

%define subinput
q3by3 = q5by5(2:4,2:4,:);

%define individual terms in the dissipation
s2east = switch2(v2,p5by5,'xi','E')*edgelength(x_abcd,y_abcd,'E')*eigenvalueDis(c,x_abcd,y_abcd,q3by3,'E')*q_opDis(q5by5,'xi',1,'E');
s2west = switch2(v2,p5by5,'xi','W')*edgelength(x_abcd,y_abcd,'W')*eigenvalueDis(c,x_abcd,y_abcd,q3by3,'W')*q_opDis(q5by5,'xi',1,'W');
s2north = switch2(v2,p5by5,'n','N')*edgelength(x_abcd,y_abcd,'N')*eigenvalueDis(c,x_abcd,y_abcd,q3by3,'N')*q_opDis(q5by5,'n',1,'N');
s2south = switch2(v2,p5by5,'n','S')*edgelength(x_abcd,y_abcd,'S')*eigenvalueDis(c,x_abcd,y_abcd,q3by3,'S')*q_opDis(q5by5,'n',1,'S');
s4east = switch4(v2,v4,p5by5,'xi','E')*edgelength(x_abcd,y_abcd,'E')*eigenvalueDis(c,x_abcd,y_abcd,q3by3,'E')*q_opDis(q5by5,'xi',3,'E');
s4west = switch4(v2,v4,p5by5,'xi','W')*edgelength(x_abcd,y_abcd,'W')*eigenvalueDis(c,x_abcd,y_abcd,q3by3,'W')*q_opDis(q5by5,'xi',3,'W');
s4north = switch4(v2,v4,p5by5,'n','N')*edgelength(x_abcd,y_abcd,'N')*eigenvalueDis(c,x_abcd,y_abcd,q3by3,'N')*q_opDis(q5by5,'n',3,'N');
s4south = switch2(v2,v4,p5by5,'n','S')*edgelength(x_abcd,y_abcd,'S')*eigenvalueDis(c,x_abcd,y_abcd,q3by3,'S')*q_opDis(q5by5,'n',3,'S');

Disout = (s2east-s2west+s2north-s2south)-(s4east-s4west+s4north-s4south);
end

