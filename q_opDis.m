function [q_op_vec_out] = q_opDis(q5by5,type,order,face)
%This function is a subroutine that calculates the delta_xi or delta_n
%functions in the dissipation calculation.

q = q5by5;
%Inputs
%a 5x5x4 q matrix centered at the cell being
    %calculated in the dissipation function
%type = 'xi' or 'n'
%order = 1 or 3
%face = 'N','E','W', or 'S'

%place reference points in the middle of the 5x5 matrix for ease of
%interpreation
i = 3;
j = 3;
q_op_vec_out = NaN;

    switch order
        case 1 %operatior of type 1
            switch type
                case 'xi' 
                    switch face
                        case 'E'
                            q_op_vec_out = q(i+1,j,:)-q(i,j,:);
                        case 'W'
                            q_op_vec_out = q(i,j,:)-q(i-1,j,:);
                    end
                case 'n'
                    switch face
                        case 'N'
                            q_op_vec_out = q(i,j+1,:)-q(i,j,:);
                        case 'S'
                            q_op_vec_out = q(i,j,:)-q(i,j-1,:);
                    end
            end
        case 3 %operator of type 3
            switch type
                case 'xi'
                    switch face
                        case 'E'
                            q_op_vec_out = (q(i+2,j,:)-2*q(i+1,j,:)+q(i,j,:))-(q(i+1,j,:)-2*q(i,j,:)+q(i-1,j,:));
                        case 'W'
                            q_op_vec_out = (q(i+1,j,:)-2*q(i,j,:)+q(i-1,j,:))-(q(i,j,:)-2*q(i-1,j,:)+q(i-2,j,:));
                    end
                case 'n'
                    switch face
                        case 'N'
                            q_op_vec_out = (q(i,j+2,:)-2*q(i,j+1,:)+q(i,j,:))-(q(i,j+1,:)-2*q(i,j,:)+q(i,j-1,:));
                        case 'S'
                            q_op_vec_out = (q(i,j+1,:)-2*q(i,j,:)+q(i,j-1,:))-(q(i,j,:)-2*q(i,j-1,:)+q(i,j-2,:));
                    end
            end 
    end


end

