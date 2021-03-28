function [switch2out] = switch2(v2,p5by5,type,face)
%This is the switching term associated with second order dissipation. 
%type should either be 'xi' or 'n', related to which direction

%The input vector p refers to the range of pressures in the nearby cells
    %P5by5 should be a 5x5 matrix of pressures near the cell
    %p should be a 1x5 vector for xi, and a 5x1 vector for n
    %if the switch is for s(2)(i+1,j), the range of input pressures would
        %be p(i+2,j),p(i+1,j),p(i,j),p(i-1,j), p(i-2,j)
    %P(1) should be the largest value, such as P(i+2,j) for a centered case

    switch type
        case 'xi'
            p = p5by5(3,:);
            switch face
                case 'E'
                	switch2out = 0.5*(v2*(abs((p(2)-2*p(3)+p(4)))/(p(2)+2*p(3)+p(4)))+v2*(abs((p(3)-2*p(4)+p(5)))/(p(3)+2*p(4)+p(5))));
                case 'W'
                    switch2out = 0.5*(v2*(abs((p(2)-2*p(3)+p(4)))/(p(2)+2*p(3)+p(4)))+v2*(abs((p(1)-2*p(2)+p(3)))/(p(1)+2*p(2)+p(3))));
            end
        case 'n'
            p = p5by5(:,3);
            switch face
                case 'N'
                    switch2out = 0.5*(v2*(abs((p(2)-2*p(3)+p(4)))/(p(2)+2*p(3)+p(4)))+v2*(abs((p(3)-2*p(4)+p(5)))/(p(3)+2*p(4)+p(5))));
                case 'S'
                    switch2out = 0.5*(v2*(abs((p(2)-2*p(3)+p(4)))/(p(2)+2*p(3)+p(4)))+v2*(abs((p(1)-2*p(2)+p(3)))/(p(1)+2*p(2)+p(3))));
            end
    end
    
    
    
    
    switch type
        case 'xi'
        p = p3by3(2,:);
        switch2out = v2*(abs((p(1)-2*p(2)+p(3)))/(p(1)+2*p(2)+p(3)));
        case 'n'
        p = p3by3(:,2);
        switch2out = v2*(abs((p(1)-2*p(2)+p(3)))/(p(1)+2*p(2)+p(3)));
    end

end

