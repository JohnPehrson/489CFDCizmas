function [switch2out] = switch2(v2,p5by5,type,face)
%This is the switching term associated with second order dissipation. 
%type should either be 'xi' or 'n', related to which direction

%The input vector p refers to the range of pressures in the nearby cells
    %P5by5 should be a 5x5 matrix of pressures near the cell
    %p should be a 1x5 vector for xi, and a 5x1 vector for n
    %if the switch is for s(2)(i+1,j), the range of input pressures would
        %be p(i+2,j),p(i+1,j),p(i,j),p(i-1,j), p(i-2,j)
    %P(1) should be the largest value, such as P(i+2,j) for a centered case

    i = 3;
    j = 3;
    p = p5by5;
    
    switch type
        case 'xi'
            %calculate s(2)xi
            s2xicentral = v2*abs(p(i+1,j)-2*p(i,j)+p(i-1,j))/(p(i+1,j)+2*p(i,j)+p(i-1,j));
            switch face
                case 'E'
                    i = i+1;
                    s2xiadj = v2*abs(p(i+1,j)-2*p(i,j)+p(i-1,j))/(p(i+1,j)+2*p(i,j)+p(i-1,j));
                    switch2out = 0.5*(s2xicentral+s2xiadj);
                case 'W'
                    i = i-1;
                    s2xiadj = v2*abs(p(i+1,j)-2*p(i,j)+p(i-1,j))/(p(i+1,j)+2*p(i,j)+p(i-1,j));
                    switch2out = 0.5*(s2xicentral+s2xiadj);
            end
        case 'n'
            %calculate s(2)n
            s2ncentral = v2*abs(p(i,j+1)-2*p(i,j)+p(i-1,j))/(p(i+1,j)+2*p(i,j)+p(i-1,j));
            switch face
                case 'N'
                    j = j+1;
                    s2nadj = v2*abs(p(i,j+1)-2*p(i,j)+p(i-1,j))/(p(i+1,j)+2*p(i,j)+p(i-1,j));
                    switch2out = 0.5*(s2ncentral+s2nadj);
                case 'S'
                    j = j-1;
                    s2nadj = v2*abs(p(i,j+1)-2*p(i,j)+p(i-1,j))/(p(i+1,j)+2*p(i,j)+p(i-1,j));
                    switch2out = 0.5*(s2ncentral+s2nadj);
            end
    end
    
    
end

