function [switch2out] = switch2(i,j,v2,p,type)
%This is the switching term associated with second order dissipation. 
%type should either be xi or n, related to which direction we are watching

%currently i,j, and p all refer to the entire matrix of p (which is . i and j denote
%cells, not nodes

%P IS A MATRIX OF ONLY THE STATIC PRESSURE, NEEDS TO BE FIXED! THIS ISN'T
%THE MATRIX q,f, or g

    switch type
        case 'xi'
        switch2out = v2*(abs((p(i+1,j)-2*p(i,j)+p(i-1,j)))/(p(i+1,j)+2*p(i,j)+p(i-1,j)));
        case 'n'
        switch2out = v2*(abs((p(i,j+1)-2*p(i,j)+p(i,j-1)))/(p(i,j+1)+2*p(i,j)+p(i,j-1)));
    end

end

