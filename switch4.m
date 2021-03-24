function [switch4out] = switch4(i,j,v4,p,type,pm)
%This is the switching term associated with fourth order dissipation. 
%type should either be xi or n, related to which direction we are watching

%currently i,j, and p all refer to the entire matrix of p. i and j denote
%cells, not nodes

%variable pm denotes if the switch is +1/2 or -1/2, becasue we will be
%dealing with integers and that is an internal kind of transformation to
%nodes
    switch type
        case 'xi'
            switch pm
                case '+'
                    switch4out = max(0,v4-1/2*(switch2(i,j,v2,p,type)+switch2(i+1,j,v2,p,type)));
                case '-'
                    switch4out = max(0,v4-1/2*(switch2(i,j,v2,p,type)+switch2(i-1,j,v2,p,type)));
            end
            
        case 'n'
            switch pm
                case '+'
                    switch4out = max(0,v4-1/2*(switch2(i,j,v2,p,type)+switch2(i,j+1,v2,p,type)));
                case '-'
                    switch4out = max(0,v4-1/2*(switch2(i,j,v2,p,type)+switch2(i,j-1,v2,p,type)));
            end
    end
end

