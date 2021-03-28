function [switch4out] = switch4(v2,v4,p5by5,type,face)
%This is the switching term associated with fourth order dissipation. 
%type should either be xi or n, related to which direction we are watching

%The input vector p refers to the range of pressures in the nearby cells
    %p5by5 should be a 5x5 vector 
    % p uses p5by5 to make a 1x5 matrix for xi, and a 5x1 vector for n
    %if the switch is for s(4)(i+1,j), the range of input pressures would
        %be p(i+2,j),p(i+1,j),p(i,j),p(i-1,j),p(i-2,j)
            %longer than s(2) because it calls two different switches, and
            %needs to be able to call 3 different switches with the pm
    %P(1) should be the largest value,  P(i+2,j) or P(i,j+2) 

%face denotes the face on which the switch is being activated in the
    %dissipation function

%     switch type
%         case 'xi'         
%             switch face
%                 case 'E'
%                     switch4out = max(0,v4-1/2*(switch2(v2,p5by5,type,face)));
%                 case 'W'
%                     switch4out = max(0,v4-1/2*(switch2(v2,p5by5,type,face)));
%             end
%             
%         case 'n'
%             switch face
%                 case 'N'
%                     switch4out = max(0,v4-1/2*(switch2(v2,p5by5,type,face)));
%                 case 'S'
%                     switch4out = max(0,v4-1/2*(switch2(v2,p5by5,type,face)));
%             end
%     end

    switch4out = max(0,v4-1/2*(switch2(v2,p5by5,type,face)));
    
end

