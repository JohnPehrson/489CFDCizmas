function [AoAtangent] = bottomwalltangent(x)
%This function will give the slope of the tangent to the wall at a distance
%along the bottom wall x. This is used to find the direction of the flow
%for the IC

%AoAtangent is reported in radians

    if x>=0 && x<=2
        AoAtangent = 0;
    elseif x>2 && x<3
        %y = 0.1*sin((x-2)*pi);
        AoAtangent = 0.1*pi*cos((x-2)*pi);
    elseif x>=3 && x<=5 
        AoAtangent = 0;
    else
        fprintf('x is outside the bounds of the geometry');
    end


end

