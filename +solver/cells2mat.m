function [x,y,z] = cells2mat(X) % X is now a 1 x T cell array with data per time unit
%Process velocities

if size(X{1}, 2) == 3 % vector quantity
    cm = cell2mat(X);
    
    x = cm(:,1:3:end);
    y = cm(:,2:3:end);
    z = cm(:,3:3:end);
else
    x = cell2mat(X);
    y = NaN;
    z = NaN;
end

end