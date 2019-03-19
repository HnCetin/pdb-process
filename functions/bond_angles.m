function [angles] = bond_angles(coords) 
% Syntax: [angles] = bond_angles(coords) 
%
% coords: A double structure contains values for X, Y and Z coordinates.

    for i=3:length(coords)

a = coords(i-2,:) - coords(i-1,:);
b = coords(i-1,:) - coords(i,:);
angles(i-1,:) = atan2d(norm(cross(a,b)),dot(a,b)); % first bond does not have an angle                          
                   
    end
end     