function [xyz] = transform_to_xyz(l, theta, phi)
% Syntax: [xyz] = transform_to_xyz(l, theta, phi)
%
% l     : lenght of bonds
% theta : angle between two bonds
% phi   : dihedral angle between three bonds
%
% Dimensions must be the same.

        % Transformation matrix
        for i=1:length(theta)
        tm{i} = [cos(theta(i)),             sin(theta(i)),             0; ...
                 sin(theta(i))*cos(phi(i)), -cos(theta(i))*cos(phi(i)), sin(phi(i)) ; ...
                 sin(theta(i))*sin(phi(i)), -cos(theta(i))*sin(phi(i)), -cos(phi(i)) ];
        end

        % Iterative tms are collected
        tm_iter{1} = tm{1};
        for i=1:length(tm)-1
            tm_iter{i+1} =  tm_iter{i} * tm{i+1};
        end

        % Forming coordinates
        xyz(1,:) = [0; 0; 0];
        for i=1:length(theta)
            xyz(i+1,:) = xyz(i,:) + transpose(tm_iter{i} * [l(i); 0; 0]);
        end 

end