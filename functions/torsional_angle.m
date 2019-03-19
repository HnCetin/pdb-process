function [t_angles] = torsional_angle(cords)
% Syntax: [t_angles] = torsional_angle(cords)
%
% coords: A double structure contains values for X, Y and Z coordinates.

    for i=4:length(cords)

        bc = (cords(i-1,:) - cords(i-2,:)) / norm(cords(i-1,:) - cords(i-2,:));  

        % Orthagonal normalized vectors
        ab_orth  = (cords(i-3,:) - cords(i-2,:) - bc * dot(cords(i-3,:)...
                    - cords(i-2,:), bc)) / norm(cords(i-3,:) ...
                    - cords(i-2,:) - bc * dot(cords(i-3,:) - cords(i-2,:), bc));
        cd_orth  = (cords(i,:) - cords(i-1,:) - bc * dot(cords(i,:) - cords(i-1,:), bc)) ...
                    / norm(cords(i,:) - cords(i-1,:) - bc * dot(cords(i,:) - cords(i-1,:), bc));

        torsional_angle = acos(dot(ab_orth, cd_orth)) * 180 / pi;
        
        sign = dot(cross(ab_orth, cd_orth), bc); 
            if (sign < 0)
                torsional_angle = -torsional_angle;
            end 
            
        t_angles(i-2,:) = torsional_angle; % first bond does not have torsional angle
    end
    t_angles(end+1,1) = 0; % last bond does not have torsional angle
end
