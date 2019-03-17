function [anglesdegree] = bond_angles(coords) 
    for i=2:length(coords)
        anglesdegree(i-1,:) = atan2d(norm(cross(coords(i-1,:),coords(i,:))), ...
                                     dot(coords(i-1,:),coords(i,:)));
    end
end     