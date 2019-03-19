function [bond_lengths] = bond_lengths(coords)
% Syntax: [bond_lengths] = bond_lengths(coords)
%
% coords: A double structure contains values for X, Y and Z coordinates.
    for i=2:length(coords)
        bond_lengths(i-1,:) = sqrt( (coords(i,1) - coords(i-1,1)).^2 + ...
                                    (coords(i,2) - coords(i-1,2)).^2 + ...
                                    (coords(i,3) - coords(i-1,3)).^2 );
    end
end