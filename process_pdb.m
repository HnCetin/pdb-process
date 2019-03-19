function [Protein] = process_pdb(pdbstruct, chain)
% Syntax    : [Protein] = process_pdb(pdbstruct, chain)
%
% pdbstruct : A struct that obtained from getpdb() function
% chain     : Chain ID in the protin (for example: 'A')

Protein.Description = char([pdbstruct.Header.idCode, ' (Chain ', chain, ')' ]);
Protein.XYZ = backbone_cordinates(pdbstruct, chain);
Protein.BondLengths = bond_lengths(Protein.XYZ);
Protein.BondAngles = bond_angles(Protein.XYZ);
Protein.TorsionalAngles = torsional_angle(Protein.XYZ);
Protein.PhiPsiAngles(:,1) = Protein.TorsionalAngles(1:2:end,:);
Protein.PhiPsiAngles(1:length(Protein.PhiPsiAngles(:,1))-1,2) = Protein.TorsionalAngles(2:2:end,:);


CA_Cords = Protein.XYZ(2:3:end,:);
Protein.ReducedModel.XYZ = CA_Cords;
Protein.ReducedModel.BondLengths = bond_lengths(CA_Cords);
Protein.ReducedModel.BondAngles = bond_angles(CA_Cords);
Protein.ReducedModel.TorsionalAngles = torsional_angle(CA_Cords);
Protein.ReducedModel.RVector = CA_Cords(end,:) - CA_Cords(1,:);
Protein.ReducedModel.RadiusG = radius_gyration(CA_Cords);


Reversed_Cords =  transform_to_xyz(Protein.ReducedModel.BondLengths, ...
                                   Protein.ReducedModel.BondAngles,...
                                   Protein.ReducedModel.TorsionalAngles);
Protein.ReversedModel.XYZ = Reversed_Cords;
Protein.ReversedModel.BondLengths = bond_lengths(Reversed_Cords);
Protein.ReversedModel.BondAngles = bond_angles(Reversed_Cords);
Protein.ReversedModel.TorsionalAngles = torsional_angle(Reversed_Cords);
Protein.ReversedModel.RVector = Reversed_Cords(end,:) - Reversed_Cords(1,:);
Protein.ReversedModel.RadiusG = radius_gyration(Reversed_Cords);


end

function [coords] = backbone_cordinates(pdbstruct, chain)
% Syntax    : backbone_cordinates(Protein, chain)
%
% pdbstruct : A struct that obtained from getpdb() function
% chain     : Chain ID in the protin (for example: 'A')
j=1;
    for i=1:length(pdbstruct.Model(1).Atom)
        if (pdbstruct.Model(1).Atom(i).chainID == chain)
            a = pdbstruct.Model(1).Atom(i).AtomName;
            if any(strcmp(a, {'N','CA','C'}))
                coords(j,1) = pdbstruct.Model(1).Atom(i).X;
                coords(j,2) = pdbstruct.Model(1).Atom(i).Y;
                coords(j,3) = pdbstruct.Model(1).Atom(i).Z;
                j = j + 1;
            end
        end
    end

end

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

function [radius_gyration] = radius_gyration(coords)
% Syntax: [radius_gyration] = radius_gyration(coords)
%
% coords: A double structure contains values for X, Y and Z coordinates.

    masscenter(1,1) = mean(coords(:,1)); 
    masscenter(1,2) = mean(coords(:,2));
    masscenter(1,3) = mean(coords(:,3));
    masscenter_cords = [coords(:,1) - masscenter(1), coords(:,2) - masscenter(2), coords(:,3) - masscenter(3)];
    clear pro2_mc
    radius_gyration = sqrt(mean(diag(masscenter_cords*masscenter_cords'))); 
    
end

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
