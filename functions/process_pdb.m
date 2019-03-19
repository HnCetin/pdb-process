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