function [Protein] = process_pdb(pdbstruct, chain)
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


end