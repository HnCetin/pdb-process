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