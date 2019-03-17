function [coords] = backbone_cordinates(Protein, chain)
j=1;
    for i=1:length(Protein.Model(1).Atom)
        if (Protein.Model(1).Atom(i).chainID == chain)
            a = Protein.Model(1).Atom(i).AtomName;
            if any(strcmp(a, {'N','CA','C'}))
                coords(j,1) = Protein.Model(1).Atom(i).X;
                coords(j,2) = Protein.Model(1).Atom(i).Y;
                coords(j,3) = Protein.Model(1).Atom(i).Z;
                j = j + 1;
            end
        end
    end

end