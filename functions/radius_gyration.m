function [radius_gyration] = radius_gyration(coords)

    masscenter(1,1) = mean(coords(:,1)); 
    masscenter(1,2) = mean(coords(:,2));
    masscenter(1,3) = mean(coords(:,3));
    masscenter_cords = [coords(:,1) - masscenter(1), coords(:,2) - masscenter(2), coords(:,3) - masscenter(3)];
    clear pro2_mc
    radius_gyration = sqrt(mean(diag(masscenter_cords*masscenter_cords'))); 
    
end