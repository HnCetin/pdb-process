%% Setup
Protein1 = getpdb('4AKE','ToFile','4AKE.pdb');
Protein2 = getpdb('1AKE','ToFile','1AKE.pdb');

%% -- All Atom Chain Model : Backbone atoms
%% Coordinates 
% Coordinates of the backbone (A chain only) is extracted

j = 1;
for i=1:length(Protein1.Model(1).Atom)
    if (Protein1.Model(1).Atom(i).chainID == 'A')
        a = Protein1.Model(1).Atom(i).AtomName;
        if any(strcmp(a, {'N','CA','C'}))
            pro1_cords(j,1) = Protein1.Model(1).Atom(i).X;
            pro1_cords(j,2) = Protein1.Model(1).Atom(i).Y;
            pro1_cords(j,3) = Protein1.Model(1).Atom(i).Z;
            j = j + 1;
        end
    end
end
j = 1;
for i=1:length(Protein2.Model(1).Atom)
    if (Protein2.Model(1).Atom(i).chainID == 'A')
        a = Protein2.Model(1).Atom(i).AtomName;
        if any(strcmp(a, {'N','CA','C'}))
            pro2_cords(j,1) = Protein2.Model(1).Atom(i).X;
            pro2_cords(j,2) = Protein2.Model(1).Atom(i).Y;
            pro2_cords(j,3) = Protein2.Model(1).Atom(i).Z;
            j = j + 1;
        end
    end
end
clear i j a

%% Bond lengths 
for i=2:length(pro1_cords)
    pro1_lengths(i-1,:) = sqrt( (pro1_cords(i,1) - pro1_cords(i-1,1)).^2 + ...
                            (pro1_cords(i,2) - pro1_cords(i-1,2)).^2 + ...
                            (pro1_cords(i,3) - pro1_cords(i-1,3)).^2 );
end
for i=2:length(pro2_cords)
    pro2_lengths(i-1,:) = sqrt( (pro2_cords(i,1) - pro2_cords(i-1,1)).^2 + ...
                            (pro2_cords(i,2) - pro2_cords(i-1,2)).^2 + ...
                            (pro2_cords(i,3) - pro2_cords(i-1,3)).^2 );
end
clear i
%% Bond angles
for i=2:length(pro1_cords)
    pro1_angles(i-1,:) = atan2d(norm(cross(pro1_cords(i-1,:),pro1_cords(i,:))), ...
                                 dot(pro1_cords(i-1,:),pro1_cords(i,:)));
end
for i=2:length(pro2_cords)
    pro2_angles(i-1,:) = atan2d(norm(cross(pro2_cords(i-1,:),pro2_cords(i,:))), ...
                                 dot(pro2_cords(i-1,:),pro2_cords(i,:)));
end
clear i

%% Torsional angles
for i=4:length(pro1_cords)

        bc = (pro1_cords(i-1,:) - pro1_cords(i-2,:)) / norm(pro1_cords(i-1,:) - pro1_cords(i-2,:));  

        % Orthagonal normalized vectors
        ab_orth  = (pro1_cords(i-3,:) - pro1_cords(i-2,:) - bc * dot(pro1_cords(i-3,:)...
                    - pro1_cords(i-2,:), bc)) / norm(pro1_cords(i-3,:) ...
                    - pro1_cords(i-2,:) - bc * dot(pro1_cords(i-3,:) - pro1_cords(i-2,:), bc));
        cd_orth  = (pro1_cords(i,:) - pro1_cords(i-1,:) - bc * dot(pro1_cords(i,:) - pro1_cords(i-1,:), bc)) ...
                    / norm(pro1_cords(i,:) - pro1_cords(i-1,:) - bc * dot(pro1_cords(i,:) - pro1_cords(i-1,:), bc));

        torsional_angle = acos(dot(ab_orth, cd_orth)) * 180 / pi;
        sign = dot(cross(ab_orth, cd_orth), bc); 
            if (sign < 0)
                torsional_angle = -torsional_angle;

            end  
         pro1_tors(i-3,:) = torsional_angle;
  end
 clear i bc ab_orth cd_orth sign torsional_angle
 
for i=4:length(pro2_cords)

        bc = (pro2_cords(i-1,:) - pro2_cords(i-2,:)) / norm(pro2_cords(i-1,:) - pro2_cords(i-2,:));  

        % Orthagonal normalized vectors
        ab_orth  = (pro2_cords(i-3,:) - pro2_cords(i-2,:) - bc * dot(pro2_cords(i-3,:)...
                    - pro2_cords(i-2,:), bc)) / norm(pro2_cords(i-3,:) ...
                    - pro2_cords(i-2,:) - bc * dot(pro2_cords(i-3,:) - pro2_cords(i-2,:), bc));
        cd_orth  = (pro2_cords(i,:) - pro2_cords(i-1,:) - bc * dot(pro2_cords(i,:) - pro2_cords(i-1,:), bc)) ...
                    / norm(pro2_cords(i,:) - pro2_cords(i-1,:) - bc * dot(pro2_cords(i,:) - pro2_cords(i-1,:), bc));

        torsional_angle = acos(dot(ab_orth, cd_orth)) * 180 / pi;
        sign = dot(cross(ab_orth, cd_orth), bc); 
            if (sign < 0)
                torsional_angle = -torsional_angle;

            end  
         pro2_tors(i-3,:) = torsional_angle;
  end
 clear i bc ab_orth cd_orth sign torsional_angle
 

% Psi-Phi Angles
pro1_tors_phi = pro1_tors(1:2:end-1,:);
pro1_tors_psi = pro1_tors(2:2:end,:);

pro2_phi = pro2_tors(1:2:end-1,:);
pro2_psi = pro2_tors(2:2:end,:);


%% -- Reduced Models : CA atoms only

%% Coordinates
pro1_ca_cords = pro1_cords(2:3:end,:); % CAs are repeating every 3rd row in the backbone coordinates
pro2_ca_cords = pro2_cords(2:3:end,:); % CAs are repeating every 3rd row in the backbone coordinates

%% Bond lengths
for i=2:length(pro1_ca_cords)
    pro1_ca_lengths(i-1,:) = sqrt( (pro1_ca_cords(i,1) - pro1_ca_cords(i-1,1)).^2 + ...
                            (pro1_ca_cords(i,2) - pro1_ca_cords(i-1,2)).^2 + ...
                            (pro1_ca_cords(i,3) - pro1_ca_cords(i-1,3)).^2 );
end
for i=2:length(pro2_ca_cords)
    pro2_ca_lengths(i-1,:) = sqrt( (pro2_ca_cords(i,1) - pro2_ca_cords(i-1,1)).^2 + ...
                            (pro2_ca_cords(i,2) - pro2_ca_cords(i-1,2)).^2 + ...
                            (pro2_ca_cords(i,3) - pro2_ca_cords(i-1,3)).^2 );
end
clear i

%% Bond angles
for i=2:length(pro1_ca_cords)
    pro1_ca_angles(i-1,:) = atan2d(norm(cross(pro1_ca_cords(i-1,:),pro1_ca_cords(i,:))), ...
                                 dot(pro1_ca_cords(i-1,:),pro1_ca_cords(i,:)));
end
for i=2:length(pro2_ca_cords)
    pro2_ca_angles(i-1,:) = atan2d(norm(cross(pro2_ca_cords(i-1,:),pro2_ca_cords(i,:))), ...
                                 dot(pro2_ca_cords(i-1,:),pro2_ca_cords(i,:)));
end
clear i

%% Torsional angles
for i=4:length(pro1_ca_cords)

        bc = (pro1_ca_cords(i-1,:) - pro1_ca_cords(i-2,:)) / norm(pro1_ca_cords(i-1,:) - pro1_ca_cords(i-2,:));  

        % Orthagonal normalized vectors
        ab_orth  = (pro1_ca_cords(i-3,:) - pro1_ca_cords(i-2,:) - bc * dot(pro1_ca_cords(i-3,:)...
                    - pro1_ca_cords(i-2,:), bc)) / norm(pro1_ca_cords(i-3,:) ...
                    - pro1_ca_cords(i-2,:) - bc * dot(pro1_ca_cords(i-3,:) - pro1_ca_cords(i-2,:), bc));
        cd_orth  = (pro1_ca_cords(i,:) - pro1_ca_cords(i-1,:) - bc * dot(pro1_ca_cords(i,:) - pro1_ca_cords(i-1,:), bc)) ...
                    / norm(pro1_ca_cords(i,:) - pro1_ca_cords(i-1,:) - bc * dot(pro1_ca_cords(i,:) - pro1_ca_cords(i-1,:), bc));

        torsional_angle = acos(dot(ab_orth, cd_orth)) * 180 / pi;
        sign = dot(cross(ab_orth, cd_orth), bc); 
            if (sign < 0)
                torsional_angle = -torsional_angle;

            end  
         pro1_ca_tors(i-3,:) = torsional_angle;
end
 clear i bc ab_orth cd_orth sign torsional_angle
  
for i=4:length(pro2_ca_cords)

        bc = (pro2_ca_cords(i-1,:) - pro2_ca_cords(i-2,:)) / norm(pro2_ca_cords(i-1,:) - pro2_ca_cords(i-2,:));  

        % Orthagonal normalized vectors
        ab_orth  = (pro2_ca_cords(i-3,:) - pro2_ca_cords(i-2,:) - bc * dot(pro2_ca_cords(i-3,:)...
                    - pro2_ca_cords(i-2,:), bc)) / norm(pro2_ca_cords(i-3,:) ...
                    - pro2_ca_cords(i-2,:) - bc * dot(pro2_ca_cords(i-3,:) - pro2_ca_cords(i-2,:), bc));
        cd_orth  = (pro2_ca_cords(i,:) - pro2_ca_cords(i-1,:) - bc * dot(pro2_ca_cords(i,:) - pro2_ca_cords(i-1,:), bc)) ...
                    / norm(pro2_ca_cords(i,:) - pro2_ca_cords(i-1,:) - bc * dot(pro2_ca_cords(i,:) - pro2_ca_cords(i-1,:), bc));

        torsional_angle = acos(dot(ab_orth, cd_orth)) * 180 / pi;
        sign = dot(cross(ab_orth, cd_orth), bc); 
            if (sign < 0)
                torsional_angle = -torsional_angle;

            end  
         pro2_ca_tors(i-3,:) = torsional_angle;
  end
 clear i bc ab_orth cd_orth sign torsional_angle
            
%% R (end-to-end vector) calculation
pro1_rvector = pro1_ca_cords(end,:) - pro1_ca_cords(1,:);
pro2_rvector = pro2_ca_cords(end,:) - pro2_ca_cords(1,:);

%% S (Radius of Gyration) calculation
pro1_mc(1,1) = mean(pro1_ca_cords(:,1)); % Center of mass:
pro1_mc(1,2) = mean(pro1_ca_cords(:,2));
pro1_mc(1,3) = mean(pro1_ca_cords(:,3));
pro1_mc_cords = [pro1_ca_cords(:,1) - pro1_mc(1), pro1_ca_cords(:,2) - pro1_mc(2), pro1_ca_cords(:,3) - pro1_mc(3)];
clear pro2_mc
pro1_rgyration = sqrt(mean(diag(pro1_mc_cords*pro1_mc_cords'))); 

pro2_mc(1,1) = mean(pro2_ca_cords(:,1)); % Center of mass:
pro2_mc(1,2) = mean(pro2_ca_cords(:,2));
pro2_mc(1,3) = mean(pro2_ca_cords(:,3));
pro2_mc_cords = [pro2_ca_cords(:,1) - pro2_mc(1), pro2_ca_cords(:,2) - pro2_mc(2), pro2_ca_cords(:,3) - pro1_mc(3)];
clear pro2_mc
pro2_rgyration = sqrt(mean(diag(pro2_mc_cords*pro2_mc_cords'))); 


%% -- Plots 
% Plots of bond properties:
figure(1)
    subplot(1,3,1);
        plot3(pro1_ca_cords(:,1), pro1_ca_cords(:,2), pro1_ca_cords(:,3) )
        grid on
        title('4AKE in 3D Space')
    subplot(1,3,2);
        plot3(pro2_ca_cords(:,1), pro2_ca_cords(:,2), pro2_ca_cords(:,3) )
        grid on
        title('1AKE in 3D Space')
    subplot(1,3,3);
        plot3(xyz(:,1), xyz(:,2), xyz(:,3) )
        grid on
        title('Reversed in 3D Space')

% Ramachandran Plots
figure(2)

    subplot(1,2,1);
        scatter(pro1_tors_phi, pro1_tors_psi, 15, 'filled', 'k');
            xticks([-180:45:+180]);
            yticks([-180:45:+180]);
            xlim([-180 +180]);
            ylim([-180 +180]);
            line([0 0], [180 -180],'Color','black','LineStyle','--');
            line([-180 180], [0 0],'Color','black','LineStyle','--');
            title('Ramachandran of 4AKE Chain A');
            xlabel('Phi (Degrees)');
            ylabel('Psi (Degrees)');

    subplot(1,2,2);
        scatter(pro2_phi, pro2_psi, 15, 'filled', 'k');
            xticks([-180:45:+180]);
            yticks([-180:45:+180]);
            xlim([-180 +180]);
            ylim([-180 +180]);
            line([0 0], [180 -180],'Color','black','LineStyle','--');
            line([-180 180], [0 0],'Color','black','LineStyle','--');
            title('Ramachandran of 1AKE Chain A');
            xlabel('Phi (Degrees)');
            ylabel('Psi (Degrees)');
 
% Plotted onto built-in ramachandran to see if it is correct
% ramachandran('4AKE.pdb','chain','A') 
%   hold on
%   scatter(pro1_tors_phi, pro1_tors_psi)
% ramachandran('1AKE.pdb','chain','A')
%   hold on
%   scatter(pro2_phi, pro2_psi)

% Plots of bond properties:
figure(3)
    subplot(3,1,1); % Bond length vs. residue index
        plot(pro1_ca_lengths)
            title('Bond length vs. Residue Index')
            xticks([1:5:213]);
            xlim([0, 214]);
                hold on 
                plot(pro2_ca_lengths)
                hold off
    subplot(3,1,2); % Bond angle vs. residue index
        plot(pro1_ca_angles)
            xticks([1:5:213]);
            xlim([0, 214]);
            title('Bond angle vs. Residue Index')
                hold on 
                plot(pro2_ca_angles)
                hold off
    subplot(3,1,3); % Torsion angle vs. residue index
        plot(pro1_ca_tors)
            title('Torsional angle vs. Residue Index')
            xticks([1:5:211]);
            xlim([0, 212]);
                hold on 
                plot(pro2_ca_tors)
                hold off

% Probability Distrubitions 
figure(4)
    subplot(2,2,1); % bond lengths
        plot(linspace(min(pro1_ca_lengths), max(pro1_ca_lengths), 1000), ...
             pdf(fitdist(pro1_ca_lengths, 'kernel'), ...
             linspace(min(pro1_ca_lengths), max(pro1_ca_lengths), 1000)))
            title('Probability Distrubition of Bond Lengths')
                hold on
                plot(linspace(min(pro2_ca_lengths), max(pro2_ca_lengths), 1000), ...
                     pdf(fitdist(pro2_ca_lengths, 'kernel'), ...
                     linspace(min(pro2_ca_lengths), max(pro2_ca_lengths), 1000)))
                hold off
    subplot(2,2,2);% bond angles
        plot(linspace(-180, 180, 1000), pdf(fitdist(pro1_ca_angles, 'kernel'), linspace(-180, 180, 1000)))
            title('Probability Distrubition of Bond Angles');
            xticks([-180:45:+180]);
                hold on
                plot(linspace(-180, 180, 1000), pdf(fitdist(pro2_ca_angles, 'kernel'), linspace(-180, 180, 1000)))
                hold off
    subplot(2,2,[3,4]); % torsion angles
        plot(linspace(-180, 180, 1000), pdf(fitdist(pro1_ca_tors, 'kernel'), linspace(-180, 180, 1000)))
            title('Probability Distrubition of Torsion Angles');
            xticks([-180:45:+180]);
                hold on
                plot(linspace(-180, 180, 1000), pdf(fitdist(pro2_ca_tors, 'kernel'), linspace(-180, 180, 1000)))
                hold off
 