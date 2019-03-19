%% Protein Selection
pro1 = getpdb('4AKE');
pro2 = getpdb('1AKE');
% save proteins.mat pro1 pro2

% load proteins.mat
pro1 = process_pdb(pro1, 'A');
pro2 = process_pdb(pro2, 'A');


%% Plots 
% Plots of bond properties:
figure(1)
    subplot(2,2,1);
        plot3(pro1.ReducedModel.XYZ(:,1), pro1.ReducedModel.XYZ(:,2), pro1.ReducedModel.XYZ(:,3) )
        grid on
        title('4AKE in 3D Space')
    subplot(2,2,2);
        plot3(pro2.ReducedModel.XYZ(:,1), pro2.ReducedModel.XYZ(:,2), pro2.ReducedModel.XYZ(:,3) )
        grid on
        title('1AKE in 3D Space')
        
    subplot(2,2,3);
        plot3(pro1.ReversedModel.XYZ(:,1), pro1.ReversedModel.XYZ(:,2), pro1.ReversedModel.XYZ(:,3) )
        grid on
        title('Reverse Coordinated 4AKE')
    subplot(2,2,4);
        plot3(pro2.ReversedModel.XYZ(:,1), pro2.ReversedModel.XYZ(:,2), pro2.ReversedModel.XYZ(:,3) )
        grid on
        title('Reverse Coordinated 1AKE')
        
% Ramachandran Plots
figure(2)

    subplot(1,2,1);
        scatter(pro1.PhiPsiAngles(:,1), pro1.PhiPsiAngles(:,2), 20, 'filled', 'k');
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
        scatter(pro2.PhiPsiAngles(:,1), pro2.PhiPsiAngles(:,2), 20, 'filled', 'k');
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
% hold on
% scatter(pro1.PhiPsiAngles(:,1), pro1.PhiPsiAngles(:,2))
%     
% ramachandran('1AKE.pdb','chain','A')
% hold on
% scatter(pro2.PhiPsiAngles(:,1), pro2.PhiPsiAngles(:,2))


% Plots of bond properties:
figure(3)
    subplot(3,1,1); % Bond length vs. residue index
        plot(pro1.ReducedModel.BondLengths)
            title('Bond Lengths vs. Residue Index')
            xticks([0:5:215]);
            xlim([-2, 215]);
                hold on 
                plot(pro2.ReducedModel.BondLengths)
                hold off
    subplot(3,1,2); % Bond angle vs. residue index
        plot(pro1.ReducedModel.BondAngles)
            xticks([0:5:215]);
            xlim([-2, 215]);
            title('Bond Angles vs. Residue Index')
                hold on 
                plot(pro2.ReducedModel.BondAngles)
                hold off
    subplot(3,1,3); % Torsion angle vs. residue index
        plot(pro1.ReducedModel.TorsionalAngles)
            title('Torsional Angles vs. Residue Index')
            xticks([0:5:215]);
            xlim([-2, 215]);
                hold on 
                plot(pro2.ReducedModel.TorsionalAngles)
                hold off

% Probability Distrubitions 
figure(4)
    subplot(2,2,1); % bond lengths
        plot(linspace(min(pro1.ReducedModel.BondLengths), max(pro1.ReducedModel.BondLengths), 1000), ...
             pdf(fitdist(pro1.ReducedModel.BondLengths, 'kernel'), ...
             linspace(min(pro1.ReducedModel.BondLengths), max(pro1.ReducedModel.BondLengths), 1000)))
            title('Probability Distrubition of Bond Lengths')
                hold on
                plot(linspace(min(pro2.ReducedModel.BondLengths), max(pro2.ReducedModel.BondLengths), 1000), ...
                     pdf(fitdist(pro2.ReducedModel.BondLengths, 'kernel'), ...
                     linspace(min(pro2.ReducedModel.BondLengths), max(pro2.ReducedModel.BondLengths), 1000)))
                hold off
    subplot(2,2,2);% bond angles
        plot(linspace(-180, 180, 1000), pdf(fitdist(pro1.ReducedModel.BondAngles, 'kernel'), linspace(-180, 180, 1000)))
            title('Probability Distrubition of Bond Angles');
            xticks([-180:45:+180]);
                hold on
                plot(linspace(-180, 180, 1000), pdf(fitdist(pro2.ReducedModel.BondAngles, 'kernel'), linspace(-180, 180, 1000)))
                hold off
    subplot(2,2,[3,4]); % torsion angles
        plot(linspace(-180, 180, 1000), pdf(fitdist(pro1.ReducedModel.TorsionalAngles, 'kernel'), linspace(-180, 180, 1000)))
            title('Probability Distrubition of Torsion Angles');
            xticks([-180:45:+180]);
                hold on
                plot(linspace(-180, 180, 1000), pdf(fitdist(pro2.ReducedModel.TorsionalAngles, 'kernel'), linspace(-180, 180, 1000)))
                hold off
 