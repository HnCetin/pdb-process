clear

%% Import coordinate data from snapshots
for i=0:200
    filename = [num2str(i), '.pdb'];
    snap_coords{i+1} = snapcoords(filename);
end

%% Calculate Autocorrelations

% Mean values
sum = snap_coords{1,1};
for i=2:201
    sum = sum + snap_coords{1,i} ;
end
means = sum / 201;

% Variances
for i=1:201
    autocorr{i} = (snap_coords{1,i} - means).^2;
end
clear i sum filename

 plot(autocorr)
 
 
%% Calculate Covariances

for j=1:201
    for i=1:201
        covars{j,i} = (snap_coords{1,j} - means) .* (snap_coords{1,i} - means);
    end
end


save workspace.mat  snap_coords means covars autocorr

clear
load workspace.mat