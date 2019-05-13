function [coords] = snapcoords(snapname)
% snapname = 'C:\Users\hncetin\OneDrive\Belgeler\MATLAB\MATLAB Drive\process-pdb\snapshots\18.pdb';
    formatSpec = '%*26s%12f%8f%8f%[^\n\r]';
    fileID = fopen(snapname,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    coords = [dataArray{1:end-1}];
    clearvars filename formatSpec fileID dataArray ans;
end