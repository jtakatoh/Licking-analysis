% This script copies files to a new folder after extracting them from the original folder.
% The initial step in creating training data is to find suitable frames and
% make masks. Because use of the same side and bottom frames is preferable for
% training data. The names of the side and bottom frames must match.
% If you made masks from the side frames, 
% you can move the png files from the original bottom folder to a new folder 
% by using names/frames used in side frames. 
%%
% change directory to the mask folder
fclose all;
fileinfo = dir('*.png');
fnames = {fileinfo.name};
%%
% define the source folder that contains all the png files
psource = '/Users/jun/Documents/Work/Project/Licking/Training Data/0';
% define the destination
pdest = '/Users/jun/Documents/Work/Project/Licking/Training Data/origina_bottom_20';

% This script moves files from the original to destination. If you need to
% keep files in the original folder, they can be copied back.

for i = 1:length(fnames)
    sourceFile = string(fullfile(psource, fnames(i)));
    destFile   = string(fullfile(pdest, fnames(i)));
    fclose all;
    movefile(sourceFile, destFile);
end