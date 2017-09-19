function [] = folderunpacker(folderpath)
%FOLDERUNPACKER Summary of this function goes here
%   Detailed explanation goes here
filelist = dir(folderpath);
for f = 3:length(filelist)
    source = [folderpath,'\',filelist(f).name]
    movefile([source,'\*'],folderpath)
    rmdir(source)
end
end

