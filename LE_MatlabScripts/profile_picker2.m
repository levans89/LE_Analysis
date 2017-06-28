clear profiles
close all
s = 1
f = height(plate)

for p=s:f %1:size(plate2.drug_names,1)
    profiles(p,:) = plate.profiles{p};
end
figure
clims = [-1 1];
imagesc(profiles(s:f,:),clims)
colormap(jet)
