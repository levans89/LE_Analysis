clear profiles2

s = 1
f = length(profiles)

for p=s:f 
    profiles2(p,:) = profiles{p,1};
end

figure
clims = [-1 1];
imagesc(profiles2(s:f,:),clims)
colormap(jet)
