clear profiles2 profiles

s = 1
f = 42

for p=s:f 
    profiles{p,:} = bioactive4lines.bioactivity_tables5{1,1}.profiles{p}
    profiles2(p,:) = profiles{p,1};
end

figure
clims = [-1 1];
imagesc(profiles2(s:f,:),clims)
colormap(jet)
