plateIDs2 = plate.plateIDs;
plateIDs2 = categorical(plateIDs2);
categories(plateIDs2)
wellcounts = countcats(plateIDs2);
if max(wellcounts)> 384
    warning ('problem with one of your plates')
else
    display (['min',{min(wellcounts)}])
    display (['max',{max(wellcounts)}])
end