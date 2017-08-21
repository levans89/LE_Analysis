% anno_updater()

% plate = Load_Batch_LE(plates,plate_IDs, profiles_folder_NBT,qc_folder); % note different profiles folder

for r = 1:length(plate.drug_categories)
plate.drug_categories{r,1} = strrep(plate.drug_categories{r,1}, 'AuroraB', 'Cell Cycle');
plate.drug_categories{r,1} = strrep(plate.drug_categories{r,1}, 'PLK1','Cell Cycle');
plate.drug_categories{r,1} = strrep(plate.drug_categories{r,1}, 'DNA Damage', 'DNA');
plate.drug_categories{r,1} = strrep(plate.drug_categories{r,1}, 'PI3K/Akt/mTOR','mTOR');
plate.drug_categories{r,1} = strrep(plate.drug_categories{r,1}, 'Proteases','Proteasome');
end

for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1}, 'JAK');
if a>=1
plate.drug_categories{r,1} = 'JAKSTAT' ;
end
end
clear a

for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1}, 'Integrin');
if a>=1
plate.drug_categories{r,1} = 'Cytoskeleton' ;
end
end
clear a

for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1}, 'HSP (e.g. HSP90)');
b= strfind(plate.targets{r,1}, 'HSP90');
if a>=1
plate.drug_categories{r,1} = 'HSP90' ;
elseif b>=1
plate.drug_categories{r,1} = 'HSP90' ;
end
end
clear a b

for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1}, 'HDAC');
if a>=1
plate.drug_categories{r,1} = 'HDAC' ;
else
end
end
clear a

for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1}, 'Kinesin');
if a>=1
plate.drug_categories{r,1} = 'Cell Cycle' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
PARP = strfind(plate.targets{r,1}, 'PARP');
if PARP>=1
plate.drug_categories{r,1} = 'PARP' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
BET = strfind(plate.targets{r,1}, 'BET');
if BET>=1
plate.drug_categories{r,1} = 'BET' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.descriptions{r,1},'glucocorticoid');
if a >=1
plate.drug_categories{r,1} = 'Glucocorticoid' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1}, 'Estrogen/progestogen Receptor');
if a>=1
plate.drug_categories{r,1} = 'EstrogenR' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1}, 'Histone Methyltransferase');
if a>=1
plate.drug_categories{r,1} = 'HistoneMT' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'Sirtuin');
if a>=1
plate.drug_categories{r,1} = 'Sirtuin' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'Epigenetic Reader Domain');
if a>=1
plate.drug_categories{r,1} = 'Epigenetic Reader Domain' ;
else
end
end
clear a


for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'Histone Acetyltransferase');
if a>=1
plate.drug_categories{r,1} = 'Histone Acetyltransferase' ;
else
end
end
clear a


for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'DNA Methyltransferase');
if a>=1
plate.drug_categories{r,1} = 'DNA Methyltransferase' ;
else
end
end
clear a

for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'HDAC');
if a >=1
plate.drug_categories{r,1} = 'HDAC' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'Histone demethylases');
if a >=1
plate.drug_categories{r,1} = 'Histone demethylases' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'COX');
if a >=1
plate.drug_categories{r,1} = 'COX' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'HDAC3');
if a >=1
plate.drug_categories{r,1} = 'HDAC' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'MEK');
if a >=1
plate.drug_categories{r,1} = 'MEK' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'HMG-CoA Reductase');
if a >=1
plate.drug_categories{r,1} = 'HMG-CoA Reductase' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'PPAR');
if a >=1
plate.drug_categories{r,1} = 'PPAR' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'Histamine Receptor');
if a >=1
plate.drug_categories{r,1} = 'Histamine Receptor' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'5-HT Receptor');
if a >=1
plate.drug_categories{r,1} = '5-HT Receptor' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'AChR');
if a>=1
plate.drug_categories{r,1} = 'AChR';
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'GSK-3');
if a >=1
plate.drug_categories{r,1} = 'GSK-3' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'PLK1');
if a >=1
plate.drug_categories{r,1} = 'Cell Cycle' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'DUB');
if a >=1
plate.drug_categories{r,1} = 'Ubiquitin' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'EGFR');
if a >=1
plate.drug_categories{r,1} = 'EGFR' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'HER2');
if a ==1
plate.drug_categories{r,1} = 'EGFR' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'ALK');
if a >=1
plate.drug_categories{r,1} = 'ALK' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'Microtubule Associated');
if a ==1
plate.drug_categories{r,1} = 'MT' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'VEGFR');
if a >=1
plate.drug_categories{r,1} = 'Angiogenesis' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'FGFR');
if a>=1
plate.drug_categories{r,1} = 'Angiogenesis' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'PDGFR');
if a>=1
plate.drug_categories{r,1} = 'Angiogenesis' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.targets{r,1},'p38');
if a >=1
plate.drug_categories{r,1} = 'MAPK' ;
else
end
end
clear a
for r = 1:length(plate.drug_categories)
a = strfind(plate.descriptions{r,1},'retino');
if a >=1
plate.drug_categories{r,1} = 'Retinoid' ;
else
end
end
clear a

for r = 1:length(plate.drug_categories)
a = strfind(plate.drug_names{r,1}, 'Geldanamycin');
if a>=1
plate.drug_categories{r,1} = 'HSP90' ;
end
end
clear a b



drug_categories = categories(categorical(plate.drug_categories))';
