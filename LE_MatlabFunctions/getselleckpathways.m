 function [targets, pathways, compound_IDs]=getselleckpathways()
[SelleckBioactivesSMDC384wellmapping] = import_compoundplates;
%% change names
Selleck = SelleckBioactivesSMDC384wellmapping;%rename
%% extract info
Selleck.PLATE_384 = categorical(Selleck.PLATE_384);%assign categorical variable
targets = cell2table(Selleck.Target);
pathways = cell2table(Selleck.Pathway);
compound_IDs = cell2table(unique(Selleck.Compound_ID));
 end