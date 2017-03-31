 function [targets, pathways, compound_IDs]=getselleckpathways()
[SelleckBioactivesSMDC384wellmapping] = import_compoundplates;
%% change names
Selleck = SelleckBioactivesSMDC384wellmapping;%rename
%% extract info
Selleck.PLATE_384 = categorical(Selleck.PLATE_384);%assign categorical variable
targets = cell2table(Selleck.Target, 'VariableNames',{'Target'});
pathways = cell2table(Selleck.Pathway, 'VariableNames',{'Pathway'});
compound_IDs = cell2table(Selleck.Compound_ID, 'VariableNames',{'drug_names'});
compound_IDs = horzcat(compound_IDs,targets,pathways);
 end