 function [targets, pathways, compound_IDs]=getselleckpathways()
 % Louise Heinrich 2017.04.04 in Altschuler and Wu Laboratories
 
 % [targets, pathways, compound_IDs]=getselleckpathways()
 % uses import_compoundplates to import the Selleck 2K library raw data
 % file
 % outputs targets, pathways, compound_IDs as tables 
 
[SelleckBioactivesSMDC384wellmapping] = import_compoundplates;
%% change names
Selleck = SelleckBioactivesSMDC384wellmapping;%rename
%% extract info
Selleck.PLATE_384 = categorical(Selleck.PLATE_384);%assign categorical variable
targets = cell2table(Selleck.Target, 'VariableNames',{'Target'});
pathways = cell2table(Selleck.Pathway, 'VariableNames',{'Pathway'});
compound_IDs = cell2table(Selleck.Compound_ID, 'VariableNames',{'drug_names'});
compound_IDs = horzcat(compound_IDs,targets,pathways);
for s=1:height(compound_IDs)
    for c=1:3
        compound_IDs{s,c}=strtrim(compound_IDs{s,c});
    end
end
 end