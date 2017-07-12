plates.CellLine = categorical(plates.CellLine);

mode = use; %use_2ch use_all
c = 1;

%for c = 1:4

cellplates = plates(plates.CellLine == pORACL{c},:);
plateIDs = cellplates.expt_plate;
plate = Load_Batch_LE(cellplates,plateIDs, profiles_folder,qc_folder);
anno_updater;

[name,category_id,channel_id,use,use_2ch,use_all]= import_featuresselection2('W:\2015_09_HTS_LE\Code_LE\LE_Analysis\LE_MatlabScripts\feature_select_NBT.txt');
use = use==1;
use_2ch = use_2ch==1;
use_all = use_all==1;

plate.feaInfo.name = plate.feaInfo.name(mode);
plate.feaInfo.channel_id = plate.feaInfo.channel_id(mode);
plate.feaInfo.category_id = plate.feaInfo.category_id(mode);

for w = 1:size(plate.profiles,1)
    plate.profiles{w,1} = plate.profiles{w,1}(mode);
end

plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);

[plate_bioactive,...
    plate3,...
    plate_ref_bioactive,...
    plate_ref,...
    plate_cpd_bioactive,...
    plate_cpd] = Get_Bioactive_Compounds(plate2, HCSparas);

is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
bioactive = Plate2Table(plate_cpd_bioactive);

bioactive_use = bioactive;
%bioactive_use_2ch = bioactive;
%bioactive_use_all = bioactive;
clear bioactive

bioactive_use(:,'profiles_1') = [];
bioactive_use_2ch(:,'profiles_1') = [];
bioactive_use_all(:,'profiles_1') = [];
uniqID = cell(1,1);
for u = 1:height(bioactive_use)
uniqID(u,1) = strcat(bioactive_use.plateIDs(u,1),bioactive_use.well_names(u,1));
end
uniqID = cell2table(uniqID);
bioactive_use = horzcat(uniqID, bioactive_use);
bioactive_use_2ch = horzcat(uniqID, bioactive_use_2ch);
bioactive_use_all = horzcat(uniqID, bioactive_use_all);

bioactive = join(bioactive_use,bioactive_use_2ch,'LeftKeys',1,'RightKeys',1);
bioactive = join(bioactive,bioactive_use_all,'LeftKeys',1,'RightKeys',1);

bioactive(:,{'plateIDs_bioactive_use','well_names_bioactive_use'}) = [];
bioactive(:,{'plateIDs_bioactive_use_2ch','well_names_bioactive_use_2ch','drug_categories_bioactive_use_2ch','drug_names_bioactive_use_2ch','dose_bioactive_use_2ch','concentrations_bioactive_use_2ch'}) = [];
bioactive(:,{'targets_bioactive_use_2ch','pathways_bioactive_use_2ch','descriptions_bioactive_use_2ch','cas_bioactive_use_2ch','cpd_usage_bioactive_use_2ch','nCells_bioactive_use_2ch','plate_types_bioactive_use_2ch','clone_names_bioactive_use_2ch','plateIDs','well_names','drug_categories','drug_names','dose','concentrations'}) = [];
bioactive(:,{'targets','pathways','descriptions','cas','cpd_usage','nCells','plate_types','clone_names'}) = [];
bioactive(:,'cas_bioactive_use') = [];
bioactive(:,'drug_categories_bioactive_use') = [];

bioactive.cpd_usage_bioactive_use=categorical(bioactive.cpd_usage_bioactive_use);
query = bioactive(bioactive.cpd_usage_bioactive_use =='query_cpd',:);
percents.bioactive.NBT = sum(query.bioactive_pVals_bioactive_use <= pVal_thr)/height(compound_IDs)*100;
percents.bioactive.CH1_3 = sum(query.bioactive_pVals_bioactive_use_2ch <= pVal_thr)/height(compound_IDs)*100;
percents.bioactive.ALL = sum(query.bioactive_pVals <= pVal_thr)/height(compound_IDs)*100;

query_table = struct2table(percents.bioactive);

%end