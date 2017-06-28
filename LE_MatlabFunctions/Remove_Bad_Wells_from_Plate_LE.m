function plate = Remove_Bad_Wells_from_Plate_LE(plates,plate,qc_folder)
% Remove bad wells from the given plate structure.
% The list of bad wells is obtained from Bad_Wells.m
% plate = Remove_Bad_Wells_from_Plate( plate )
% Chien-Hsiang Hsu, 2016.04.20

%% 1. Get plateID_bad_well labels
bws = LE_Bad_Wells(plates,qc_folder);
all_bad_wells = strcat(bws.PlateID,{'_'},bws.Well);
plate_wells = strcat(plate.plateIDs,{'_'},plate.well_names);
%% 2. Remove bad wells
is_bad = ismember(plate_wells,all_bad_wells);
plate = Select_Rows(plate,~is_bad);

end

