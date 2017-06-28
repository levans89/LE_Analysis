function plate = Load_Batch_LE(plates,plate_IDs, profiles_folder,qc_folder)
% Load and merge multiple saved plates.
% plate = Load_Batch( plateIDs)
% Input:
%     plate_IDs: Cell array of plate numbers, e.g. {'2013002090','2013002091',...,'2013002097'}
%     profiles_folder: path to the folder containing saved profiles.
% Output:
%     plate is a plate struct containing all plate data.
% Author: Chien-Hsiang Hsu at the Altschuler and Wu Laboratories, modified by Louise Heinrich 2017.03.29
plate = load([profiles_folder filesep 'plate_' plate_IDs{1} '.mat']);
plate = plate.plate;

if length(plate_IDs) > 1
    for p = 2:length(plate_IDs)
        try
            p_tmp = load([profiles_folder filesep 'plate_' plate_IDs{p} '.mat']);
        catch
            print('Pbl with plate' + plate_IDs{p})
            continue;
        end
        p_tmp = p_tmp.plate;
        plate = Merge_Plates(plate,p_tmp);
    end
end

% Remove bad wells
plate = Remove_Bad_Wells_from_Plate_LE(plates,plate,qc_folder);
end

