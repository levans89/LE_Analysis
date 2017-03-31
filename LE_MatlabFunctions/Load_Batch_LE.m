function plate = Load_Batch_LE( plateIDs, profile_folder)
% Load and merge multiple saved plates.
% 
% plate = Load_Batch( plateIDs)
% 
% Input:
%     plateIDs: Cell array of plate numbers, e.g. {'2013002090','2013002091',...,'2013002097'}
%     profile_folder: path to the folder containing saved profiles.
% 
% Output:
%     plate is a plate struct containing all plate data.
% 
% Note:
%     The plate must be saved to have the "plate_" as prefix,for example, "plate_2013002097.mat"
% 
% Chien-Hsiang Hsu 2015.12.16
% 
% Copyright: Altschuler and Wu laboratories
% Author: Chien-Hsiang Hsu at the Altschuler and Wu Laboratories, modified
% by Louise Heinrich 2017.03.29
% For latest updates, check: <http://www.altschulerwulab.org/>.
%
% All rights reserved.
% Louise altered to remove duplication of input plateIDs{2}, and no
% plateIDs {1}.



plate = load([profile_folder filesep 'plate_' plateIDs{1} '.mat']);
plate = plate.plate;

if length(plateIDs) > 1
    for p = 2:length(plateIDs)
        try
            p_tmp = load([profile_folder filesep 'plate_' plateIDs{p} '.mat']);
        catch
            print('Pbl with plate '+plateIDs{p})
            continue;
        end
        p_tmp = p_tmp.plate;
        plate = Merge_Plates(plate,p_tmp);
    end
end

% Remove bad wells
plate = Remove_Bad_Wells_from_Plate(plate);
end

