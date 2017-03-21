function bad_wells = LE_Bad_Wells(plates, qc_folder)
% Record all bad/weird wells
% plates is a table made from currexp in main_4lines.m which contains the
% plate names of all plates in this case (cell line) being analyzed
% bad_wells looks for an excel file containing sheets describing bad wells
% and fluo drugs in 0/1 format. each plate should have one file. plate is
% laid out in physical plate format.
% qc_folder = 'W:/2015_09_HTS_LE/QC/4lines/'
% Louise Evans 2017.03.13 (modified from Chien-Hsiang Hsu Bad_Wells.m)

platemap = txt_importfile('W:/2015_09_HTS_LE/Code/platemap_384.xlsx', 'WELLS', 'B2:Y17'); % standard layout of plate with wellname

col_names = {'PlateID','Well','Comment'};
bad_wells = cell(1,3);
bad_wells = cell2table(bad_wells,'VariableNames',col_names); % make table ready for writing to

for x = 1: size(plates.plate_name,1)
    try
        bad_QC = importfile(char(fullfile(qc_folder,strcat(plates.expt_plate(x),'.xlsx'))),'BADWELLS','B2:Y17'); % import QC files marking BADWELLS and FLUODRUGS from qc_folder
        bad_QC = logical(cell2mat(bad_QC));% find wells containing 1 (positive)
        bad_QC = platemap(bad_QC~=0);% extract wellnames for those wells
        for y = 1:size(bad_QC,1)
            bad_QC{y,2}='bad_QC';% insert comment for each well describing issue
        end;
        
        fluodrugs = importfile(char(fullfile(qc_folder,strcat(plates.expt_plate(x),'.xlsx'))),'FLUODRUGS','B2:Y17');
        fluodrugs=logical(cell2mat(fluodrugs));
        fluodrugs = platemap(fluodrugs~=0);
        for z = 1:size(fluodrugs,1)
            fluodrugs{z,2}='fluo_drug';
        end;
        
        clear y z
        
        prebw = cell2table(vertcat(bad_QC,fluodrugs));% precursor table
        prebw.Properties.VariableNames{1} = 'Well';% make prebw have the same headers as bad_wells so can join together
        prebw.Properties.VariableNames{2} = 'Comment';
        PlateID=repmat(plates.expt_plate(x),height(prebw),1);% assign plateID to each well
        PlateID = cell2table(PlateID);%reformat
        prebw = horzcat(PlateID,prebw);%join
        bad_wells = [bad_wells;prebw];% join new bad wells for current plate onto existing list
    catch
        warning('QC file for plate not found')   % stops function crashing if missing
    end
end
bad_wells(1,:) = []; % remove placeholder first row

