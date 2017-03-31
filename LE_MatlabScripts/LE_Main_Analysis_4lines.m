% Louise Evans March 2017
% Purpose: Perform analogous analysis to Charles' NBT paper for 4lines
% experiment x Selleck 2K set 'pORACL' determination.
%% Set paths
%====Environment setup====%
addpath(genpath('W:\Lab_Analysis\common\plate_annotation')) %plate annotation from labtoolbox
addpath(genpath('W:\2015_09_HTS_LE\Code\CD_Tag_Drug_Screen')) % Charles code for NBT paper containing useful functions
addpath(genpath('W:\2015_09_HTS_LE\Code\LE_analysis')) % Louise code afer cluster operations
addpath(genpath('W:\2015_09_HTS_LE\Code\LE_Preprocessing')) % Louise code before cluster operations
profiles_folder = 'W:\2015_09_HTS_LE\data\profiles';% place to save profiles to and load from
qc_folder = 'W:\2015_09_HTS_LE\QC\4lines\';% store any QC analysis Python, Matlab, or manual
results_folder = 'W:\2015_09_HTS_LE\results\matlab_results\4lines'; % Louise results for 4lines analysis
load('W:\2015_09_HTS_LE\data\currexp.mat');% load in currexp.mat which details all plates in this screen (extracted from exp_DB)

%Mac Paths
addpath(genpath('../../Lab_Analysis/common/plate_annotation')) %plate annotation from labtoolbox
addpath(genpath('../CD_Tag_Drug_Screen')) %plate annotation from labtoolbox
addpath(genpath('../LE_Analysis')) %plate annotation from labtoolbox
addpath(genpath('../LE_Preprocessing')) %plate annotation from labtoolbox
profiles_folder = '../data/profiles';% place to save profiles to and load from
qc_folder = '../QC/4lines';% store any QC analysis Python, Matlab, or manual
results_folder = '../results/matlab_results/4lines'; % Louise results for 4lines analysis
load('../data/currexp.mat');% load in currexp.mat which details all plates in this screen (extracted from exp_DB)

currexp.Properties.VariableNames{8} = 'ProductNoPlate'; % change unclear name to avoid confusion with plate_type
currexp.CellLine = categorical(currexp.CellLine); %categorical to select by cell line later
currexp.batch = categorical(currexp.batch); % categorical to select by batch (1-8) later
currexp.cpdAdditionProtocol = categorical(currexp.cpdAdditionProtocol); % later used to remove unwanted plate
currexp.plate_type = categorical(currexp.plate_type); % use to switch between reference_plate and screen_plate
currexp = currexp(currexp.CellLine ~='MEDIA' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7' & currexp.cpdAdditionProtocol ~= '50nL Echo Transfer',:); % remove bad batches, illum correction plate and extra plate
clc
%% Profile a batch of plates Charles
%====Set parameters and enable switching between cell lines "pORACLs"====%
pORACL = 'XRCC5'; % we screened 4 lines - here we switch profile generation between the lines to deal with marker encoding in feature names
switch pORACL %potential ORACL
    case 'NQO1' % S \ _Alice method
        features_to_use = feature_selector('NBT_NQO1'); %{'all'};%Get_Default_Features(); % different for cell lines as plate annotation encodes reporter name in feature
        plates = currexp(currexp.CellLine=='A549_NQO1',:); %select cell line, exclude dud batches
        plateIDs = plates.expt_plate; % extract plate IDs
        bad_wells = LE_Bad_Wells(plates,qc_folder); % Get all bad wells from per-plate manual Excel annotation
        %writetable(bad_wells,'W:\2015_09_HTS_LE\QC\4lines\bw_NQO1.xlsx'); % write to a single file can therefore choose not to use later
        plate_maps = cell(height(plates),1); % preallocate
        for m = 1:height(plates) % for all plates in this cell line
            plate_maps{m} = [num2str(plates.expt_plate{m}),'.xlsx']; % find the platemap name
        end
        plate_types = plates.plate_type; % assign the plate type e.g. reference_plate, screen_plate
        paras = Set_Paras(); % set parameters (path,bioactive,screen,cluster)
        plate_maps_folder = paras.path.plate_maps; % where to find plate maps
        h = waitbar(0,num2str(0)); % display waitbar
        for p = 1:length(plateIDs) % for all plates
            waitbar(p\length(plateIDs),h,[num2str(p),'\',num2str(length(plateIDs))]) % change the waitbar
            
            % 1. Get all inputs for Get_KS_Profiles
            plateID = plateIDs{p};% ID
            plate_type = plate_types{p}; % e.g. reference_plate
            plate_map_file = [plate_maps_folder filesep plate_maps{p}]; % expt plate file
            plate_layout = Read_Plate_Annotation(plate_map_file); % read in expt plate file
            
            % 2. Gather information for profiling
            well_names = {plate_layout.wellData.wellName}';
            drug_names = {plate_layout.wellData.Compound_ID}';
            drug_categories = {plate_layout.wellData.Compound_Category}';
            clone_names = {plate_layout.wellData.Cell_Line}';
            cpd_usage = {plate_layout.wellData.Compound_Usage}';
            dose = {plate_layout.wellData.Dose}';
            cas = {plate_layout.wellData.CAS_Number}';
            targets = {plate_layout.wellData.Target}';
            pathways = {plate_layout.wellData.Pathway}';
            descriptions = {plate_layout.wellData.Description}';
            concentrations = str2double({plate_layout.wellData.Dose_Category})';
            concentrations(isnan(concentrations))=0;
            
            % 3. Define your KS control %%%
            is_ks_ctrl = ismember(well_names(:),{'A2','A23',...
                'B2','B23',...
                'C2','C23',...
                'D2','D23',...
                'E2','E23',...
                'F2','F23',...
                'G2','G23',...
                'H2','H23',...
                'I2','I23',...
                'J2','J23',...
                'K2','K23',...
                'L2','L23',...
                'M2','M23',...
                'N2','N23',...
                'O2','O23',...
                'P2','P23'});
            
            % 4. Eliminate bad wells
            plate_bad_wells = bad_wells.Well(strcmp(bad_wells.PlateID,plateID));
            is_bad = ismember(well_names,plate_bad_wells);
            
            well_names      = well_names(~is_bad);
            drug_names      = drug_names(~is_bad);
            drug_categories = drug_categories(~is_bad);
            clone_names     = clone_names(~is_bad);
            cpd_usage       = cpd_usage(~is_bad);
            dose            = dose(~is_bad);
            cas             = cas(~is_bad);
            targets         = targets(~is_bad);
            pathways        = pathways(~is_bad);
            descriptions    = descriptions(~is_bad);
            concentrations  = concentrations(~is_bad);
            is_ks_ctrl      = is_ks_ctrl(~is_bad);
            
            % 5. Calculate profiles
            Get_KS_Profiles(plateID, plate_type, well_names, drug_names, drug_categories, concentrations,...
                dose, clone_names, cas, targets, pathways, descriptions, ...
                is_ks_ctrl, cpd_usage, paras,...
                'save_destination',profiles_folder,'features_to_use',features_to_use);
            
        end
        waitbar(p\length(plateIDs),h,'Done')
        delete(h)
        
    case 'XRCC5'% S \ _Alice method
        features_to_use = feature_selector('NBT_XRCC5');
        plates = currexp(currexp.CellLine=='A549_XRCC5' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plateIDs = plates.expt_plate;
        bad_wells = LE_Bad_Wells(plates,qc_folder);
        %writetable(bad_wells,'W:\2015_09_HTS_LE\QC\4lines\bw_XRCC5.xlsx');
        plate_maps = cell(height(plates),1);
        for m = 1:height(plates)
            plate_maps{m} = [num2str(plates.expt_plate{m}),'.xlsx'];
        end
        plate_types = plates.plate_type;
        paras = Set_Paras();
        paras.path.features = ['features' filesep 'cbfeatures_Alice' filesep 'cbfeatures-']; % switch feature folder here to get S method
        plate_maps_folder = paras.path.plate_maps;
        h = waitbar(0,num2str(0));
        for p = 2:length(plateIDs)
            waitbar(p\length(plateIDs),h,[num2str(p),'\',num2str(length(plateIDs))])
            
            % 1. Get all inputs for Get_KS_Profiles
            plateID = plateIDs{p};
            plate_type = plate_types{p};
            plate_map_file = [plate_maps_folder filesep plate_maps{p}];
            plate_layout = Read_Plate_Annotation(plate_map_file);
            
            % 2. Gather information for profiling
            well_names = {plate_layout.wellData.wellName}';
            drug_names = {plate_layout.wellData.Compound_ID}';
            drug_categories = {plate_layout.wellData.Compound_Category}';
            clone_names = {plate_layout.wellData.Cell_Line}';
            cpd_usage = {plate_layout.wellData.Compound_Usage}';
            dose = {plate_layout.wellData.Dose}';
            cas = {plate_layout.wellData.CAS_Number}';
            targets = {plate_layout.wellData.Target}';
            pathways = {plate_layout.wellData.Pathway}';
            descriptions = {plate_layout.wellData.Description}';
            concentrations = str2double({plate_layout.wellData.Dose_Category})';
            concentrations(isnan(concentrations))=0;
            
            % 3. Define your KS control %%%
            is_ks_ctrl = ismember(well_names(:),{'A2','A23',...
                'B2','B23',...
                'C2','C23',...
                'D2','D23',...
                'E2','E23',...
                'F2','F23',...
                'G2','G23',...
                'H2','H23',...
                'I2','I23',...
                'J2','J23',...
                'K2','K23',...
                'L2','L23',...
                'M2','M23',...
                'N2','N23',...
                'O2','O23',...
                'P2','P23'});
            
            % 4. Eliminate bad wells
            plate_bad_wells = bad_wells.Well(strcmp(bad_wells.PlateID,plateID));
            is_bad = ismember(well_names,plate_bad_wells);
            
            well_names      = well_names(~is_bad);
            drug_names      = drug_names(~is_bad);
            drug_categories = drug_categories(~is_bad);
            clone_names     = clone_names(~is_bad);
            cpd_usage       = cpd_usage(~is_bad);
            dose            = dose(~is_bad);
            cas             = cas(~is_bad);
            targets         = targets(~is_bad);
            pathways        = pathways(~is_bad);
            descriptions    = descriptions(~is_bad);
            concentrations  = concentrations(~is_bad);
            is_ks_ctrl      = is_ks_ctrl(~is_bad);
            
            % 5. Calculate profiles
            
            Get_KS_Profiles(plateID, plate_type, well_names, drug_names, drug_categories, concentrations,...
                dose, clone_names, cas, targets, pathways, descriptions, ...
                is_ks_ctrl, cpd_usage, paras,...
                'save_destination',profiles_folder,'features_to_use',features_to_use);
            
        end
        waitbar(p\length(plateIDs),h,'Done')
        delete(h)
        
    case 'SET'
        features_to_use = feature_selector('NBT_SET');
        plates = currexp(currexp.CellLine=='A549_SET' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plateIDs = plates.expt_plate;
        bad_wells = LE_Bad_Wells(plates,qc_folder);
        %writetable(bad_wells,'W:\2015_09_HTS_LE\QC\4lines\bw_SET.xlsx');
        plate_maps = cell(height(plates),1);
        for m = 1:height(plates)
            plate_maps{m} = [num2str(plates.expt_plate{m}),'.xlsx'];
        end
        plate_types = plates.plate_type;
        paras = Set_Paras();
        plate_maps_folder = paras.path.plate_maps;
        h = waitbar(0,num2str(0));
        for p = 1:length(plateIDs)
            waitbar(p\length(plateIDs),h,[num2str(p),'\',num2str(length(plateIDs))])
            
            % 1. Get all inputs for Get_KS_Profiles
            plateID = plateIDs{p};
            plate_type = plate_types{p};
            plate_map_file = [plate_maps_folder filesep plate_maps{p}];
            plate_layout = Read_Plate_Annotation(plate_map_file);
            
            % 2. Gather information for profiling
            well_names = {plate_layout.wellData.wellName}';
            drug_names = {plate_layout.wellData.Compound_ID}';
            drug_categories = {plate_layout.wellData.Compound_Category}';
            clone_names = {plate_layout.wellData.Cell_Line}';
            cpd_usage = {plate_layout.wellData.Compound_Usage}';
            dose = {plate_layout.wellData.Dose}';
            cas = {plate_layout.wellData.CAS_Number}';
            targets = {plate_layout.wellData.Target}';
            pathways = {plate_layout.wellData.Pathway}';
            descriptions = {plate_layout.wellData.Description}';
            concentrations = str2double({plate_layout.wellData.Dose_Category})';
            concentrations(isnan(concentrations))=0;
            
            % 3. Define your KS control %%%
            is_ks_ctrl = ismember(well_names(:),{'A2','A23',...
                'B2','B23',...
                'C2','C23',...
                'D2','D23',...
                'E2','E23',...
                'F2','F23',...
                'G2','G23',...
                'H2','H23',...
                'I2','I23',...
                'J2','J23',...
                'K2','K23',...
                'L2','L23',...
                'M2','M23',...
                'N2','N23',...
                'O2','O23',...
                'P2','P23'});
            
            % 4. Eliminate bad wells
            plate_bad_wells = bad_wells.Well(strcmp(bad_wells.PlateID,plateID));
            is_bad = ismember(well_names,plate_bad_wells);
            
            well_names      = well_names(~is_bad);
            drug_names      = drug_names(~is_bad);
            drug_categories = drug_categories(~is_bad);
            clone_names     = clone_names(~is_bad);
            cpd_usage       = cpd_usage(~is_bad);
            dose            = dose(~is_bad);
            cas             = cas(~is_bad);
            targets         = targets(~is_bad);
            pathways        = pathways(~is_bad);
            descriptions    = descriptions(~is_bad);
            concentrations  = concentrations(~is_bad);
            is_ks_ctrl      = is_ks_ctrl(~is_bad);
            
            % 5. Calculate profiles
            Get_KS_Profiles(plateID, plate_type, well_names, drug_names, drug_categories, concentrations,...
                dose, clone_names, cas, targets, pathways, descriptions, ...
                is_ks_ctrl, cpd_usage, paras,...
                'save_destination',profiles_folder,'features_to_use',features_to_use);
            
        end
        waitbar(p\length(plateIDs),h,'Done')
        delete(h)
        
    case 'S100A11' % S \ _Alice method
        features_to_use = feature_selector('NBT_S100A11');
        plates = currexp(currexp.CellLine=='A549_S100A11' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plateIDs = plates.expt_plate;
        bad_wells = LE_Bad_Wells(plates,qc_folder);
        %writetable(bad_wells,'W:\2015_09_HTS_LE\QC\4lines\bw_S100A11.xlsx');
        plate_maps = cell(height(plates),1);
        for m = 1:height(plates)
            plate_maps{m} = [num2str(plates.expt_plate{m}),'.xlsx'];
        end
        plate_types = plates.plate_type;
        paras = Set_Paras();
        plate_maps_folder = paras.path.plate_maps;
        h = waitbar(0,num2str(0));
        for p = 1:length(plateIDs)
            waitbar(p\length(plateIDs),h,[num2str(p),'\',num2str(length(plateIDs))])
            
            % 1. Get all inputs for Get_KS_Profiles
            plateID = plateIDs{p};
            plate_type = plate_types{p};
            plate_map_file = [plate_maps_folder filesep plate_maps{p}];
            plate_layout = Read_Plate_Annotation(plate_map_file);
            
            % 2. Gather information for profiling
            well_names = {plate_layout.wellData.wellName}';
            drug_names = {plate_layout.wellData.Compound_ID}';
            drug_categories = {plate_layout.wellData.Compound_Category}';
            clone_names = {plate_layout.wellData.Cell_Line}';
            cpd_usage = {plate_layout.wellData.Compound_Usage}';
            dose = {plate_layout.wellData.Dose}';
            cas = {plate_layout.wellData.CAS_Number}';
            targets = {plate_layout.wellData.Target}';
            pathways = {plate_layout.wellData.Pathway}';
            descriptions = {plate_layout.wellData.Description}';
            concentrations = str2double({plate_layout.wellData.Dose_Category})';
            concentrations(isnan(concentrations))=0;
            
            % 3. Define your KS control %%%
            is_ks_ctrl = ismember(well_names(:),{'A2','A23',...
                'B2','B23',...
                'C2','C23',...
                'D2','D23',...
                'E2','E23',...
                'F2','F23',...
                'G2','G23',...
                'H2','H23',...
                'I2','I23',...
                'J2','J23',...
                'K2','K23',...
                'L2','L23',...
                'M2','M23',...
                'N2','N23',...
                'O2','O23',...
                'P2','P23'});
            
            % 4. Eliminate bad wells
            plate_bad_wells = bad_wells.Well(strcmp(bad_wells.PlateID,plateID));
            is_bad = ismember(well_names,plate_bad_wells);
            
            well_names      = well_names(~is_bad);
            drug_names      = drug_names(~is_bad);
            drug_categories = drug_categories(~is_bad);
            clone_names     = clone_names(~is_bad);
            cpd_usage       = cpd_usage(~is_bad);
            dose            = dose(~is_bad);
            cas             = cas(~is_bad);
            targets         = targets(~is_bad);
            pathways        = pathways(~is_bad);
            descriptions    = descriptions(~is_bad);
            concentrations  = concentrations(~is_bad);
            is_ks_ctrl      = is_ks_ctrl(~is_bad);
            
            % 5. Calculate profiles
            Get_KS_Profiles(plateID, plate_type, well_names, drug_names, drug_categories, concentrations,...
                dose, clone_names, cas, targets, pathways, descriptions, ...
                is_ks_ctrl, cpd_usage, paras,...
                'save_destination',profiles_folder,'features_to_use',features_to_use);
            
        end
        waitbar(p\length(plateIDs),h,'Done')
        delete(h)
end

%% Check cell number
%====Generate heatmaps of cell number for all plates to check for pattern====%
cd ('W:\2015_09_HTS_LE\QC\4lines\') % move to the location to do all the saving in to improve speed and make checking for files quicker
for b = 55:132 %plate IDs last 3 digits
    pid = [b]; % first check blanks, ref plates, compound plates
    plateIDs = arrayfun(@(x)sprintf('2017018%03d',x),pid,'unif',false)'; % did not trouble to change to my plate_IDs method for one-off task
    filename =char(strcat(plateIDs(1),'.bmp'));
    display(filename)
    if exist(filename,'file')==0
        try % in case file is missing
            plate = Load_Batch_LE(plateIDs,profiles_folder); % load plate
            Show_Cell_Number(plate); % generate heatmap of cell number
            saveas(gcf,plateIDs{1,1},'bmp'); % save
        catch
            warning (char(strcat(filename, ' cannot be generated'))) % 
        end
    else
        warning (char(strcat(filename, ' already exists')))
    end
end
close all % close all the figures
cd ('W:\2015_09_HTS_LE\Code\') % move back
clc
%% Plot P Values on plate to view where the cutoff could or should be defined
time_to_use     = 1; % only T48 here
merge_method    = 'concatenate'; % 'concatenate' | 'pool' % doesnt matter only 1 timepoint
pVal_thr        = 1; % p-value threshold bioactivity
HCSparas = Set_Paras('pVal_thr',pVal_thr); % set parameteres
        
pORACL = 'S100A11';
switch pORACL
    case 'NQO1'
        plates = currexp(currexp.CellLine=='A549_NQO1'& currexp.plate_type == 'reference_plate',:); % select cell line & currexp.plate_type == 'reference_plate'
        %plates = currexp(currexp.plate_type == 'reference_plate',:); % select reference plates only
        plate_IDs = plates.expt_plate;
        for p = 1:length(plate_IDs)
            plate = Load_Batch_LE(plate_IDs(p),profiles_folder);
            plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
            [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
            plate2.logbioactive_pVals = log(plate2.bioactive_pVals);
            plate2.logbioactive_pVals(isnan(plate2.logbioactive_pVals))=1;%A(isnan(A)) = 0
            Show_PVALUE(plate2)
            title = ['../results/matlab_results/4lines/',plate_IDs{p},' Log P Value Plate'];
            print('-f1',title,'-dpng','-r500')
            close all
        end
        
    case 'SET'
        plates = currexp(currexp.CellLine=='A549_SET'& currexp.plate_type == 'reference_plate',:); % select cell line & currexp.plate_type == 'reference_plate'
        %plates = currexp(currexp.plate_type == 'reference_plate',:); % select reference plates only
        plate_IDs = plates.expt_plate;
        for p = 1:length(plate_IDs)
            plate = Load_Batch_LE(plate_IDs(p),profiles_folder);
            plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
            [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
            plate2.logbioactive_pVals = log(plate2.bioactive_pVals);
            plate2.logbioactive_pVals(isnan(plate2.logbioactive_pVals))=1;%A(isnan(A)) = 0
            Show_PVALUE(plate2)
            title = ['../results/matlab_results/4lines/',plate_IDs{p},' Log P Value Plate'];
            print('-f1',title,'-dpng','-r500')
            close all
        end        
        
    case 'XRCC5'
        plates = currexp(currexp.CellLine=='A549_XRCC5'& currexp.plate_type == 'reference_plate',:); % select cell line & currexp.plate_type == 'reference_plate'
        %plates = currexp(currexp.plate_type == 'reference_plate',:); % select reference plates only
        plate_IDs = plates.expt_plate;
        for p = 1:length(plate_IDs)
            plate = Load_Batch_LE(plate_IDs(p),profiles_folder);
            plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
            [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
            plate2.logbioactive_pVals = log(plate2.bioactive_pVals);
            plate2.logbioactive_pVals(isnan(plate2.logbioactive_pVals))=1;%A(isnan(A)) = 0
            Show_PVALUE(plate2)
            title = ['../results/matlab_results/4lines/',plate_IDs{p},' Log P Value Plate'];
            print('-f1',title,'-dpng','-r500')
            close all
        end
        
    case 'S100A11'
        plates = currexp(currexp.CellLine=='A549_S100A11'& currexp.plate_type == 'reference_plate',:); % select cell line & currexp.plate_type == 'reference_plate'
        %plates = currexp(currexp.plate_type == 'reference_plate',:); % select reference plates only
        plate_IDs = plates.expt_plate;
        for p = 1:length(plate_IDs)
            plate = Load_Batch_LE(plate_IDs(p),profiles_folder);
            plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
            [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
            plate2.logbioactive_pVals = log(plate2.bioactive_pVals);
            plate2.logbioactive_pVals(isnan(plate2.logbioactive_pVals))=1;%A(isnan(A)) = 0
            Show_PVALUE(plate2)
            title = ['../results/matlab_results/4lines/',plate_IDs{p},' Log P Value Plate'];
            print('-f1',title,'-dpng','-r500')
            close all
        end        
end
%% Visualize plates and specify bioactives
clc
close all
% Specify parameters
time_to_use     = 1; % only T48 here
merge_method    = 'concatenate'; % 'concatenate' | 'pool' % doesnt matter only 1 timepoint
pVal_thr        = 10^-6; % p-value threshold bioactivity
HCSparas = Set_Paras('pVal_thr',pVal_thr); % set parameteres
pORACLs = {'NQO1','SET','S100A11','XRCC5'}';

% Specify plot options
[targets,pathways, compound_IDs]=getselleckpathways();
uniquepathways_S = table2cell(unique(pathways));
uniquetargets_S = table2cell(unique(targets));
plot_opt_query = {'cate_to_show',horzcat(uniquepathways_S',{'DMSO'})};
plot_opt_reference = {'cate_to_show',{'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR'}}; %plot_opt_reference = {'cate_to_show',{'DMSO','Proteasome','mTOR'}};
plot_opt_ALL = {'cate_to_show',horzcat(plot_opt_query{1,2},plot_opt_reference{1,2})};

% Select subset of data to analyze
      
subset = 'all'; % query % all
switch subset
    case 'reference'
        plot_opt = plot_opt_reference;
        for c = 1:4
            pORACL = pORACLs{c};
            plates = currexp(currexp.CellLine== strcat('A549_',pORACL) & currexp.plate_type == 'reference_plate',:); % select plateset %plates = currexp(currexp.CellLine== pORACL & currexp.plate_type == 'reference_plate',:);
            plate_IDs = plates.expt_plate;
            
            plate = Load_Batch_LE(plate_IDs,profiles_folder);
            checknumplatewells()
            plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
            [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
            is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
            plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
            bioactive = Plate2Table(plate3);
            %writetable(bioactive, fullfile(results_folder,strcat(pORACL,'_bioactivepvalsR.xlsx')))
            clear bioactive wellcounts
            
            Visualize_Plate_LE(plate3,plot_opt{:}); % or plot_opt_query - % can only use con_trace with just reference plates
            labelpcplot()
            title(strcat(pORACL,subset));
            saveas(gcf,strcat('../results/matlab_results/4lines/',pORACL,subset))
            %print('-f1',strcat('../results/matlab_results/4lines/',pORACL,subset),'-dpng','-r500');
            close all
        end
        
    case 'query'
        plot_opt = plot_opt_query;
        for c = 1:4
            pORACL = pORACLs{c};
            plates = currexp(currexp.CellLine== strcat('A549_',pORACL),:); % select plateset %plates = currexp(currexp.CellLine== pORACL & currexp.plate_type == 'reference_plate',:);
            plate_IDs = plates.expt_plate;
            
            plate = Load_Batch_LE(plate_IDs,profiles_folder);
            checknumplatewells()
            plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
            [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
            is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
            plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
            bioactive = Plate2Table(plate3);
            %writetable(bioactive, fullfile(results_folder,strcat(pORACL,'_bioactivepvalsQ.xlsx')))
            clear bioactive wellcounts
            
            Visualize_Plate_LE(plate3,plot_opt{:}); % or plot_opt_query - % can only use con_trace with just reference plates
            labelpcplot()
            title(strcat(pORACL,subset));
            saveas(gcf,strcat('../results/matlab_results/4lines/',pORACL,subset))
            %print('-f1',strcat('../results/matlab_results/4lines/',pORACL,subset),'-dpng','-r500');
            close all
        end
        
    case 'all'
        %endrows=cell(4,1)
        plot_opt = plot_opt_ALL;
        for c = 1:4
            pORACL = pORACLs{c};
            plates = currexp(currexp.CellLine== strcat('A549_',pORACL),:); % select plateset %plates = currexp(currexp.CellLine== pORACL & currexp.plate_type == 'reference_plate',:);
            plate_IDs = plates.expt_plate;
            
            plate = Load_Batch_LE(plate_IDs,profiles_folder);
            checknumplatewells()
            plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
            [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
            is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
            plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
            bioactive = Plate2Table(plate3);
            endrows(c,1) = height(bioactive)+1;
            writetable(bioactive, fullfile(results_folder,strcat(pORACL,'_bioactivepvals.xlsx')))
            clear wellcounts
            
            Visualize_Plate_LE(plate3,plot_opt{:}); % or plot_opt_query - % can only use con_trace with just reference plates
            labelpcplot()
            title(strcat(pORACL,subset))
            saveas(gcf,strcat('../results/matlab_results/4lines/',pORACL,subset))
            %print('-f1',strcat('../results/matlab_results/4lines/',pORACL,subset),'-dpng','-r500');
            close all
        end
end


%% pvalue plots
pORACLs = {'NQO1','SET','S100A11','XRCC5'}';
%endrows = [3529;3532;3539;3529];        
NQO1_bioactive=bioactives_importfile('W:\2015_09_HTS_LE\results\matlab_results\4lines\NQO1_bioactivepvals.xlsx',1,2,endrows(1,1));
SET_bioactive=bioactives_importfile('W:\2015_09_HTS_LE\results\matlab_results\4lines\SET_bioactivepvals.xlsx',1,2,endrows(2,1));
S100A11_bioactive=bioactives_importfile('W:\2015_09_HTS_LE\results\matlab_results\4lines\S100A11_bioactivepvals.xlsx',1,2,endrows(3,1));
XRCC5_bioactive = bioactives_importfile('W:\2015_09_HTS_LE\results\matlab_results\4lines\XRCC5_bioactivepvals.xlsx',1,2,endrows(4,1));

pvalue_plotter % plot P values by cell line, and by batch

%% Cross validation for reference plates
% Specify parameters
time_to_use         = 1;
merge_mode          = 'concatenate'; % 'concatenate' | 'pool'
pVal_thr            = 10^-6; % p-value threshold for calling bioactive compounds
classifier_type     = 'KLFDA_knn'; % 'KLFDA_knn' | 'KLFDA_nearest_centroid' % means use KNN on LDA space
NumNeighbors        = 1; % Number of neighbors for knn classifier, only be used for KLFDA_knn
FDA_type            = 0; % 0: FDA, -1: no transformation, -2: RCA
g                   = 0:0.1:1; % regularization strength for within-class variance
kernel_type         = 'linear'; % linear rbf polynomial
kernel_para         = 0; % no effect if using linear kernel
kfold               = 10;
paras = Set_Paras('pVal_thr',pVal_thr);

% Pick out AW reference
aw_ref = {'DMSO','Actin','AuroraB','DNA',...
    'ER','HDAC','HSP90','MT','PLK','Proteasome',...
    'mTOR'};
pORACL = 'XRCC5';
switch pORACL
    case 'NQO1'
        %Get plates
        HCSparas.path.features = 'features/cbfeatures/cbfeatures-'; % XRCC5 has different location for S in this case
        currexp.cpd_plate_1 = categorical(currexp.cpd_plate_1);
        % plates should only be reference plates
        currexp.plate_type = categorical(currexp.plate_type);
        plates = currexp(currexp.plate_type=='reference_plate' & currexp.CellLine=='A549_NQO1' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7' & currexp.cpd_plate_1~='ECHO REFERENCE',:);
        plate_IDs = plates.expt_plate;
        plate = Load_Batch_LE(plate_IDs,profiles_folder);
        plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);
        plate2 = Get_Bioactive_Compounds(plate2, paras);
        % Get reference compounds
        plate2 = Select_Rows(plate2,ismember(plate2.drug_categories,aw_ref));
        % Prepare training data
        trainData = cell2mat(plate2.profiles);
        trainLabels = plate2.drug_categories;
        trainDataInfo = strcat(plate2.drug_names,'_',cellstr(num2str(plate2.concentrations)));
        % Cross-validation
        [clfy_obj, accuracy, cv_report, scan_scheme] = Classifier(trainData, trainLabels, ...
            classifier_type,'XtrainNames',trainDataInfo,'cv_fold',kfold,...
            'FDA_type',FDA_type,'gamma',g,...
            'kernel_type',kernel_type,'kernel_para',kernel_para,...
            'NumNeighbors',NumNeighbors);
        % Visualize the FDA space
        plate_fda = plate2;
        plate_fda.profiles = cellfun(@clfy_obj.KLFDA_proj_fun,plate2.profiles,'unif',false);
        Visualize_Plate(plate_fda,plot_opt_reference{:});
        % PREDICTIONS
        % Specify parameters
        g                   = 0:0.05:0.5; % regularization strength for within-class variance
        classifier_paras    = {'NumNeighbors',1}; % e.g. 'NumNeighbors',5
        HCSparas = Set_Paras('pVal_thr',pVal_thr,'classifier_type',classifier_type,...
            'FDA_type',FDA_type,...
            'gamma',g,'kernel_type',kernel_type,'kernel_para',kernel_para,...
            'classifier_paras',classifier_paras,'kfold',kfold);
        % Predict Output
        predict_output = Screen_Plates( plate, HCSparas);
end

%% Clustering
plates = currexp(currexp.CellLine=='A549_XRCC5' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7' & currexp.cpd_plate_1~='ECHO REFERENCE',:);
plate_IDs = plates.expt_plate;
batches{1,1} = plate_IDs;
% Set parameters
FDA_type            = -1; % 0: FDA, -1: no transformation, -2: RCA
g                   = 0; % -1: use theoretical value to regularize the within-class Variance matrix for LDA
linkage_method      = 'average'; % 'average' | 'centroid' | 'complete' | 'single'
D_within_prob_thr   = 0.5; % probability cutoff for choosing diatance threshold to define clusters
col_order = {'Cluster','Size','Significance','plateIDs','well_names','drug_names','drug_categories',...
    'concentrations','bioactive_pVals','targets','pathways','descriptions','cas','compound',...
    'cpd_usage','nCells','clone_names','plate_types'};
sort_order = {'Size','Cluster','pathways','targets','drug_names','concentrations'};
%=========================================================================================
paras = Set_Paras('pVal_thr',pVal_thr,'classifier_type',classifier_type,...
    'FDA_type',FDA_type,...
    'gamma',g,'kernel_type',kernel_type,'kernel_para',kernel_para,...
    'classifier_paras',classifier_paras,'kfold',kfold,...
    'linkage_method',linkage_method,'D_within_prob_thr',D_within_prob_thr);
Nbatches = length(batches);
cluster_info = cell(Nbatches,1);
profiles = cell(Nbatches,1);
tree = cell(Nbatches,1);
for b = 1:Nbatches
    disp(b)
    % 1. Get plates
    plate_to_cluster = Load_Batch(batches{b},profiles_folder);
    
    % 2. Clustering
    [cluster_info{b}, profiles{b}, tree{b}] = Cluster_Compounds(plate_to_cluster,paras);
    
    % 3. Add compound names for L1700-Selleck-Bioactive-1836-SMDC
    %[is_found,loc] = ismember(cluster_info{b}.drug_names,smdc_lib.CatalogNo_);
    cluster_info{b}.compound = cluster_info{b}.drug_names;
    %cluster_info{b}.compound(is_found) = smdc_lib.ProductName(loc(is_found));
    
    cluster_info{b} = cluster_info{b}(:,col_order);
    [cluster_info{b},idx] = sortrows(cluster_info{b},sort_order);
    profiles{b} = profiles{b}(idx,:);
    tree{b} = Update_NodeIDs(tree{b},idx);
end
c2 = sortrows(cluster_info{b},{'pathways','drug_categories','drug_names','targets','Cluster','plateIDs'});
c3 = sortrows(cluster_info{b},{'drug_categories','Cluster','plateIDs','drug_names','concentrations'});
disp('Done')

%% Visualize clusters
batch_to_plot = 1;
clusters_to_color = [];%[13,38,80,83,86,87];
if isempty(clusters_to_color)
    %------ get all significant clusters ------%
    sig_thr = 0.05;
    sig_clusters = unique(cluster_info{batch_to_plot}.Cluster(...
        cluster_info{batch_to_plot}.Significance<sig_thr));
    %------ get DMSO cluster(s) ------%
    [m,n,g] = grpstats(strcmp(cluster_info{batch_to_plot}.drug_categories,'DMSO'),...
        cluster_info{batch_to_plot}.Cluster,{'mean','numel','gname'});
    % DMSO cluster must contain more than 10 compounds and more than 50% of them are DMSO
    dmso_cluster = str2double(g( (m>0.5) & (n>10) ));
    % color all significant clusters and DMSO cluster(s)
    clusters_to_color = [sig_clusters(:);dmso_cluster(:)];
end
PlotPhylogeneticTree(cluster_info{batch_to_plot},tree{batch_to_plot},clusters_to_color);
%% Analysis of clusters
fh = figure;
for i = 1:length(cluster_info)
    ref_cluster(i).categories = unique(cluster_info{i}.drug_categories);
    ref_cluster(i).cluster_ID = unique(cluster_info{i}.Cluster);
    
    ref_cluster(i).stat = zeros(length(ref_cluster(i).categories),length(ref_cluster(i).cluster_ID));
    for j = 1:size(ref_cluster(i).stat,1)
        is_class = strcmpi(cluster_info{i}.drug_categories,ref_cluster(i).categories{j});
        tb = tabulate(cluster_info{i}.Cluster(is_class));
        
        [is_in,loc] = ismember(tb(:,1),ref_cluster(i).cluster_ID);
        ref_cluster(i).stat(j,loc(is_in)) = tb(:,2);
    end
    
    ax(i) = subplot(2,3,i,'Parent',fh);
    imagesc(ax(i),ref_cluster(i).stat),colormap(ax(i),hot)%,caxis(ax(i),[0,1])
    ax(i).YTick = 1:length(ref_cluster(i).categories);
    ax(i).YTickLabel = ref_cluster(i).categories;
    ax(i).TickLabelInterpreter = 'none';
    
    %     ax(i).XTick = 1:length(ref_cluster(i).cluster_ID);
    %     ax(i).XTickLabel = GN_con;
    title(num2str(i))
end