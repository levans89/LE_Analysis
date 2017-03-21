%% Set paths
%restoredefaultpath;
addpath(genpath('W:/Lab_Analysis/common/plate_annotation'))
addpath(genpath('W:\2015_09_HTS_LE\Code\CD_Tag_Drug_Screen'))
save_destination = 'W:\2015_09_HTS_LE\data\profiles';% place to save profiles to
results_folder = 'W:\2015_09_HTS_LE\results\matlab_results\4lines';
profile_folder = 'W:\2015_09_HTS_LE\data\profiles';% place to load profiles from
load('currexp.mat');% load in currexp.mat which details all plates in this screen (extracted from exp_DB)
currexp.CellLine = categorical(currexp.CellLine);
currexp.batch = categorical(currexp.batch);
qc_folder = 'W:/2015_09_HTS_LE/QC/4lines/';

%% Profile a batch of plates Charles
%%% Set parameters %%%
pORACL = 'NQO1'; % we screened 4 lines - here we switch profile generation between the lines to deal with marker encoding in feature names
switch pORACL %potential ORACL
    case 'NQO1'
        features_to_use = feature_selector('NBT_NQO1'); %{'all'};%Get_Default_Features('LE_147_to_150_fixed');
        plates = currexp(currexp.CellLine=='A549_NQO1' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plateIDs = plates.expt_plate;
        bad_wells = LE_Bad_Wells(plates,qc_folder); % Get all bad wells
        writetable(bad_wells,'W:\2015_09_HTS_LE\QC\4lines\bw_NQO1.xlsx');
        
        for m = 1:height(plates)
            plate_maps{m} = [num2str(plates.expt_plate{m}),'.xlsx'];
        end
        plate_types = plates.plate_type;
        paras = Set_Paras();
        plate_maps_foler = paras.path.plate_maps;
        
        h = waitbar(0,num2str(0));
        for p = 1:length(plateIDs)
            waitbar(p/length(plateIDs),h,[num2str(p),'/',num2str(length(plateIDs))])
            
            % 1. Get all inputs for Get_KS_Profiles
            plateID = plateIDs{p};
            plate_type = plate_types{p};
            plate_map_file = [plate_maps_foler filesep plate_maps{p}];
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
                'save_destination',save_destination,'features_to_use',features_to_use);
            
        end
        waitbar(p/length(plateIDs),h,'Done')
        delete(h)
      
    case 'XRCC5'%%% S / _Alice method
        features_to_use = feature_selector('NBT_XRCC5'); %{'all'};%Get_Default_Features('LE_147_to_150_fixed');
        plates = currexp(currexp.CellLine=='A549_XRCC5' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plateIDs = plates.expt_plate;
        bad_wells = LE_Bad_Wells(plates,qc_folder); % Get all bad wells
        writetable(bad_wells,'W:\2015_09_HTS_LE\QC\4lines\bw_XRCC5.xlsx');
        
        for m = 1:height(plates)
            plate_maps{m} = [num2str(plates.expt_plate{m}),'.xlsx'];
        end
        plate_types = plates.plate_type;
        paras = Set_Paras();
        paras.path.features = ['features' filesep 'cbfeatures_Alice' filesep 'cbfeatures-']; % switch feature folder here to get S method
        plate_maps_foler = paras.path.plate_maps;
        
        h = waitbar(0,num2str(0));
        for p = 2:length(plateIDs)
            waitbar(p/length(plateIDs),h,[num2str(p),'/',num2str(length(plateIDs))])
            
            % 1. Get all inputs for Get_KS_Profiles
            plateID = plateIDs{p};
            plate_type = plate_types{p};
            plate_map_file = [plate_maps_foler filesep plate_maps{p}];
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
                'save_destination',save_destination,'features_to_use',features_to_use);
            
        end
        waitbar(p/length(plateIDs),h,'Done')
        delete(h)
        
    case 'SET'
        features_to_use = feature_selector('NBT_SET'); %{'all'};%Get_Default_Features('LE_147_to_150_fixed');
        plates = currexp(currexp.CellLine=='A549_SET' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plateIDs = plates.expt_plate;
        bad_wells = LE_Bad_Wells(plates,qc_folder); % Get all bad wells
        writetable(bad_wells,'W:\2015_09_HTS_LE\QC\4lines\bw_SET.xlsx');
        
        for m = 1:height(plates)
            plate_maps{m} = [num2str(plates.expt_plate{m}),'.xlsx'];
        end
        plate_types = plates.plate_type;
        paras = Set_Paras();
        plate_maps_foler = paras.path.plate_maps;
        
        h = waitbar(0,num2str(0));
        for p = 1:length(plateIDs)
            waitbar(p/length(plateIDs),h,[num2str(p),'/',num2str(length(plateIDs))])
            
            % 1. Get all inputs for Get_KS_Profiles
            plateID = plateIDs{p};
            plate_type = plate_types{p};
            plate_map_file = [plate_maps_foler filesep plate_maps{p}];
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
                'save_destination',save_destination,'features_to_use',features_to_use);
            
        end
        waitbar(p/length(plateIDs),h,'Done')
        delete(h)
        
    case 'S100A11'
        features_to_use = feature_selector('NBT_S100A11'); %{'all'};%Get_Default_Features('LE_147_to_150_fixed');
        plates = currexp(currexp.CellLine=='A549_S100A11' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plateIDs = plates.expt_plate;
        bad_wells = LE_Bad_Wells(plates,qc_folder); % Get all bad wells
        writetable(bad_wells,'W:\2015_09_HTS_LE\QC\4lines\bw_S100A11.xlsx');
        
        for m = 1:height(plates)
            plate_maps{m} = [num2str(plates.expt_plate{m}),'.xlsx'];
        end
        plate_types = plates.plate_type;
        paras = Set_Paras();
        plate_maps_foler = paras.path.plate_maps;
        
        h = waitbar(0,num2str(0));
        for p = 1:length(plateIDs)
            waitbar(p/length(plateIDs),h,[num2str(p),'/',num2str(length(plateIDs))])
            
            % 1. Get all inputs for Get_KS_Profiles
            plateID = plateIDs{p};
            plate_type = plate_types{p};
            plate_map_file = [plate_maps_foler filesep plate_maps{p}];
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
                'save_destination',save_destination,'features_to_use',features_to_use);
            
        end
        waitbar(p/length(plateIDs),h,'Done')
        delete(h)
end
%% Check cell number
clc
cd ('W:\2015_09_HTS_LE\QC\4lines\')
for b = 55:132
pid = [b]; % first check blanks, ref plates, compound plates
plateIDs = arrayfun(@(x)sprintf('2017018%03d',x),pid,'unif',false)'; %appears that only reference plates can be checked but not true
try
plate = Load_Batch(plateIDs,profile_folder); 
Show_Cell_Number(plate);
saveas(gcf,plateIDs{1,1},'bmp');
catch    
end
end
close all
cd ('W:\2015_09_HTS_LE\Code\')
%pp = Plate2PhysicalPlate(plate);% Convert a plate from a column-oriented data format to a physical plate format.
%% Visualize plates
% Specify parameters
time_to_use     = 1;
merge_method    = 'concatenate'; % 'concatenate' | 'pool'
pVal_thr        = 10^-6; % p-value threshold bioactivity
plot_opt_reference = {'cate_to_show',{'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR'}};
[targets,pathways, compound_IDs]=getselleckpathways();
uniquepathways = table2cell(unique(pathways));
plot_opt_query = {'cate_to_show',uniquepathways'};
HCSparas = Set_Paras('pVal_thr',pVal_thr);

% to_use = (~strcmp(plate.cpd_usage,'')|strcmp(plate.pathways,'Apoptosis'));
% plate = Select_Rows(plate,to_use);

pORACL = 'NQO1'
switch pORALC
    case 'NQO1'
        HCSparas.path.features = 'features\cbfeatures\cbfeatures_-'
        plates = currexp(currexp.CellLine=='A549_NQO1' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plate_IDs = plates.expt_plate;
        plate = Load_Batch(plate_IDs,profile_folder);
        plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
        [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
        is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
        plate3.drug_categories(is_inactive) = {'Nonbioactive'};
        Visualize_Plate(plate3,plot_opt_reference{:}) % or plot_opt_query
        NQO1_bioactive = Plate2Table(plate3);
        writetable(NQO1_bioactive, fullfile(results_folder,'NQO1_bioactive.xlsx'))
        
    case 'XRCC5'
        HCSparas.path.features = 'features\cbfeatures\cbfeatures_Alice-' % XRCC5
        currexp.cpd_plate_1 = categorical(currexp.cpd_plate_1);
        plates = currexp(currexp.CellLine=='A549_XRCC5' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7' & currexp.cpd_plate_1~='ECHO REFERENCE',:);
        plate_IDs = plates.expt_plate;
        plate = Load_Batch(plate_IDs,profile_folder);
        plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
        [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
        is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
        plate3.drug_categories(is_inactive) = {'Nonbioactive'};
        Visualize_Plate(plate3,plot_opt_reference{:})
        XRCC5_bioactive = Plate2Table(plate3);
        writetable(XRCC5_bioactive, fullfile(results_folder,'XRCC5_bioactive.xlsx'));
        
    case 'SET'
        HCSparas.path.features = 'features\cbfeatures\cbfeatures_-'
        plates = currexp(currexp.CellLine=='A549_SET' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plate_IDs = plates.expt_plate;
        plate = Load_Batch(plate_IDs,profile_folder);
        plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
        [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
        is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
        plate3.drug_categories(is_inactive) = {'Nonbioactive'};
        Visualize_Plate(plate3,plot_opt_reference{:})
        SET_bioactive = Plate2Table(plate3);
        writetable(SET_bioactive, fullfile(results_folder,'SET_bioactive.xlsx'))
        
    case 'S100A11'
        HCSparas.path.features = 'features\cbfeatures\cbfeatures_-'
        plates = currexp(currexp.CellLine=='A549_S100A11' & currexp.batch~='4L2K_8' & currexp.batch~='4L2K_7',:);
        plate_IDs = plates.expt_plate;
        plate = Load_Batch(plate_IDs,profile_folder);
        plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
        [plate2,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
        is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
        plate3.drug_categories(is_inactive) = {'Nonbioactive'};
        Visualize_Plate(plate3,plot_opt_reference{:})
        S100A11_bioactive = Plate2Table(plate3);
        writetable(S100A11_bioactive, fullfile(results_folder,'S100A11_bioactive.xlsx'))
        
end
%% Cross validation for reference plates
% clc
% Specify parameters
%pid                 = [59]; % plate ID, last 2 digits
time_to_use         = 1:2;
merge_mode          = 'concatenate'; % 'concatenate' | 'pool'
pVal_thr            = 10^-2; % p-value threshold for calling bioactive compounds
classifier_type     = 'KLFDA_knn'; % 'KLFDA_knn' | 'KLFDA_nearest_centroid'
NumNeighbors        = 1; % Number of neighbors for knn classifier, only be used for KLFDA_knn
FDA_type            = 0; % 0: FDA, -1: no transformation, -2: RCA 
g                   = 0:0.1:1; % regularization strength for within-class variance
kernel_type         = 'linear'; % linear rbf polynomial
kernel_para         = 0; % no effect if using linear kernel
kfold               = 10;
plot_opt_reference            = {};
%=========================================================================================
paras = Set_Paras('pVal_thr',pVal_thr);
% Get plates and call bioactive compounds
plate_IDs = arrayfun(@(x)sprintf('2017018%03d',x),pid,'unif',false)'; 
plate = Load_Batch(plate_IDs,profile_folder); 
plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);
plate2 = Get_Bioactive_Compounds(plate2, paras);
% Pick out AW reference
aw_ref = {'DMSO','Actin','AuroraB','DNA',...
          'ER','HDAC','HSP90','MT','PLK','Proteasome',...
          'mTOR'};
plate2 = Select_Rows(plate2,ismember(plate2.drug_categories,aw_ref));
% Prepare training data
trainData = cell2mat(plate2.profiles);
trainLabels = plate2.drug_categories;
trainDataInfo = strcat(plate2.drug_names,'_',cellstr(num2str(plate2.concentrations)));
% Cross-validation
% rng(1065551) % for reproducibility
[clfy_obj, accuracy, cv_report, scan_scheme] = Classifier( trainData, trainLabels, ...
                       classifier_type,'XtrainNames',trainDataInfo,'cv_fold',kfold,...
                       'FDA_type',FDA_type,'gamma',g,...
                       'kernel_type',kernel_type,'kernel_para',kernel_para,...
                       'NumNeighbors',NumNeighbors);
% Visualize the FDA space
plate_fda = plate2;
plate_fda.profiles = cellfun(@clfy_obj.KLFDA_proj_fun,plate2.profiles,'unif',false);
Visualize_Plate(plate_fda,plot_opt_reference{:});
%% Make Predictions
clc
%-----------------------------------------------------------------------------------------
% Specify parameters
pids                = [59]; % plate ID, last 2 digits 
pVal_thr            = 10^-6; % p-value threshold for calling bioactive compounds
classifier_type     = 'KLFDA_knn'; % 'KLFDA_knn' | 'KLFDA_nearest_centroid'
FDA_type            = 0; % 0: FDA, -1: no transformation, -2: RCA 
g                   = 0:0.05:0.5; % regularization strength for within-class variance
kernel_type         = 'linear'; % linear, rbf, polynomial
kernel_para         = 0; % no effect if using linear kernel
classifier_paras    = {'NumNeighbors',1}; % e.g. 'NumNeighbors',5
kfold               = 10; % fold of cross validation when training the classifier
%=========================================================================================
HCSparas = Set_Paras('pVal_thr',pVal_thr,'classifier_type',classifier_type,...
                     'FDA_type',FDA_type,...
                     'gamma',g,'kernel_type',kernel_type,'kernel_para',kernel_para,...
                     'classifier_paras',classifier_paras,'kfold',kfold);
plateIDs = arrayfun(@(x)sprintf('2017018%03d',x),pids,'unif',false)'; 
plate = Load_Batch(plateIDs,profile_folder);
predict_output = Screen_Plates( plate, HCSparas);
%% Clustering
pids=[];
batches = {arrayfun(@(x)sprintf('2017018%03d',x),pids,'unif',false)'};
%-----------------------------------------------------------------------------------------
% Set parameters
pVal_thr            = 10^-3; % p-value threshold for calling bioactive compounds
classifier_type     = 'KLFDA_knn'; % 'KLFDA_knn' | 'KLFDA_nearest_centroid'
FDA_type            = -1; % 0: FDA, -1: no transformation, -2: RCA 
g                   = 0; % -1: use theoretical value to regularize the within-class Variance matrix for LDA
kernel_type         = 'linear';
kernel_para         = 0;
classifier_paras    = {'NumNeighbors',1};
kfold               = 10;
linkage_method      = 'average'; % 'average' | 'centroid' | 'complete' | 'single'
D_within_prob_thr   = 0.5; % probability cutoff for choosing diatance threshold to define clusters
%-----------------------------------------------------------------------------------------
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
    plate_to_cluster = Load_Batch(batches{b},profile_folder);
    
    % 2. Clustering
    [cluster_info{b}, profiles{b}, tree{b}] = Cluster_Compounds(plate_to_cluster,paras);
    
    % 3. Add compound names for L1700-Selleck-Bioactive-1836-SMDC
    %[is_found,loc] = ismember(cluster_info{b}.drug_names,smdc_lib.CatalogNo_);
    cluster_info{b}.compound = cluster_info{b}.drug_names;
    cluster_info{b}.compound(is_found) = smdc_lib.ProductName(loc(is_found));
    
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
clusters_to_color = []%[13,38,80,83,86,87];
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