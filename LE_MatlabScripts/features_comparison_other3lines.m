% Clean Features_comparison
% Louise Heinrich July 2017 in Altschuler & Wu Lab
% Code to determine classification accuracy of pORACL with different feature
% sets
%% Set Paths
addpath(genpath('W:\Lab_Analysis\common\plate_annotation')) %plate annotation from labtoolbox
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\CD_Tag_Drug_Screen')) % Charles code for NBT paper containing useful functions
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_analysis')) % Louise code afer cluster operations
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_Preprocessing')) % Louise code before cluster operations
addpath(genpath('W:\Lab_Analysis\common\image_browser')) % Image_Browser to callback images

%% Define Variables
code_folder = 'W:\2015_09_HTS_LE\Code_LE\';
profiles_folder = 'W:\2015_09_HTS_LE\data\profiles';% place to save profiles to and load from
qc_folder = 'W:\2015_09_HTS_LE\QC\bigscreen';% store any QC analysis Python, Matlab, or manual
results_folder = 'W:\2015_09_HTS_LE\results\feature_sets\20170901\'; % Louise results for 4lines analysis
plate_DB_path ='W:\2015_09_HTS_LE\project_database\'; % directory containing plate database file
plate_DB = 'Plate_database_latest2.xlsx'; % file describing relationship between compound and experiment plates

% Load Constants
Exp_DB = readtable(fullfile(plate_DB_path,plate_DB),'Sheet',1); % full table of all plates

year = 2017;
theseyearplates = ismember(Exp_DB.year,year);
plates = Exp_DB(theseyearplates,:); % subset 2017 plates

batch='4L2K_8';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:); % subset 4lines2K plates b1

batch='4L2K_7';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);% subset 4lines2K plates b2

cellplates = plates(strcmpi(plates.plate_type,'reference_plate'),:); %select only ref plates
cellplates = cellplates(strcmpi(cellplates.Experiment_,'LE_4'),:);

cellplates.CellLine = categorical(cellplates.CellLine);
clear thesebatchplates theseyearplates batch Exp_DB plates_a plates_b plates_c year plate_DB_path plate_DB

% Get Analysis Parameters
% Get plotting options
getplotoptions() % generate standard plot options
plot_opt = plot_opt_reference; % specific plot option
clear plot_opt_query plot_opt_ALL plot_opt_reference compound_IDs pathways targets uniquepathways_S uniquetargets_S

%bioactive_paras
pVal_thr = 10^-6; %10^-6
max_nCtrl = 2500; %2500
var_pct_to_keep = 0.95; %0.95
k_fold = 10; %10

%screen_paras
time_to_use = 1; %1 | 2 | 'all'
merge_mode = 'concatenate';% 'concatenate' | 'pool'
classifier_type = 'KLFDA_knn'; % 'KLFDA_knn' | 'KLFDA_nearest_centroid'
FDA_type = 0; % 0: FDA | -1: no transformation | -2: RCA
gamma = 0:0.1:1; % regularization strength for within-class variance
kernel_type = 'linear';
kernel_para = 0;
classifier_paras = {'NumNeighbors',3};
kfold = 10;
conf_thr = 0.1;
do_view = true;
KNN_to_list = 10;

%cluster_paras
linkage_method = 'average';
ref_flag = {'reference_cpd','negative_ctrl'};
nDis_max = 50000;
D_within_prob_thr= 0.9;
min_dmso_pct_to_be_DMSO_cluster = 0.5;
min_dmso_number_to_be_DMSO_cluster = 10;

% Set Parameters
paras = Set_Paras('pVal_thr',pVal_thr,'max_nCtrl',max_nCtrl,'var_pct_to_keep',var_pct_to_keep,...
    'k_fold',k_fold,'time_to_use',time_to_use,'merge_mode',merge_mode,'classifier_type',classifier_type,...
    'FDA_type',FDA_type,'gamma',gamma,'kernel_type',kernel_type,'kernel_para',kernel_para,...
    'classifier_paras',classifier_paras,'kfold',kfold,'conf_thr',conf_thr,'do_view',do_view,...
    'linkage_method',linkage_method,'ref_flag',ref_flag,'nDis_max',nDis_max,'D_within_prob_thr',D_within_prob_thr,...
    'min_dmso_pct_to_be_DMSO_cluster',min_dmso_pct_to_be_DMSO_cluster,'min_dmso_number_to_be_DMSO_cluster',min_dmso_number_to_be_DMSO_cluster, 'KNN_to_list',KNN_to_list);% set parameters

%% Generate Full Profiles - all features all channels - subset from full profiles internally in code
pORACL = {'A549_XRCC5','A549_NQO1','A549_SET','A549_S100A11'};
for o = 1:4
    oraclplates = cellplates(cellplates.CellLine == pORACL{o},:);
    LE_Gen_Profiles(oraclplates,pORACL{o},profiles_folder)
end
%% Import and subset features
% Import indexing of features to select
[name,category_id,channel_id,NBT,use_2ch,use_all]= import_featuresselection2('W:\2015_09_HTS_LE\Code_LE\LE_Analysis\LE_MatlabScripts\feature_select_NBT.txt'); % manually indexed feature list based on NBT feature set from supplementary tables
feature_table = table(name,category_id,channel_id,NBT,use_2ch,use_all); % comparison of possible feature subsets in table format
writetable(feature_table,[results_folder,'\feature_table.xlsx'])

% create logicals for subsetting features
NBT = NBT==1; % NBT set 234 features focusing on YFP channel
NBT = NBT';
use_2ch = use_2ch==1; % Ch1 Ch3 443 features not using any YFP features
use_2ch = use_2ch';
use_all = use_all==1; % Ch1 Ch2 Ch3 659 all features all channels
use_all = use_all';

clear name category_id channel_id w

%% Classification Accuracy
pORACL = {'A549_NQO1','A549_SET','A549_S100A11','A549_XRCC5'};
% Pick 6 class or 10 class model for classification accuracy
%aw_ref = {'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR'};
aw_ref = {'DMSO','DNA','HDAC','HSP90','MT','Proteasome','mTOR'};

% Feature Set Selection
featureset ={NBT,use_2ch,use_all}; % use_2ch use_all NBT
featurenames = {'NBT','use_2ch','use_all'};

for f = 1:3
    for o = 1:4 
        cellplates2 = cellplates(cellplates.CellLine == pORACL{o},:);
        plateIDs = cellplates2.expt_plate;
        % Import profiles
        plate = Load_Batch_LE(cellplates2,plateIDs, profiles_folder,qc_folder);
        % Subset features
        plate.feaInfo.name = plate.feaInfo.name(featureset{f});
        plate.feaInfo.channel_id = plate.feaInfo.channel_id(featureset{f});
        plate.feaInfo.category_id = plate.feaInfo.category_id(featureset{f});
        for w = 1:size(plate.profiles,1)
            plate.profiles{w,1} = plate.profiles{w,1}(featureset{f});
        end
        % Calc CA % for all ref plates pooled together
        plate2 = Merge_Time_Points(plate,paras.screen.time_to_use,paras.screen.merge_mode); % legacy
        [~,plate3,~,~,~,~] = Get_Bioactive_Compounds(plate2, paras);% Charles function - bioactivity
        is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO')); %inactive cpds
        plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
        bioactive_use_2ch = Plate2Table(plate3);% make table
        clear plate2 is_inactive
        % Get reference compounds
        plate3 = Select_Rows(plate3,ismember(plate3.drug_categories,aw_ref));
        plate3 = Select_Rows(plate3,ismember(plate3.plate_types,{'reference_plate'}));
        % Prepare training data
        trainData = cell2mat(plate3.profiles);
        trainLabels = plate3.drug_categories;
        trainDataInfo = strcat(plate3.drug_names,'_',cellstr(num2str(plate3.concentrations)));
        % Cross-validation
        [clfy_obj, accuracy, cv_report, scan_scheme] = Classifier(trainData, trainLabels,classifier_type,'XtrainNames',trainDataInfo,'cv_fold',kfold,...
            'FDA_type',FDA_type,'gamma',gamma,'kernel_type',kernel_type,'kernel_para',kernel_para,...
            'NumNeighbors',3);
        classifier.accuracy{o} = accuracy;
        classifier.clfy_obj{o} = clfy_obj;
        classifier.cv_report{o} = cv_report;
        classifier.scan_scheme{o} = scan_scheme;
        save([results_folder,pORACL{o},featurenames{f},'classifier.mat'],'classifier')
        clear accuracy clfy_obj cv_report scan_scheme
        
        % Calculate per plate CA%
        for p=1:height(cellplates2)
            % Load single plate
            plate_IDs = cellplates2.expt_plate;
            plate_ID = plate_IDs(p);
            display(plate_ID);
            plate = Load_Batch_LE(cellplates2,plate_ID, profiles_folder,qc_folder);
            % subset single plate
            plate.feaInfo.name = plate.feaInfo.name(featureset{f});
            plate.feaInfo.channel_id = plate.feaInfo.channel_id(featureset{f});
            plate.feaInfo.category_id = plate.feaInfo.category_id(featureset{f});
            for w = 1:size(plate.profiles,1)
                plate.profiles{w,1} = plate.profiles{w,1}(featureset{f});
            end
            % Calc CA%
            plate2 = Merge_Time_Points(plate,paras.screen.time_to_use,paras.screen.merge_mode);
            plate2 = Get_Bioactive_Compounds(plate2, paras);
            % Get reference compounds
            plate2 = Select_Rows(plate2,ismember(plate2.drug_categories,aw_ref));
            plate2 = Select_Rows(plate2,ismember(plate2.plate_types,{'reference_plate'}));
            % Prepare training data
            trainData = cell2mat(plate2.profiles);
            trainLabels = plate2.drug_categories;
            trainDataInfo = strcat(plate2.drug_names,'_',cellstr(num2str(plate2.concentrations)));
            % Cross-validation
            [clfy_obj, accuracy, cv_report, scan_scheme] = Classifier(trainData, trainLabels,classifier_type,'XtrainNames',trainDataInfo,'cv_fold',kfold,...
                'FDA_type',FDA_type,'gamma',gamma,'kernel_type',kernel_type,'kernel_para',kernel_para,...
                'NumNeighbors',3);
            % Record Results
            classifier_pp.accuracy{p,o} = accuracy;
            classifier_pp.clfy_obj{p,o} = clfy_obj;
            classifier_pp.cv_report{p,o} = cv_report;
            classifier_pp.scan_scheme{p,o} = scan_scheme;
        end
        save([results_folder,pORACL{o},featurenames{f},'classifier_pp.mat'],'classifier_pp')
        
        % Calc CA for any 2 plates pooled (so example training data in dose/class)
        % define all combos lof 6 plates to pool
        combos = cell(15,1);
        for x = 1:5
            combos{x,1} = cellplates2([1 x+1],{'expt_plate','plate_type'});
        end
        for y = 3:6
            combos{y+3,1} = cellplates2([2 y],{'expt_plate','plate_type'});
        end
        for z = 4:6
            combos{z+6,1} = cellplates2([3 z],{'expt_plate','plate_type'});
        end
        for w = 5:6
            combos{w+8,1} = cellplates2([4 w],{'expt_plate','plate_type'});
        end
        combos{15,1} = cellplates2([5 6],{'expt_plate','plate_type'});
        
        % Calc CA% for pairwise combinations of reference plates
        for c=1:length(combos)
            cellplates2 = combos{c};
            plateIDs = cellplates2.expt_plate;
            display(plateIDs)
            % Load 2 selected plates from 6
            plate = Load_Batch_LE(cellplates,plateIDs, profiles_folder,qc_folder);
            % subset
            plate.feaInfo.name = plate.feaInfo.name(featureset{f});
            plate.feaInfo.channel_id = plate.feaInfo.channel_id(featureset{f});
            plate.feaInfo.category_id = plate.feaInfo.category_id(featureset{f});
            for w = 1:size(plate.profiles,1)
                plate.profiles{w,1} = plate.profiles{w,1}(featureset{f});
            end
            % calc bioactivity
            plate2 = Merge_Time_Points(plate,paras.screen.time_to_use,paras.screen.merge_mode);
            [~,plate3,~,~,~,~] = Get_Bioactive_Compounds(plate2, paras);
            is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
            plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
            bioactive_use_2ch = Plate2Table(plate3);
            clear plate2 is_inactive
            % Get reference compounds
            plate3 = Select_Rows(plate3,ismember(plate3.drug_categories,aw_ref));
            plate3 = Select_Rows(plate3,ismember(plate3.plate_types,{'reference_plate'}));
            % Prepare training data
            trainData = cell2mat(plate3.profiles);
            trainLabels = plate3.drug_categories;
            trainDataInfo = strcat(plate3.drug_names,'_',cellstr(num2str(plate3.concentrations)));
            % Cross-validation
            [clfy_obj, accuracy, cv_report, scan_scheme] = Classifier(trainData, trainLabels,classifier_type,'XtrainNames',trainDataInfo,'cv_fold',kfold,...
                'FDA_type',FDA_type,'gamma',gamma,'kernel_type',kernel_type,'kernel_para',kernel_para,...
                'NumNeighbors',3);
            classifier_pw.accuracy{c,o} = accuracy;
            classifier_pw.clfy_obj{c,o} = clfy_obj;
            classifier_pw.cv_report{c,o} = cv_report;
            classifier_pw.scan_scheme{c,o} = scan_scheme;
        end
        save([results_folder,pORACL{o},featurenames{f},'classifier_pw.mat'],'classifier_pw')
        
    end
end
%% Visualisation in PC space
% Choose class model
aw_ref = {'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR','Unknown'};
%aw_ref = {'DMSO','DNA','HDAC','HSP90','MT','Proteasome','mTOR'};

featureset{f} = use_all; %NBT %use_2ch %use_all


% Bring back full plateset
plates.plate_type = categorical(plates.plate_type);
cellplates = plates(plates.plate_type == 'reference_plate',:); %select only ref plates
cellplates.CellLine = categorical(cellplates.CellLine);
cellplates = cellplates(cellplates.CellLine == 'A549_NQO1',:); %try on one cell line
plateIDs = cellplates.expt_plate;

% Load all plates
plate = Load_Batch_LE(cellplates,plateIDs, profiles_folder,qc_folder);

% Subset features
plate.feaInfo.name = plate.feaInfo.name(featureset{f});
plate.feaInfo.channel_id = plate.feaInfo.channel_id(featureset{f});
plate.feaInfo.category_id = plate.feaInfo.category_id(featureset{f});
for w = 1:size(plate.profiles,1)
    plate.profiles{w,1} = plate.profiles{w,1}(featureset{f});
end
plate2 = Merge_Time_Points(plate,paras.screen.time_to_use,paras.screen.merge_mode);
[plate_bioactive,plate3,~,~,~,~] = Get_Bioactive_Compounds(plate2, paras);
is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.

% Plot PCA of all plates
Visualize_Plate_LE(plate_bioactive,plot_opt{:});
labelpcplot()
saveas(gcf,[results_folder,'\','NQO1','use_all','PCA']); % change naming by featureset{f}!

% Visualize FDA space
% Choose class model
aw_ref_extended = {'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR'};
%aw_ref = {'DMSO','DNA','HDAC','HSP90','MT','Proteasome','mTOR'};

% Load all plates
plate = Load_Batch_LE(cellplates,plateIDs, profiles_folder,qc_folder);

% Subset features
plate.feaInfo.name = plate.feaInfo.name(featureset{f});
plate.feaInfo.channel_id = plate.feaInfo.channel_id(featureset{f});
plate.feaInfo.category_id = plate.feaInfo.category_id(featureset{f});
for w = 1:size(plate.profiles,1)
    plate.profiles{w,1} = plate.profiles{w,1}(featureset{f});
end
plate2 = Merge_Time_Points(plate,paras.screen.time_to_use,paras.screen.merge_mode);
[plate_bioactive,plate3,~,~,~,~] = Get_Bioactive_Compounds(plate2, paras);
is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.

% Train FDA space
trainData = cell2mat(plate3.profiles);
trainLabels = plate3.drug_categories;
trainDataInfo = strcat(plate3.drug_names,'_',cellstr(num2str(plate3.concentrations)));
[clfy_obj, accuracy, cv_report, scan_scheme] = Classifier(trainData, trainLabels,classifier_type,'XtrainNames',trainDataInfo,'cv_fold',kfold,...
    'FDA_type',FDA_type,'gamma',gamma,'kernel_type',kernel_type,'kernel_para',kernel_para,...
    'NumNeighbors',3); % used as input to project into FDA
plate_fda = plate3; % rename for clarity
plate_fda.profiles = cellfun(@clfy_obj.KLFDA_proj_fun,plate3.profiles,'unif',false); % transform to FDA
Visualize_Plate(plate_fda,plot_opt{:});
saveas(gcf,[results_folder,'\','NQO1','use_all','FDA']);

% Clustering
% Hierarchical clustering of bioactive compounds based on profiles

% Cluster info table layout
col_order = {'Cluster','Size','Significance',...
    'plateIDs','well_names','drug_names',...
    'drug_categories','concentrations',...
    'bioactive_pVals','targets','pathways',...
    'descriptions','cas','compound',...
    'cpd_usage','nCells','clone_names','plate_types'};

% For a batch of plates
batches = cell(1,1); % legacy from other code where clustered in batches
batches{1,1} = plateIDs;
Nbatches = length(batches);
cluster_info = cell(Nbatches,1);
profiles = cell(Nbatches,1);
tree = cell(Nbatches,1);

% Perform HC
for b = 1:Nbatches
    disp(b)
    [cluster_info{b}, profiles{b}, tree{b}] = Cluster_Compounds(plate,paras);
    cluster_info{b}.compound = cluster_info{b}.drug_names;
    cluster_info{b} = cluster_info{b}(:,col_order);
end

% Visualize clusters
batch_to_plot = 1;
clusters_to_color = [];

if isempty(clusters_to_color)
    %-- get significant clusters --%
    sig_thr = 0.05;
    sig_clusters = unique(cluster_info{batch_to_plot}.Cluster(...
        cluster_info{batch_to_plot}.Significance<sig_thr));
    %-- get DMSO cluster(s) --%
    [m,n,g] = grpstats(strcmp(cluster_info{batch_to_plot}.drug_categories,'DMSO'),...
        cluster_info{batch_to_plot}.Cluster,{'mean','numel','gname'});
    dmso_cluster = str2double(g( (m>0.5) & (n>10) ));
    
    clusters_to_color = [sig_clusters(:);dmso_cluster(:)];
end

PlotPhylogeneticTree(cluster_info{batch_to_plot},tree{batch_to_plot},clusters_to_color);
saveas(gcf,[results_folder,'\','NQO1','use_all','clustering'])

% Cluster Purity
% Analysis of clusters - what percent of a drug class is in a cluster
fh = figure;
for i = 1:length(cluster_info)
    ref_cluster(i).categories = unique(cluster_info{i}.drug_categories); % find drug classes
    ref_cluster(i).cluster_ID = unique(cluster_info{i}.Cluster); % find cluster IDs
    ref_cluster(i).stat = zeros(length(ref_cluster(i).categories),length(ref_cluster(i).cluster_ID)); % make statistics table
    
    for j = 1:size(ref_cluster(i).stat,1)
        is_class = strcmpi(cluster_info{i}.drug_categories,ref_cluster(i).categories{j});
        tb = tabulate(cluster_info{i}.Cluster(is_class));
        [is_in,loc] = ismember(tb(:,1),ref_cluster(i).cluster_ID); % loc therefore sorts this matrix by cluster ID
        ref_cluster(i).stat(j,loc(is_in)) = tb(:,2);
    end
    
    % Sort the matrix by the dendrogram
    [c_id, ia] = unique(cluster_info{i}.Cluster); % ia gives the sort order of the
    [~,idx] = sort(ia);
    num_cpds = ref_cluster(i).stat(:,idx); % all the rows are sorted by idx, which is the order the columns need to go in to match the dendrogram
    c_id = c_id(idx);
    c_id_sig = ismember(c_id,[sig_clusters;dmso_cluster]);
    row_totals = sum(num_cpds,2);
    percs = num_cpds./repmat(sum(num_cpds,2),1,size(num_cpds,2));
    
    % Plot figure
    ax(i) = subplot(1,1,i,'Parent',fh);
    imagesc(percs(:,c_id_sig)),colormap(ax(i),hot),caxis(ax(i),[0,1])
    ax(i).YTick = 1:length(ref_cluster(i).categories);
    ax(i).YTickLabel = ref_cluster(i).categories;
    ax(i).TickLabelInterpreter = 'none';
end

saveas(gcf,[results_folder,'\','NQO1','use_all','cluster_purity'])

%% Checking profile subsetting

% Load all plates
plate = Load_Batch_LE(cellplates,plateIDs, profiles_folder,qc_folder);
for p = 1:length(plate.profiles)
    for w=1:length(plate.profiles{p})
        plate.profiles{p,2}(1,w) = w;
    end
end

plate_original.name = plate.feaInfo.name';
plate_original.category_id = plate.feaInfo.category_id';
plate_original.channel_id = plate.feaInfo.channel_id';
plate_original = struct2table(plate_original);

% Subsetting optiomns
[name,category_id,channel_id,NBT,use_2ch,use_all]= import_featuresselection2('W:\2015_09_HTS_LE\Code_LE\LE_Analysis\LE_MatlabScripts\feature_select_NBT.txt');
feature_table = table(name,category_id,channel_id,NBT,use_2ch,use_all);
use_2ch = use_2ch==1; % Ch1 Ch3
use_2ch = use_2ch';
use_all = use_all==1; % Ch1 Ch3
use_all = use_all';
NBT = NBT==1; % Ch1 Ch3
NBT = NBT';

plate_NBT=plate;
plate_2ch=plate;
plate_all=plate;

% Using an example profile, subset in different ways
plate_NBT.feaInfo.name = plate_NBT.feaInfo.name(NBT);
plate_NBT.feaInfo.channel_id = plate_NBT.feaInfo.channel_id(NBT);
plate_NBT.feaInfo.category_id = plate_NBT.feaInfo.category_id(NBT);
for w = 1:size(plate_NBT.profiles,1)
    plate_NBT.profiles{w,1} = plate_NBT.profiles{w,1}(NBT);
    plate_NBT.profiles{w,2} = plate_NBT.profiles{w,2}(NBT);
end

plate_2ch.feaInfo.name = plate_2ch.feaInfo.name(use_2ch);
plate_2ch.feaInfo.channel_id = plate_2ch.feaInfo.channel_id(use_2ch);
plate_2ch.feaInfo.category_id = plate_2ch.feaInfo.category_id(use_2ch);
for w = 1:size(plate_2ch.profiles,1)
    plate_2ch.profiles{w,1} = plate_2ch.profiles{w,1}(use_2ch);
    plate_2ch.profiles{w,2} = plate_2ch.profiles{w,2}(use_2ch);
end

plate_all.feaInfo.name = plate_all.feaInfo.name(use_all);
plate_all.feaInfo.channel_id = plate_all.feaInfo.channel_id(use_all);
plate_all.feaInfo.category_id = plate_all.feaInfo.category_id(use_all);
for w = 1:size(plate_all.profiles,1)
    plate_all.profiles{w,1} = plate_all.profiles{w,1}(use_all);
    plate_all.profiles{w,2} = plate_all.profiles{w,2}(use_all);
end

% Table name wrangling before joining
plate_NBT2.name = plate_NBT.feaInfo.name';
plate_NBT2.category_id = plate_NBT.feaInfo.category_id';
plate_NBT2.channel_id = plate_NBT.feaInfo.channel_id';
plate_NBT2.profile = plate_NBT.profiles{1,1}';
plate_NBT2.profileidx = plate_NBT.profiles{1,2}';
plate_NBT2 = struct2table(plate_NBT2);
plate_NBT2.Properties.VariableNames{1} = 'name_NBT2';
plate_NBT2.Properties.VariableNames{2} = 'category_id_NBT2';
plate_NBT2.Properties.VariableNames{3} = 'channel_id_NBT2';
plate_NBT2.Properties.VariableNames{4} = 'profile_NBT2';

plate_2ch2.name = plate_2ch.feaInfo.name';
plate_2ch2.category_id = plate_2ch.feaInfo.category_id';
plate_2ch2.channel_id = plate_2ch.feaInfo.channel_id';
plate_2ch2.profiles = plate_2ch.profiles{1,1}';
plate_2ch2.profileidx = plate_2ch.profiles{1,2}';
plate_2ch2 = struct2table(plate_2ch2);
plate_2ch2.Properties.VariableNames{1} = 'name_2ch';
plate_2ch2.Properties.VariableNames{2} = 'category_id_2ch';
plate_2ch2.Properties.VariableNames{3} = 'channel_id_2ch';
plate_2ch2.Properties.VariableNames{4} = 'profiles_2ch';

plate_all2.name = plate_all.feaInfo.name';
plate_all2.category_id = plate_all.feaInfo.category_id';
plate_all2.channel_id = plate_all.feaInfo.channel_id';
plate_all2.profile = plate_all.profiles{1,1}';
plate_all2.profileidx = plate_all.profiles{1,2}';
plate_all2 = struct2table(plate_all2);
plate_all2.Properties.VariableNames{1} = 'name_all2';
plate_all2.Properties.VariableNames{2} = 'category_id_All2';
plate_all2.Properties.VariableNames{3} = 'channel_id_all2';
plate_all2.Properties.VariableNames{4} = 'profile_all2';

clear plate_NBT plate_2ch plate_all

% Join the tables to show indexing finds the right number for a feature
bigtb = outerjoin(plate_2ch2,plate_all2,'Keys', 'profileidx');
bigtb.Properties.VariableNames{10} = 'profileidx';
bigtb = outerjoin(plate_NBT2,bigtb,'Keys', 'profileidx');

writetable(bigtb,[results_folder,'\profilesubsetting.xlsx'])