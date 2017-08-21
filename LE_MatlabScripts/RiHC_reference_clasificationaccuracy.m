% Feature Comparison RiHC line (pSeg only)
% Louise Heinrich in Altschuler & Wu Labs
% August 15th 2017

%% Set Paths
addpath(genpath('W:\Lab_Analysis\common\plate_annotation')) %plate annotation from labtoolbox
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\CD_Tag_Drug_Screen')) % Charles code for NBT paper containing useful functions
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_analysis')) % Louise code afer cluster operations
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_Preprocessing')) % Louise code before cluster operations
addpath(genpath('W:\Lab_Analysis\common\image_browser')) % Image_Browser to callback images

%% Environment & Variables
code_folder = 'W:\2015_09_HTS_LE\Code_LE\';
profiles_folder = 'W:\2015_09_HTS_LE\data\profiles';% place to save profiles to and load from
qc_folder = 'W:\2015_09_HTS_LE\QC\bigscreen';% store any QC analysis Python, Matlab, or manual
results_folder = 'W:\2015_09_HTS_LE\results\feature_sets\RIHC'; % Louise results for 4lines analysis
plate_DB_path ='W:\2015_09_HTS_LE\project_database\'; % directory containing plate database file
plate_DB = 'Plate_database_latest2.xlsx'; % file describing relationship between compound and experiment plates

% Read in plate table for plates of interest
Exp_DB = readtable(fullfile(plate_DB_path,plate_DB),'Sheet',1); % full table of all plates
cellplates = Exp_DB(strcmpi(Exp_DB.Experiment_,'LE_8'),:); % Select RiHC x Reference set plates 2017018405,403

% Get Analysis Parameters
% Get plotting options
getplotoptions() % generate standard plot options
plot_opt = plot_opt_reference; % specific plot option

%bioactive_paras
pVal_thr = 10^-2; %10^-6
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

% Subset Setup
% Import indexing of features to select
[name,category_id,channel_id,NBT,use_2ch,use_all]= import_featuresselection2('W:\2015_09_HTS_LE\Code_LE\LE_Analysis\LE_MatlabScripts\feature_select_NBT.txt'); % manually indexed feature list based on NBT feature set from supplementary tables

% create logicals for subsetting features
NBT = NBT==1; % NBT set 234 features focusing on YFP channel
NBT = NBT';
use_2ch = use_2ch==1; % Ch1 Ch3 443 features not using any YFP features
use_2ch = use_2ch';
use_all = use_all==1; % Ch1 Ch2 Ch3 659 all features all channels
use_all = use_all';

clear name category_id channel_id w
clc

%% Generate Full Profiles
cellline = 'A549_RiHC';
LE_Gen_Profiles(cellplates,cellline,profiles_folder);
%% Choose Feature Set

featureset ={NBT,use_2ch,use_all}; % use_2ch use_all NBT
featurenames = {'NBT','use_2ch','use_all'};
for f = 1
    % Calc. Bioactivity
    % Pick 6 class or 10 class model for classification accuracy
    aw_ref = {'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR'};
    %aw_ref = {'DMSO','DNA','HDAC','HSP90','MT','Proteasome','mTOR'};
    
    % Import profiles
    plateIDs = cellplates.expt_plate;
    plate = Load_Batch_LE(cellplates,plateIDs, profiles_folder,qc_folder);
    
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
    %%% outputs : [plate_bioactive,plate,plate_ref_bioactive,plate_ref,plate_cpd_bioactive,plate_cpd]
    is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO')); %inactive cpds
    plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
    bioactive = Plate2Table(plate3);% make table
    clear plate2 is_inactive
    writetable(bioactive,[results_folder, '\', cellline,featurenames{f},'_bioactive'])
    
    % Bioactivity Percents
    bioactive.cpd_usage = categorical(bioactive.cpd_usage);
    bioactive.Properties.VariableNames{7} = 'bioactive_pVals';
    reference = bioactive(bioactive.cpd_usage =='reference_cpd',:);
    percents.ref = sum(reference.bioactive_pVals <= pVal_thr)/height(reference)*100;
    %percents_summary = struct2table(percents.ref);
    %percents_summary.Properties.RowNames{1}='Reference';
    display(percents.ref)
    %writetable(percents.ref,[results_folder,'\',featurenames{f},'percents_bwDMSO.xlsx'])
    
    % PCA FDA Visualisation
    % Plot PCA of all plates
    Visualize_Plate_LE(plate3,plot_opt_ALL{:});
    labelpcplot()
    saveas(gcf,[results_folder,'\',cellline,featurenames{f},'PCA','_2K']); % change naming by featureset!
    
    % Calc FDA space
    % Get reference compounds
    plate3_r = Select_Rows(plate3,ismember(plate3.drug_categories,aw_ref));
    %plate3_r = Select_Rows(plate3_r,ismember(plate3_r.plate_types,{'reference_plate'}));
    % Train FDA space
    trainData = cell2mat(plate3_r.profiles);
    trainLabels = plate3_r.drug_categories;
    trainDataInfo = strcat(plate3_r.drug_names,'_',cellstr(num2str(plate3_r.concentrations)));
    [clfy_obj, accuracy, cv_report, scan_scheme] = Classifier(trainData, trainLabels,classifier_type,'XtrainNames',trainDataInfo,'cv_fold',kfold,...
        'FDA_type',FDA_type,'gamma',gamma,'kernel_type',kernel_type,'kernel_para',kernel_para,...
        'NumNeighbors',3); % used as input to project into FDA
    plate_fda = plate3; % rename for clarity
    plate_fda.profiles = cellfun(@clfy_obj.KLFDA_proj_fun,plate3.profiles,'unif',false); % transform to FDA
    
    % Project FDA
    Visualize_Plate(plate_fda,plot_opt_ALL{:});
    saveas(gcf,[results_folder,'\',cellline,featurenames{f},'FDA']);
    
    
    % Clustering & Cluster Analysis
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
        [cluster_info{b}, profiles{b}, tree{b}] = Cluster_Compounds(plate3,paras);
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
    save([results_folder,'XRCC5',featurenames{f},'clustering','_2K'],'cluster_info')
    saveas(gcf,[results_folder,'\',cellline,featurenames{f},'clustering'])
    
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
        ax(i) = subplot(2,1,1,'Parent',fh);
        dendrogram(tree{1,1}(1,1))
%       imagesc(percs(:,c_id_sig)),colormap(ax(i),hot),caxis(ax(i),[0,1])
        subplot(2,1,2,'Parent',fh)        
        imagesc(percs)

        ax(i).YTick = 1:length(ref_cluster.categories);
        ax(i).YTickLabel = ref_cluster.categories;
        ax(i).TickLabelInterpreter = 'none';
        
    end
    
    saveas(gcf,[results_folder,'\',cellline,featurenames{f},'cluster_purity'])
end

