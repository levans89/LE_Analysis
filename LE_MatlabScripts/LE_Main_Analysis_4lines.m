% Louise Evans March 2017 Purpose: Perform analogous analysis to Charles'
% NBT paper for 4lines experiment x Selleck 2K set 'pORACL' determination.

%% Set Paths
addpath(genpath('W:\Lab_Analysis\common\plate_annotation')) %plate annotation from labtoolbox
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\CD_Tag_Drug_Screen')) % Charles code for NBT paper containing useful functions
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_analysis')) % Louise code afer cluster operations
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_Preprocessing')) % Louise code before cluster operations

%% Define Variables
code_folder = 'W:\2015_09_HTS_LE\Code_LE';
profiles_folder = 'W:\2015_09_HTS_LE\data\profiles';% place to save profiles to and load from
profiles_folder_RB = '..\data\profiles_4lines_XRCC5_RB'; % for comparison of ill corr methods on C.A. %
profiles_folder_NBT = 'W:\2015_09_HTS_LE\data\profiles_NBT';% subset of features used in NBT paper
qc_folder = 'W:\2015_09_HTS_LE\QC\4lines\stringent';% store any QC analysis Python, Matlab, or manual
results_folder = 'W:\2015_09_HTS_LE\results\matlab_results\4lines'; % Louise results for 4lines analysis
plate_DB_path ='W:\2015_09_HTS_LE\project_database\'; % directory containing plate database file
plate_DB = 'Plate_database_latest2.xlsx'; % file describing relationship between compound and experiment plates

%% Load Constants
Exp_DB = readtable(fullfile(plate_DB_path,plate_DB),'Sheet',1);

year = 2017;
theseyearplates = ismember(Exp_DB.year,year);
plates = Exp_DB(theseyearplates,:);

expt = 'LE_4';
thesebatchplates = ismember(plates.Experiment_,expt);
plates = plates(thesebatchplates,:);

batch='4L2K_7';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);

batch='4L2K_8';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);

referenceplates = plates((ismember(plates.plate_type,'reference_plate')),:);

pORACL = {'A549_XRCC5';'A549_NQO1';'A549_SET';'A549_S100A11'};
clear theseyearplates thesebatchplates thesereferenceplates

%% Generate Profiles
% Generate profile for full feature set
LE_Gen_Profiles(plates,'A549_XRCC5',profiles_folder)
LE_Gen_Profiles(plates,'A549_SET',profiles_folder)
LE_Gen_Profiles(plates,'A549_S100A11',profiles_folder)
LE_Gen_Profiles(plates,'A549_NQO1',profiles_folder)

% Subset Profiles to NBT features
for p = 1:height(plates)
    plateID = plates.expt_plate{p};
    load(fullfile(profiles_folder,strcat('plate_',plateID,'.mat')));
    [~,~,~,use,~,~] = import_featuresselection('W:\2015_09_HTS_LE\Code_LE\LE_Analysis\LE_MatlabScripts\feature_select_NBT.txt');
    for w = 1:size(plate.profiles,1)
        plate.profiles{w,1} = plate.profiles{w,1}(use);
    end
    plate.feaInfo.name = plate.feaInfo.name(use);
    plate.feaInfo.channel_id = plate.feaInfo.channel_id(use);
    plate.feaInfo.category_id = plate.feaInfo.category_id(use);
    save(fullfile(profiles_folder_NBT,strcat('plate_',plateID,'.mat')),'plate')
end
clear use p plateID plate

%% Get Analysis Parameters
% Bioactive Paras
pVal_thr         = 10^-6; %10^-6
max_nCtrl        = 2500; %2500
var_pct_to_keep  = 0.95; %0.95
k_fold           = 10; %10

% Screen Paras
time_to_use      = 1; %1 | 2 | 'all'
merge_mode       = 'concatenate';% 'concatenate' | 'pool'
classifier_type  = 'KLFDA_knn'; % 'KLFDA_knn' | 'KLFDA_nearest_centroid'
FDA_type         = 0; % 0: FDA | -1: no transformation | -2: RCA
gamma            = 0:0.05:0.5; % regularization strength for within-class variance
kernel_type      = 'linear';
kernel_para      = 0;
NumNeighbors     = 3;
kfold            = 10;
conf_thr         = 0.1;
do_view          = true;
KNN_to_list      = 10;

% Cluster Paras
linkage_method   = 'average';
ref_flag         = {'reference_cpd','negative_ctrl'};
nDis_max         = 50000;
D_within_prob_thr= 0.5;
min_dmso_pct     = 0.5;
min_dmso_num     = 10;

% Set Paras
classifier_paras = {'NumNeighbors',NumNeighbors};
HCSparas         = Set_Paras('pVal_thr',pVal_thr,'max_nCtrl',max_nCtrl,'var_pct_to_keep',var_pct_to_keep,...
                  'k_fold',k_fold,'time_to_use',time_to_use,'merge_mode',merge_mode,'classifier_type',classifier_type,...
                  'FDA_type',FDA_type,'gamma',gamma,'kernel_type',kernel_type,'kernel_para',kernel_para,...
                  'classifier_paras',classifier_paras,'kfold',kfold,'conf_thr',conf_thr,'do_view',do_view,...
                  'linkage_method',linkage_method,'ref_flag',ref_flag,'nDis_max',nDis_max,'D_within_prob_thr',D_within_prob_thr,...
                  'min_dmso_pct_to_be_DMSO_cluster',min_dmso_pct,'min_dmso_number_to_be_DMSO_cluster',min_dmso_num,'KNN_to_list',KNN_to_list);% set parameters

% now plates can be categorical
plates.CellLine = categorical(plates.CellLine);
referenceplates.CellLine = categorical(referenceplates.CellLine);
plates.CellLine = categorical(plates.CellLine);            
getplotoptions()              
%% Visualize Plates & Specify Bioactives
% Generate standard plot options
plot_opt = plot_opt_ALL;

% Visualize plates
for c = 1%:4
    cellplates = plates(plates.CellLine == pORACL{c},:);
    plateIDs = cellplates.expt_plate;
    plate = Load_Batch_LE(cellplates,plateIDs, profiles_folder,qc_folder);
    
    anno_updater;
    
    checknumplatewells()
    plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);
 
    [plate_bioactive,...
     plate3,...
     plate_ref_bioactive,...
     plate_ref,...
     plate_cpd_bioactive,...
     plate_cpd] = Get_Bioactive_Compounds(plate2, HCSparas);
 
    is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
    plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
    bioactive = Plate2Table(plate3);
    %writetable(bioactive,fullfile(results_folder,strcat(pORACL{c},'_bioactivepvals_DMSObw_NBT.xlsx')))
    
    Visualize_Plate_LE(plate3,plot_opt{:});
    labelpcplot()
    title(pORACL{c});
    %saveas(gcf,fullfile(results_folder,pORACL{c}))
end

%%  Extract bioactivity files for RB XRCC5
plates = plates(plates.CellLine== strcat('A549_',pORACL{1}),:); % select plateset 
plateIDs = plates.expt_plate;
plate = Load_Batch_LE(plates,plateIDs, profiles_folder_RB,qc_folder);%note different profiles folder accessed here from S method
checknumplatewells()
plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
[~,plate3] = Get_Bioactive_Compounds(plate2, HCSparas);
is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
bioactive = Plate2Table(plate3);
writetable(bioactive,fullfile(results_folder,strcat(pORACL,'_bioactivepvals.xlsx')))
%writetable(bioactive, fullfile(results_folder,strcat(pORACL,'_bioactivepvals_nonbioactivenotlabelled_RB.xlsx')))

%% Comparison pORACL Bioactivity
NQO1_bioactive=bioactives_importfile(fullfile(results_folder,'A549_NQO1_bioactivepvals_DMSObw_NBT.xlsx'),'Sheet1');
SET_bioactive=bioactives_importfile(fullfile(results_folder,'A549_SET_bioactivepvals_DMSObw_NBT.xlsx'),'Sheet1');
S100A11_bioactive=bioactives_importfile(fullfile(results_folder,'A549_S100A11_bioactivepvals_DMSObw_NBT.xlsx'),'Sheet1');
XRCC5_bioactive = bioactives_importfile(fullfile(results_folder,'A549_XRCC5_bioactivepvals_DMSObw_NBT.xlsx'),'Sheet1');

NQO1_bioactive(3639:end,:) = [];
XRCC5_bioactive(3640:end,:) = [];
S100A11_bioactive(3646:end,:) = [];
SET_bioactive(3638:end,:) = [];

pvalue_plotter % plot P values by cell line, and by batch

%% Classification Accuracy
% Pick out AW reference
%aw_ref = {'DMSO','DNA','HDAC','HSP90','MT','Proteasome','mTOR'}; % 6 classes
aw_ref = {'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR'}; %10 classes

ref_plate = {'reference_plate'};

% Calculate classification accuracy for reference plates
for c = 1%:4
plateaccuracy = cell(6,2);    
refplates = referenceplates(referenceplates.CellLine== pORACL{c},:);
for p=1%:height(refplates)
    plate_IDs = refplates.expt_plate;
    plateIDs = plate_IDs(p);
    plate = Load_Batch_LE(refplates,plateIDs, profiles_folder,qc_folder);
    plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);
    plate2 = Get_Bioactive_Compounds(plate2, HCSparas);
    % Get reference compounds
    plate2 = Select_Rows(plate2,ismember(plate2.drug_categories,aw_ref));
    plate2 = Select_Rows(plate2,ismember(plate2.plate_types,{'reference_plate'}));
    % Prepare training data
    trainData = cell2mat(plate2.profiles);
    trainLabels = plate2.drug_categories;
    trainDataInfo = strcat(plate2.drug_names,'_',cellstr(num2str(plate2.concentrations)));
    % Cross-validation
    [clfy_obj, accuracy, cv_report, scan_scheme] = Classifier(trainData, trainLabels, ...
                                                   classifier_type,'XtrainNames',trainDataInfo,'cv_fold',kfold,...
                                                   'FDA_type',FDA_type,'gamma',gamma,...
                                                   'kernel_type',kernel_type,'kernel_para',kernel_para,...
                                                   'NumNeighbors',NumNeighbors);
    % record results
    plateaccuracy{p,2} = accuracy;
    plateaccuracy{p,1}= plateIDs{1,1};
end
plateaccuracy = cell2table(plateaccuracy);
%writetable(plateaccuracy,fullfile(results_folder,strcat(pORACL{c},'_Accuracies_6.xlsx')))
end

%% Save classifier results
classifier.XRCC5 = XRCC5_classifier;
classifier.SET = SET_classifier;
classifier.S100A11 = S100A11_classifier;
classifier.NQO1 = NQO1_classifier;
classifier.XRCC5_RB = XRCC5_classifier_RB;
accuracies = [NQO1_classifier.accuracy,XRCC5_classifier.accuracy,SET_classifier.accuracy,S100A11_classifier.accuracy, XRCC5_classifier_RB.accuracy]';
accuraciesh = {'NQO1','XRCC5','SET','S100A11', 'XRCC5_RB'}';
accuracies = table(accuracies);
CellLine = table(accuraciesh);
Accuracies = [CellLine accuracies];
Accuracies.Properties.VariableNames{1} = 'CellLine';
Accuracies.Properties.VariableNames{2} = 'Accuracy';
Accuracies = sortrows(Accuracies,'Accuracy','descend');
%writetable(Accuracies,fullfile(results_folder,'Accuracies.xlsx')) % HOW TO GET SD??? CURRENTLY NEED TO MANUALLY SAVE - COME BACK TO THIS

%% Visualize the FDA space 
%drug_categories = aw_ref;
drug_category_strimmer

plates = plates(plates.CellLine=='A549_XRCC5',:);
plateIDs = cellplates.expt_plate;
plate = Load_Batch_LE(cellplates,plateIDs, profiles_folder,qc_folder); % note different profiles folder

anno_updater;

plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);
plate2 = Get_Bioactive_Compounds(plate2, HCSparas);
plate2 = Select_Rows(plate2,ismember(plate2.drug_categories,drug_categories));

% Prepare training data
trainData = cell2mat(plate2.profiles);
trainLabels = plate2.drug_categories;
trainDataInfo = strcat(plate2.drug_names,'_',cellstr(num2str(plate2.concentrations)));

[clfy_obj, accuracy, cv_report, scan_scheme] = Classifier(trainData, trainLabels, ...
                                               classifier_type,'XtrainNames',trainDataInfo,'cv_fold',kfold,...
                                               'FDA_type',FDA_type,'gamma',gamma,...
                                               'kernel_type',kernel_type,'kernel_para',kernel_para,...
                                               'NumNeighbors',NumNeighbors);

plate_fda = plate2;
plate_fda.profiles = cellfun(@clfy_obj.KLFDA_proj_fun,plate2.profiles,'unif',false);
getplotoptions();
plot_opt_ALL{1,2}= categories(categorical(plate.drug_categories))';
Visualize_Plate(plate_fda,plot_opt_ALL{:});
%saveas(gcf,fullfile(results_folder,pORACL));

%% Predictions
% Load plate
plates = plates(plates.CellLine=='A549_XRCC5',:);
plateIDs = plates.expt_plate;
plate = Load_Batch_LE(plates,plateIDs, profiles_folder_NBT,qc_folder); % note different profiles folder

% Predict Output
predict_output = Screen_Plates(plate, HCSparas);
save(fullfile(results_folder,'predict_output_XRCC5.mat'),'predict_output')

%% Clustering

col_order = {'Cluster','Size','Significance',...
             'plateIDs','well_names','drug_names',...
             'drug_categories','concentrations',...
             'bioactive_pVals','targets','pathways',...
             'descriptions','cas','compound',...
             'cpd_usage','nCells','clone_names','plate_types'};

sort_order = {'Size','Cluster','pathways','targets',...
              'drug_names','concentrations'};

%batches = cell(4,1);
for c = 1%:4
plates = plates(plates.CellLine==pORACL{c},:);
plateIDs = plates.expt_plate;
batches{c,1} = plateIDs;
end

Nbatches = length(batches);
cluster_info = cell(Nbatches,1);
profiles = cell(Nbatches,1);
tree = cell(Nbatches,1);
for b = 1:Nbatches
    disp(b)
    plate = Load_Batch_LE(plates,plate_IDs, profiles_folder_NBT,qc_folder);
    [cluster_info{b}, profiles{b}, tree{b}] = Cluster_Compounds(plate_fda,HCSparas);
    cluster_info{b}.compound = cluster_info{b}.drug_names;
    cluster_info{b} = cluster_info{b}(:,col_order);
%     [cluster_info{b},idx] = sortrows(cluster_info{b},sort_order); % #
%     profiles{b} = profiles{b}(idx,:); % #
%     tree{b} = Update_NodeIDs(tree{b},idx); % #
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
%saveas(gcf,strcat(resultsfolder,'clustering'))
%% Analysis of clusters
fh = figure;
ref_cluster = struct();
for i = 1:length(cluster_info)
    ref_cluster(i).categories = unique(cluster_info{i}.drug_categories);
    ref_cluster(i).cluster_ID = unique(cluster_info{i}.Cluster);
    ref_cluster(i).stat = zeros(length(ref_cluster(i).categories),length(ref_cluster(i).cluster_ID));
    
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
    
    % Plot
    ax(i) = subplot(1,1,i,'Parent',fh);
    imagesc(percs(:,c_id_sig)),colormap(ax(i),hot),caxis(ax(i),[0,1])
    ax(i).YTick = 1:length(ref_cluster(i).categories);
    ax(i).YTickLabel = ref_cluster(i).categories;
    ax(i).TickLabelInterpreter = 'none';
    ax(i).XTick = [];
end

%% Louise Percents
[targets, pathways, compound_IDs]=getselleckpathways();
pathways.Pathway = categorical(pathways.Pathway);
pathwaysize = tabulate(pathways.Pathway);
query_cluster = struct();
for i = 1:height(cluster_info{1,1})
    query_cluster(i).categories = pathwaysize(:,1);
    query_cluster(i).cluster_ID = unique(cluster_info{i}.Cluster);
    query_cluster(i).stat = zeros(length(query_cluster(i).categories),length(query_cluster(i).cluster_ID));
    for j = 1:size(query_cluster(i).stat,1)
        is_class = strcmpi(cluster_info{i}.drug_categories,query_cluster(i).categories{j});
        tb = tabulate(cluster_info{i}.Cluster(is_class));

        [is_in,loc] = ismember(tb(:,1),query_cluster(i).cluster_ID);
        query_cluster(i).stat(j,loc(is_in)) = tb(:,2);
    end
end

ref_cluster.categories2 = cell2table(ref_cluster.categories);
ref_cluster.counts2 = cell2table(tabulate(cluster_info{1}.drug_categories));
ref_cluster.drugcounts = join(ref_cluster.categories2,ref_cluster.counts2);
ref_cluster.clustercounts=tabulate(cluster_info{1}.Cluster);

drugcounts = table2array(ref_cluster.drugcounts(:,2));
drugcounts = repmat(drugcounts,1,117);
x = ref_cluster.stat./drugcounts;

clustercounts = ref_cluster.clustercounts(:,2)';
clustercounts = repmat(clustercounts,32,1);
y = ref_cluster.stat./clustercounts;

x_ax = x(:);
y_ax = y(:);

i =1;
fh = figure;
ax(i) = subplot(1,1,i,'Parent',fh);

imagesc(x)
colormap(ax(i),'hot')
colorbar(ax(i))
ax(i).YTick = 1:length(ref_cluster(i).categories);
ax(i).YTickLabel = ref_cluster(i).categories;
ax(i).TickLabelInterpreter = 'none';

xlabel('Cluster #')
ylabel('Compound Category (Pathway)')
print('-f1','XRCC5clusterper','-dpng','-r500')

%% Print Physical Plate Cell Number
% Generate heatmaps of cell number for all plates to check for pattern
cd (qc_folder) % move to the location to do all the saving in to improve speed and make checking for files quicker
for c = 1:4
    plates = plates(plates.CellLine==pORACL{c},:);
    plateIDs = plates.expt_plate;
    for p = 1:length(plateIDs)
        filename =char(strcat(plateIDs(p),'.bmp'));
        display(filename)
        plate = Load_Batch_LE(plates,plateIDs(p), profiles_folder_NBT,qc_folder); % load plate
        Show_Cell_Number(plate); % generate heatmap of cell number
        saveas(gcf,plateIDs{1,1},'bmp'); % save
    end
end
cd (code_folder) % move back

%% Print Physical Plate P Values
cd (qc_folder)
for c = 1:4
    plates = plates(plates.CellLine==pORACL{c},:);
    plateIDs = plates.expt_plate;
    for p = 1:length(plateIDs)
        plate = Load_Batch_LE(plates,plateIDs(p), profiles_folder_NBT,qc_folder);
        plate2 = Merge_Time_Points(plate,time_to_use,merge_method);
        [plate2,~] = Get_Bioactive_Compounds(plate2, HCSparas);
        plate2.logbioactive_pVals = log(plate2.bioactive_pVals);
        plate2.logbioactive_pVals(isnan(plate2.logbioactive_pVals))=1;%A(isnan(A)) = 0
        Show_PVALUE(plate2)
        title = ['../results/matlab_results/4lines/',plateIDs{p},' Log P Value Plate'];
        print('-f1',title,'-dpng','-r500')
        close all
    end
end
cd (code_folder) % move back
%% Concatenate Profiles
% theseplates = plates(plates.CellLine=='A549_XRCC5',:);
% plateIDs = theseplates.expt_plate;
% XRCC5 = Load_Batch_LE(plates,plateIDs,profiles_folder_NBT,qc_folder);
% 
% theseplates = plates(plates.CellLine=='A549_NQO1',:);
% plateIDs = theseplates.expt_plate;
% NQO1 = Load_Batch_LE(plates,plateIDs,profiles_folder_NBT,qc_folder);
% 
% theseplates = plates(plates.CellLine=='A549_SET',:);
% plateIDs = theseplates.expt_plate;
% SET = Load_Batch_LE(plates,plateIDs,profiles_folder_NBT,qc_folder);
% 
% theseplates = plates(plates.CellLine=='A549_S100A11',:);
% plateIDs = theseplates.expt_plate;
% S100A11 = Load_Batch_LE(plates,plateIDs,profiles_folder_NBT,qc_folder);

all4          = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_all4.mat');
all3.XSetS100 = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_XRCC5_SET_S100A11.mat');
all3.XSetN    = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_XRCC5_SET_NQO1.mat');
all3.SSN      = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_SET_NQO1_S100A11.mat');
all3. XNS100  = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_XRCC5_NQO1_S100A11.mat');


all2.XSet     = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_XRCC5_SET.mat');
all2. XS100   = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_XRCC5_S100A11.mat');
all2.XN       = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_XRCC5_NQO1.mat');
all2.SS       = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_SET_S100A11.mat');
all2.SetN     = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_SET_NQO1.mat');
all2.S100N    = load('W:\2015_09_HTS_LE\data\profiles_concatenations\catProfile_NQO1_S100A11.mat');

[plate_bioactive,...
    plate3,...
    plate_ref_bioactive,...
    plate_ref,...
    plate_cpd_bioactive,...
    plate_cpd] = Get_Bioactive_Compounds(all4.catProfile., HCSparas);

is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
bioactive = Plate2Table(plate3);
writetable(bioactive,fullfile(results_folder,strcat(pORACL{c},'_bioactivepvals_DMSObw_NBT.xlsx')))

Visualize_Plate_LE(all4.catProfile,plot_opt_query{:});
labelpcplot()
title(pORACL{c});

% Clustering

col_order = {'Cluster','Size','Significance',...
             'plateIDs','well_names','drug_names',...
             'drug_categories','concentrations',...
             'bioactive_pVals','targets','pathways',...
             'descriptions','cas','compound',...
             'cpd_usage','nCells','clone_names','plate_types'};

sort_order = {'Size','Cluster','pathways','targets',...
              'drug_names','concentrations'};

%batches = cell(4,1);
for c = 1%:4
plates = plates(plates.CellLine==pORACL{c},:);
plateIDs = plates.expt_plate;
batches{c,1} = plateIDs;
end

Nbatches = length(batches);
cluster_info = cell(Nbatches,1);
profiles = cell(Nbatches,1);
tree = cell(Nbatches,1);
for b = 1:Nbatches
    disp(b)
    %plate = Load_Batch_LE(plates,plate_IDs, profiles_folder_NBT,qc_folder);
    [cluster_info{b}, profiles{b}, tree{b}] = Cluster_Compounds(plate_fda,HCSparas);
    cluster_info{b}.compound = cluster_info{b}.drug_names;
    cluster_info{b} = cluster_info{b}(:,col_order);
%     [cluster_info{b},idx] = sortrows(cluster_info{b},sort_order); % #
%     profiles{b} = profiles{b}(idx,:); % #
%     tree{b} = Update_NodeIDs(tree{b},idx); % #
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
%saveas(gcf,strcat(resultsfolder,'clustering'))

% Analysis of clusters
fh = figure;
ref_cluster = struct();
for i = 1:length(cluster_info)
    ref_cluster(i).categories = unique(cluster_info{i}.drug_categories);
    ref_cluster(i).cluster_ID = unique(cluster_info{i}.Cluster);
    ref_cluster(i).stat = zeros(length(ref_cluster(i).categories),length(ref_cluster(i).cluster_ID));
    
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
    
    % Plot
    ax(i) = subplot(1,1,i,'Parent',fh);
    imagesc(percs(:,c_id_sig)),colormap(ax(i),hot),caxis(ax(i),[0,1])
    ax(i).YTick = 1:length(ref_cluster(i).categories);
    ax(i).YTickLabel = ref_cluster(i).categories;
    ax(i).TickLabelInterpreter = 'none';
    ax(i).XTick = [];
end

% Louise Percents
[targets, pathways, compound_IDs]=getselleckpathways();
pathways.Pathway = categorical(pathways.Pathway);
pathwaysize = tabulate(pathways.Pathway);
query_cluster = struct();
for i = 1:height(cluster_info{1,1})
    query_cluster(i).categories = pathwaysize(:,1);
    query_cluster(i).cluster_ID = unique(cluster_info{i}.Cluster);
    query_cluster(i).stat = zeros(length(query_cluster(i).categories),length(query_cluster(i).cluster_ID));
    for j = 1:size(query_cluster(i).stat,1)
        is_class = strcmpi(cluster_info{i}.drug_categories,query_cluster(i).categories{j});
        tb = tabulate(cluster_info{i}.Cluster(is_class));

        [is_in,loc] = ismember(tb(:,1),query_cluster(i).cluster_ID);
        query_cluster(i).stat(j,loc(is_in)) = tb(:,2);
    end
end

ref_cluster.categories2 = cell2table(ref_cluster.categories);
ref_cluster.counts2 = cell2table(tabulate(cluster_info{1}.drug_categories));
ref_cluster.drugcounts = join(ref_cluster.categories2,ref_cluster.counts2);
ref_cluster.clustercounts=tabulate(cluster_info{1}.Cluster);

drugcounts = table2array(ref_cluster.drugcounts(:,2));
drugcounts = repmat(drugcounts,1,117);
x = ref_cluster.stat./drugcounts;

clustercounts = ref_cluster.clustercounts(:,2)';
clustercounts = repmat(clustercounts,32,1);
y = ref_cluster.stat./clustercounts;

x_ax = x(:);
y_ax = y(:);

i =1;
fh = figure;
ax(i) = subplot(1,1,i,'Parent',fh);

imagesc(x)
colormap(ax(i),'hot')
colorbar(ax(i))
ax(i).YTick = 1:length(ref_cluster(i).categories);
ax(i).YTickLabel = ref_cluster(i).categories;
ax(i).TickLabelInterpreter = 'none';

xlabel('Cluster #')
ylabel('Compound Category (Pathway)')
print('-f1','XRCC5clusterper','-dpng','-r500')



























