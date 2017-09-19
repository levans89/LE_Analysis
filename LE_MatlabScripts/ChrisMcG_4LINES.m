% Louise Evans March 2017 Purpose: Perform analogous analysis to Charles'
% NBT paper for 4lines experiment x Selleck 2K set 'pORACL' determination.

%% Set Paths
addpath(genpath('W:\Lab_Analysis\common\plate_annotation')) %plate annotation from labtoolbox
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\CD_Tag_Drug_Screen')) % Charles code for NBT paper containing useful functions
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_analysis')) % Louise code afer cluster operations
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_Preprocessing')) % Louise code before cluster operations

%% Define Variables
code_folder = 'W:\2015_09_HTS_LE\Code_LE';
profiles_folder_NBT = 'W:\2015_09_HTS_LE\data\profiles_NBT';% subset of features used in NBT paper
qc_folder = 'W:\2015_09_HTS_LE\QC\4lines\stringent';% store any QC analysis Python, Matlab, or manual
results_folder = 'W:\2015_09_HTS_LE\results\TOPHER_FINAL';
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

pORACL = {'A549_XRCC5';'A549_NQO1';'A549_SET';'A549_S100A11'};
clear theseyearplates thesebatchplates thesereferenceplates plate_DB plate_DB_path batch expt Exp_DB year

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
clear pVal_thr max_nCtrl var_pct_to_keep k_fold time_to_use merge_mode classifier_type FDA_type gamma kernel_type kernel_para NumNeighbors kfold...            = 10;
conf_thr do_view KNN_to_list linkage_method ref_flag nDis_max D_within_prob_thr min_dmso_pct min_dmso_num classifier_paras

%% Visualize Plates & Specify Bioactives
c=1; % Set c to XRCC5
plates = plates(plates.CellLine == pORACL{1},:);
plates(4,:) = [];
plates(4,:) = [];
plates(6:7,:) = [];
plates(end,:) = [];

plateIDs = plates.expt_plate;
plate = Load_Batch_LE(plates,plateIDs, profiles_folder_NBT,qc_folder);
anno_updater;
plate = Merge_Time_Points(plate,HCSparas.screen.time_to_use,HCSparas.screen.merge_mode);

[plate_bioactive,...
 plate,...
 plate_ref_bioactive,...
 plate_cpd_bioactive] = Get_Bioactive_Compounds(plate, HCSparas);

% Record these variables
plateT = Plate2Table(plate)
bioactive = Plate2Table(plate_bioactive);
bioactive_query = Plate2Table(plate_cpd_bioactive);
bioactive_ref = Plate2Table(plate_ref_bioactive);
writetable(plate,[results_folder,'\', pORACL{c},'plate.xlsx']);
writetable(bioactive,[results_folder,'\', pORACL{c},'_bioactive.xlsx']);
writetable(bioactive_query,[results_folder,'\', pORACL{c},'_bioactive_query.xlsx']);
writetable(bioactive_ref,[results_folder,'\', pORACL{c},'_bioactive_ref.xlsx']);

% Do PCA
plot_opt = plot_opt_ALL;
is_inactive = plate.bioactive_pVals>=HCSparas.bioactive.pVal_thr & (~strcmp(plate.drug_categories,'DMSO'));
plate.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
Visualize_Plate_LE(plate,plot_opt{:},'dim',2);
labelpcplot()
title pORACL{c};
saveas(gcf,[results_folder,'\',pORACL{c},'PCA.fig'])

%% Clustering

col_order = {'Cluster','Size','Significance',...
             'plateIDs','well_names','drug_names',...
             'drug_categories','concentrations',...
             'bioactive_pVals','targets','pathways',...
             'descriptions','cas','compound',...
             'cpd_usage','nCells','clone_names','plate_types'};

sort_order = {'Size','Cluster','pathways','targets',...
              'drug_names','concentrations'};

batches = cell(1,1);
batches{c,1} = plateIDs;


clusterplate = {plate}%, {plate_fda};
clusterplaten = {'plate'}%, {'plate_fda'};
for p = 1:2
Nbatches = length(batches);
cluster_info = cell(Nbatches,1);
profiles = cell(Nbatches,1);
tree = cell(Nbatches,1);
for b = 1%:Nbatches
    [cluster_info{b}, profiles{b}, tree{b}] = Cluster_Compounds(clusterplate{p},HCSparas);
    cluster_info{b}.compound = cluster_info{b}.drug_names;
    cluster_info{b} = cluster_info{b}(:,col_order);
    writetable(cluster_info{1,1},[results_folder,'\',pORACL{c},clusterplaten{p},'clustering'])
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
saveas(gcf,[results_folder,'\',cellline{c},clusterplaten{p},'clustering'])
end

%% Output
array2table(profiles{1,1})
cluster_table = horzcat(cluster_info{1,1},array2table(profiles{1,1}))
bioactive = Plate2Table(plate_bioactive);
bioactive_query = Plate2Table(plate_cpd_bioactive);
bioactive_ref = Plate2Table(plate_ref_bioactive);
bioactive(:,'profiles_1') = [];
bioactive.uniID = strcat(bioactive.plateIDs,bioactive.well_names)
cluster_table.uniID = strcat(cluster_table.plateIDs,cluster_table.well_names)
CMG_Table = join(bioactive,cluster_table,'Keys','uniID');
writetable(CMG_Table,[results_folder,'\','CMG_Table.xlsx'])
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
    
   
    % Plot 2
    ax(i) = subplot(2,1,2,'Parent',fh);
    %imagesc(percs(:,c_id_sig)),colormap(ax(i),hot),caxis(ax(i),[0,1])
    imagesc(percs)
    colormap(ax(i),hot),caxis(ax(i),[0,1])
    ax(i).YTick = 1:length(ref_cluster(i).categories);
    ax(i).YTickLabel = ref_cluster(i).categories;
    ax(i).TickLabelInterpreter = 'none';
    ax(i).XTick = [];
    ax(i) = subplot(2,1,1,'Parent',fh);
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



%% Classification Accuracy
% Pick out AW reference
%aw_ref = {'DMSO','DNA','HDAC','HSP90','MT','Proteasome','mTOR'}; % 6 classes
aw_ref = {'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR'}; %10 classes

ref_plate = {'reference_plate'};

% Calculate classification accuracy for reference plates
for c = 1%:4
plateaccuracy = cell(6,2);    
refplates = referenceplates(referenceplates.CellLine== pORACL{c},:);
for p=1:height(refplates)
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





















