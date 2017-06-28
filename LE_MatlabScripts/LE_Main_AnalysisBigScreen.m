% Louise Evans March 2017 Purpose: Perform analogous analysis to Charles'
% NBT paper for BigScreen.

%% Set Paths
addpath(genpath('W:\Lab_Analysis\common\plate_annotation')) %plate annotation from labtoolbox
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\CD_Tag_Drug_Screen')) % Charles code for NBT paper containing useful functions
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_analysis')) % Louise code afer cluster operations
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_Preprocessing')) % Louise code before cluster operations
addpath(genpath('W:\Lab_Analysis\common\image_browser')) % Image_Browser to callback images

%% Define Variables
code_folder = 'W:\2015_09_HTS_LE\Code_LE\';
profiles_folder = 'W:\2015_09_HTS_LE\data\profiles';% place to save profiles to and load from
profiles_folder_NBT = 'W:\2015_09_HTS_LE\data\profiles_NBT';% subset of features used in NBT paper
qc_folder = 'W:\2015_09_HTS_LE\QC\bigscreen';% store any QC analysis Python, Matlab, or manual
results_folder = 'W:\2015_09_HTS_LE\results\matlab_results\bigscreen\'; % Louise results for 4lines analysis
plate_DB_path ='W:\2015_09_HTS_LE\project_database\'; % directory containing plate database file
plate_DB = 'Plate_database_latest2.xlsx'; % file describing relationship between compound and experiment plates

%% Load Constants
Exp_DB = readtable(fullfile(plate_DB_path,plate_DB),'Sheet',1);

year = 2017;
theseyearplates = ismember(Exp_DB.year,year);
plates = Exp_DB(theseyearplates,:);

cellline = 'A549_XRCC5';
thesecellplates = ismember(plates.CellLine,cellline);
plates = plates(thesecellplates,:);

batch = '2017018_001';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);

batch='4L2K_7';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);

batch='4L2K_9';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);

batch='4L2K_1';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);

batch='4L2K_3';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);

batch='4L2K_5';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);

% batch='2017018_002';
% thesebatchplates = ~ismember(plates.batch,batch);
% plates = plates(thesebatchplates,:);
% 
% batch='2017018_003';
% thesebatchplates = ~ismember(plates.batch,batch);
% plates = plates(thesebatchplates,:);

% batch='2017018_004';
% thesebatchplates = ~ismember(plates.batch,batch);
% plates = plates(thesebatchplates,:);

batch='2017018_005';
thesebatchplates = ~ismember(plates.batch,batch);
plates = plates(thesebatchplates,:);

cpdorder='blank';
thesecpdplates = ~ismember(plates.Cpd_Addition_order,cpdorder);
plates = plates(thesecpdplates,:);

plate_type='carryover_ctrl';
thesetypeplates = ~ismember(plates.plate_type,plate_type);
plates = plates(thesetypeplates,:);

thesereferenceplates = ismember(plates.plate_type,'reference_plate');
referenceplates = plates(thesereferenceplates,:);

clear thesebatchplates thesereferenceplates thesecellplates thesetypeplates theseyearplates thesecpdplates

%% Generate Profiles
LE_Gen_ProfilesInCell(plates,profiles_folder)

% Subset Profiles to NBT features
[~,~,~,NBT] = import_featuresselection('W:\2015_09_HTS_LE\Code_LE\LE_Analysis\LE_MatlabScripts\feature_select_NBT_bigscreen.txt');   
for p = 1:height(plates)
    plate_ID = plates.expt_plate{p};
    load(fullfile(profiles_folder,strcat('plate_',plate_ID,'.mat')));
    for w = 1:size(plate.profiles,1)
        plate.profiles{w,1} = plate.profiles{w,1}(NBT);
    end
    plate.feaInfo.name = plate.feaInfo.name(NBT);
    plate.feaInfo.channel_id = plate.feaInfo.channel_id(NBT);
    plate.feaInfo.category_id = plate.feaInfo.category_id(NBT);
    save(fullfile(profiles_folder_NBT,strcat('plate_',plate_ID,'.mat')),'plate') 
end

%% Get Analysis Parameters
% Get plotting options
getplotoptions() % generate standard plot options
plot_opt = plot_opt_reference; %plot_opt = {'cate_to_show',{'DMSO','Proteasome'}}; % specific plot option
plot_opt{1,2}{1,12}= 'Unknown';

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

% now plates can be categorical
plates.CellLine = categorical(plates.CellLine);
referenceplates.CellLine = categorical(referenceplates.CellLine);
plates.CellLine = categorical(plates.CellLine);            
%getplotoptions()  

%% Specify Bioactives

plate_IDs = plates.expt_plate;
%plate_IDs = plates.expt_plate;
%plate_IDs = referenceplates.expt_plate;

plate = Load_Batch_LE(plates,plate_IDs, profiles_folder_NBT,qc_folder);
checknumplatewells()
plate = Merge_Time_Points(plate,time_to_use,merge_mode);
[plate_bioactive,plate3,plate_ref_bioactive,plate_ref,plate_cpd_bioactive,plate_cpd] = Get_Bioactive_Compounds(plate, paras);

plate_T=Plate2Table(plate);
plate_bioactive_T=Plate2Table(plate_bioactive);
plate2_T=Plate2Table(plate3);
plate_ref_T=Plate2Table(plate_ref);
plate_cpd_bioactive_T=Plate2Table(plate_cpd_bioactive);
plate_cpd_T=Plate2Table(plate_cpd);
plate_ref_bioactive_T = Plate2Table(plate_ref_bioactive);

writetable(plate_cpd_bioactive_T,strcat(results_folder,'screen_cpds_bioactive.xlsx'));

is_inactive = plate3.bioactive_pVals>=pVal_thr & (~strcmp(plate3.drug_categories,'DMSO'));
plate3.drug_categories(is_inactive) = {'Nonbioactive'}; % not using Nonbioactive to do the PCA.
bioactive = Plate2Table(plate3);
writetable(bioactive,fullfile(results_folder,'_bioactivepvals_NBT.xlsx'))
plot_opt{1,2}{1,12}='Unknown';
plot_opt_ALL{1,2}{1,35}='Unknown';
Visualize_Plate_LE(plate_bioactive,plot_opt{:}); % or plot_opt_query -
labelpcplot()
% 
%% Classification Accuracy
%aw_ref = {'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR','Unknown'};
aw_ref = {'DMSO','DNA','HDAC','HSP90','MT','Proteasome','mTOR'};

for p=1:height(referenceplates)
    % Load plate
    plate_IDs = referenceplates.expt_plate;
    plate_ID = plate_IDs(p);
    display(plate_ID);
    plate = Load_Batch_LE(referenceplates,plate_ID, profiles_folder_NBT,qc_folder);
    plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);
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
    classifier.accuracy{p} = accuracy;
    classifier.clfy_obj{p} = clfy_obj;
    classifier.cv_report{p} = cv_report;
    classifier.scan_scheme{p} = scan_scheme;
    save(fullfile(results_folder,'classifier.mat'),'classifier')
end

% Calc CA for multiple plates together
plate_IDs = referenceplates.expt_plate;
plate_IDs = referenceplates.expt_plate([1 3 4],:);
plate = Load_Batch_LE(referenceplates,plate_IDs, profiles_folder_NBT,qc_folder);
plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);
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
classifier.accuracy = accuracy;
classifier.clfy_obj = clfy_obj;
classifier.cv_report = cv_report;
classifier.scan_scheme = scan_scheme;
save(fullfile(results_folder,'classifier_all.mat'),'classifier')
    
%% Visualize FDA space
aw_ref_extended = {'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR','Unknown'};
plate_IDs = plates.expt_plate;
plate = Load_Batch_LE(plates,plate_IDs, profiles_folder_NBT,qc_folder);
plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);
plate2 = Get_Bioactive_Compounds(plate2, paras);
plate2 = Select_Rows(plate2,ismember(plate2.drug_categories,aw_ref_extended));
plate2 = Select_Rows(plate2,ismember(plate2.plate_types,{'reference_plate'}));
plate_fda = plate2;
plate_fda.profiles = cellfun(@clfy_obj.KLFDA_proj_fun,plate2.profiles,'unif',false);
Visualize_Plate(plate_fda,plot_opt_ALL{:});
saveas(gcf,fullfile(results_folder,'FDA'));

%% PREDICTIONS
% Predict Output
plate_IDs = plates.expt_plate;
plate = Load_Batch_LE(plates,plate_IDs, profiles_folder_NBT,qc_folder);
predict_output = Screen_Plates(plate, paras);
save(fullfile(results_folder,'predict_output.mat'),'predict_output')

%% Clustering
% Clustering parameters
FDA_type            = 0; % 0 FDA | -1 no transformation | -2 RCA
gamma               = 0:0.1:1; % -1: use theoretical value to regularize the within-class Variance matrix for LDA
linkage_method      = 'average'; % 'average' | 'centroid' | 'complete' | 'single'
D_within_prob_thr   = 0.5; % probability cutoff for choosing diStance threshold to define clusters

col_order           = {'Cluster','Size','Significance','plateIDs','well_names','drug_names','drug_categories',...
                       'concentrations','bioactive_pVals','targets','pathways','descriptions','cas','compound',...
                       'cpd_usage','nCells','clone_names','plate_types'};
                   
sort_order          = {'Size','Cluster','pathways','targets','drug_names','concentrations'};

paras               = Set_Paras('pVal_thr',pVal_thr,'classifier_type',classifier_type,...
                                'FDA_type',FDA_type,...
                                'gamma',gamma,'kernel_type',kernel_type,'kernel_para',kernel_para,...
                                'classifier_paras',classifier_paras,'kfold',kfold,...
                                'linkage_method',linkage_method,'D_within_prob_thr',D_within_prob_thr);
batches{1,1} = plate_IDs;
Nbatches = length(batches);
cluster_info = cell(Nbatches,1);
profiles = cell(Nbatches,1);
tree = cell(Nbatches,1);

for b = 1:Nbatches
    disp(b)
    % 1. Get plates
    plate_IDs = plates.expt_plate;
    plate = Load_Batch_LE(plates,plate_IDs, profiles_folder_NBT,qc_folder);
    % 2. Clustering
    [cluster_info{b}, profiles{b}, tree{b}] = Cluster_Compounds(plate_bioactive,paras);
    cluster_info{b}.compound = cluster_info{b}.drug_names;
    cluster_info{b} = cluster_info{b}(:,col_order);
end

% Visualize clusters
batch_to_plot = 1;
clusters_to_color = [];

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
%saveas(gcf,fullfile(results_folder,batch));

%% Analysis of clusters
fh = figure;
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
end

%% Louise Percents
% for i = 1:length(cluster_info)
%     query_cluster(i).categories = pathwaysize(:,1)
%     query_cluster(i).cluster_ID = unique(cluster_info{i}.Cluster);
%     query_cluster(i).stat = zeros(length(query_cluster(i).categories),length(query_cluster(i).cluster_ID));
%     for j = 1:size(query_cluster(i).stat,1)
%         is_class = strcmpi(cluster_info{i}.drug_categories,query_cluster(i).categories{j});
%         tb = tabulate(cluster_info{i}.Cluster(is_class));
%
%         [is_in,loc] = ismember(tb(:,1),query_cluster(i).cluster_ID);
%         query_cluster(i).stat(j,loc(is_in)) = tb(:,2);
%     end
% end

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

% scatter(x_ax,y_ax)
i =1;

fh = figure;
ax(i) = subplot(1,1,i,'Parent',fh);
imagesc(x)
colormap(ax(i),'hot');
colorbar(ax(i));
ax(i).YTick = 1:length(ref_cluster(i).categories);
ax(i).YTickLabel = ref_cluster(i).categories;
ax(i).TickLabelInterpreter = 'none';
xlabel('Cluster #');
ylabel('Compound Category (Pathway)');
title('Percentage Drugs in a Compound Category Appearing in each Cluster');
print('-f1','XRCC5clusterper','-dpng','-r500');

%% Print Physical Plate Cell Number

cd (results_folder) % move to the location to do all the saving in to improve speed and make checking for files quicker
for b = 55:132 %plate IDs last 3 digits
    pid = b; % first check blanks, ref plates, compound plates
    plateIDs = arrayfun(@(x)sprintf('2017018%03d',x),pid,'unif',false)'; % did not trouble to change to my plate_IDs method for one-off task
    filename =char(strcat(plateIDs(1),'.bmp'));
    display(filename)
    if exist(filename,'file')==0
        try % in case file is missing
            plate = Load_Batch_LE(plates,plate_IDs, profiles_folder_NBT,qc_folder); % load plate
            Show_Cell_Number(plate); % generate heatmap of cell number
            saveas(gcf,plateIDs{1,1},'bmp'); % save
        catch
            warning (char(strcat(filename, ' cannot be generated'))) %
        end
    else
            warning (char(strcat(filename, ' already exists')))
    end
end
cd (code_folder)

%% Print Physical Plate P Values

cd (results_folder) % move to the location to do all the saving in to improve speed and make checking for files quicker
plate_IDs = plates.expt_plate;
for p = 1:length(plate_IDs)
    plate = Load_Batch_LE(plates,plate_IDs(p), profiles_folder_NBT,qc_folder);
    plate2 = Merge_Time_Points(plate,time_to_use,merge_mode);
    [plate2,~] = Get_Bioactive_Compounds(plate2, paras);
    plate2.logbioactive_pVals = log(plate2.bioactive_pVals);
    plate2.logbioactive_pVals(isnan(plate2.logbioactive_pVals))=1;%A(isnan(A)) = 0
    Show_PVALUE(plate2)
    title = [plate_IDs{p},' Log P Value Plate'];
    print('-f1',title,'-dpng','-r500')
    close all
end
cd (code_folder)

%% P Value QC

pvalue_plotter % plot P values by cell line, and by batch
