% Louise Heinrich 2017.06.30 in Altschuler and Wu Labs
% Analysis of data generated for EMT / compound screen / phenopushing

%% Set Paths
addpath(genpath('W:\Lab_Analysis\common\plate_annotation')) %plate annotation from labtoolbox
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\CD_Tag_Drug_Screen')) % Charles code for NBT paper containing useful functions
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_analysis')) % Louise code afer cluster operations
addpath(genpath('W:\2015_09_HTS_LE\Code_LE\LE_Preprocessing')) % Louise code before cluster operations

%% Define Variables
code_folder = 'W:\2015_09_HTS_LE\Code_LE';
profiles_folder = 'W:\2015_09_HTS_LE\data\profiles';% place to save profiles to and load from
profiles_folder_NBT = 'W:\2015_09_HTS_LE\data\profiles_NBT';% subset of features used in NBT paper
qc_folder = 'W:\2015_09_HTS_LE\QC\4lines\stringent'; % QC folder for 4lines data
qc_folder_fibrocl = 'W:\2015_09_HTS_LE\QC\fibrocl';% store any QC analysis Python, Matlab, or manual
results_folder = 'W:\2015_09_HTS_LE\results\matlab_results\fibrocl'; % Louise results for 4lines analysis
plate_DB_path ='W:\2015_09_HTS_LE\project_database\'; % directory containing plate database file
plate_DB = 'Plate_database_latest2.xlsx'; % file describing relationship between compound and experiment plates

%% Load Constants
Exp_DB = readtable(fullfile(plate_DB_path,plate_DB),'Sheet',1);
load('fibrocl_DB.mat');

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
pORACL_plates = plates(thesebatchplates,:);

pORACL = {'A549_XRCC5';'A549_NQO1';'A549_SET';'A549_S100A11',};
clear theseyearplates thesebatchplates thesereferenceplates year batch expt plate_DB_path plate_DB

expt = 'LE_6';
thesebatchplates = ismember(Exp_DB.Experiment_,expt);
fibrocl_plates_NQO1 = Exp_DB(thesebatchplates,:);

expt = 'LE_C_3';
thesebatchplates = ismember(Exp_DB.Experiment_,expt);
fibrocl_plates = Exp_DB(thesebatchplates,:);

expt = '#LE_C_3';
thesebatchplates = ismember(Exp_DB.Experiment_,expt);
fibrocl_plates2 = Exp_DB(thesebatchplates,:);

fibrocl_plates = vertcat(fibrocl_plates,fibrocl_plates2);
clear fibrocl_plates2 expt thesebatchplates