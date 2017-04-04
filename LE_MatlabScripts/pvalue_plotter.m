%% Set P Val Thr
pVal_thr            = 10^-6;
%% set categorical
SET_bioactive.cpd_usage = categorical(SET_bioactive.cpd_usage);
NQO1_bioactive.cpd_usage = categorical(NQO1_bioactive.cpd_usage);
S100A11_bioactive.cpd_usage = categorical(S100A11_bioactive.cpd_usage);
XRCC5_bioactive.cpd_usage = categorical(XRCC5_bioactive.cpd_usage);

%% name P Values by cell line
XRCC5_bioactive.Properties.VariableNames{7} = 'XRCC5_bioactive_pVals'; 
SET_bioactive.Properties.VariableNames{7} = 'SET_bioactive_pVals';
S100A11_bioactive.Properties.VariableNames{7} = 'S100A11_bioactive_pVals';
NQO1_bioactive.Properties.VariableNames{7} = 'NQO1_bioactive_pVals';

%% define subsets of compounds
SET_b_query = SET_bioactive(SET_bioactive.cpd_usage =='query_cpd',:);
SET_b_reference = SET_bioactive(SET_bioactive.cpd_usage =='reference_cpd',:);
%SET_b_remainders = SET_bioactive(SET_bioactive.cpd_usage ~='reference_cpd' & SET_bioactive.cpd_usage ~='query_cpd',:);
XRCC5_b_query = XRCC5_bioactive(XRCC5_bioactive.cpd_usage =='query_cpd',:);
XRCC5_b_reference = XRCC5_bioactive(XRCC5_bioactive.cpd_usage =='reference_cpd',:);
%XRCC5_b_remainders = XRCC5_bioactive(XRCC5_bioactive.cpd_usage ~='reference_cpd' & XRCC5_bioactive.cpd_usage ~='query_cpd',:);
S100A11_b_query = S100A11_bioactive(S100A11_bioactive.cpd_usage =='query_cpd',:);
S100A11_b_reference = S100A11_bioactive(S100A11_bioactive.cpd_usage =='reference_cpd',:);
%S100A11_b_remainders = S100A11_bioactive(S100A11_bioactive.cpd_usage ~='reference_cpd' & S100A11_bioactive.cpd_usage ~='query_cpd',:);
NQO1_b_query = NQO1_bioactive(NQO1_bioactive.cpd_usage =='query_cpd',:);
NQO1_b_reference = NQO1_bioactive(NQO1_bioactive.cpd_usage =='reference_cpd',:);
%NQO1_b_remainders = NQO1_bioactive(NQO1_bioactive.cpd_usage ~='reference_cpd' & NQO1_bioactive.cpd_usage ~='query_cpd',:);
clc

%% Calculate Bioactive Percentages
% 1044 = max # ref wells (6 plates)
% height(compound_IDs) = max # query wells (Selleck 2K library)

percents.bioactive.NQO1 = sum(NQO1_b_query.bioactive_pVals <= pVal_thr)/height(compound_IDs)*100; %height(NQO1_b_query)*100;
percents.bioactive.SET = sum(SET_b_query.bioactive_pVals <= pVal_thr)/height(compound_IDs)*100; %height(SET_b_query)*100;
percents.bioactive.XRCC5 = sum(XRCC5_b_query.bioactive_pVals <= pVal_thr)/height(compound_IDs)*100; %height(XRCC5_b_query)*100;
percents.bioactive.S100A11 = sum(S100A11_b_query.bioactive_pVals <= pVal_thr)/height(compound_IDs)*100; %height(S100A11_b_query)*100;
percents.ref.NQO1 = sum(NQO1_b_reference.bioactive_pVals <= pVal_thr)/1044*100;%height(NQO1_b_reference)*100;
percents.ref.SET = sum(SET_b_reference.bioactive_pVals <= pVal_thr)/1044*100;%height(SET_b_reference)*100;
percents.ref.XRCC5 = sum(XRCC5_b_reference.bioactive_pVals <= pVal_thr)/1044*100;%height(XRCC5_b_reference)*100;
percents.ref.S100A11 = sum(S100A11_b_reference.bioactive_pVals <= pVal_thr)/1044*100;%height(S100A11_b_reference)*100;
%percentsunion = length(unionqueryhits)/height(compound_IDs)*100;
a = struct2table(percents.ref);
b = struct2table(percents.bioactive);
percents = vertcat(a,b);
percents.Properties.RowNames{1}='Reference';
percents.Properties.RowNames{2}='Selleck 2K';
display(percents)
writetable(percents,fullfile(results_folder,'percents_updated.xlsx'))
%% calculate intersects and unions of query hits
% get hit compounds
SET_hits = SET_b_query(SET_b_query.bioactive_pVals<=pVal_thr,:);
S100A11_hits = S100A11_b_query(S100A11_b_query.bioactive_pVals<=pVal_thr,:);
NQO1_hits = NQO1_b_query(NQO1_b_query.bioactive_pVals<=pVal_thr,:);
XRCC5_hits = XRCC5_b_query(XRCC5_b_query.bioactive_pVals<=pVal_thr,:);

% get hit reference compounds
SET_Rhits = SET_b_reference(SET_b_reference.bioactive_pVals<=pVal_thr,:);
S100A11_Rhits = S100A11_b_reference(S100A11_b_reference.bioactive_pVals<=pVal_thr,:);
NQO1_Rhits = NQO1_b_reference(NQO1_b_reference.bioactive_pVals<=pVal_thr,:);
XRCC5_Rhits = XRCC5_b_reference(XRCC5_b_reference.bioactive_pVals<=pVal_thr,:);

% calculate bioactive numbers
counts.XRCC5.BA_n_Q = height(XRCC5_hits);
counts.XRCC5.BA_n_R = height(XRCC5_Rhits);
counts.NQO1.BA_n_Q = height(NQO1_hits);
counts.NQO1.BA_n_R = height(NQO1_Rhits);
counts.SET.BA_n_Q = height(SET_hits);
counts.SET.BA_n_R = height(SET_Rhits);
counts.S100A11.BA_n_Q = height(S100A11_hits);
counts.S100A11.BA_n_R = height(S100A11_Rhits);

% calculate bioactive percentage
counts.XRCC5.BA_per_Q = height(XRCC5_hits)/height(compound_IDs)*100;
counts.XRCC5.BA_per_R = height(XRCC5_Rhits)/1044*100;
counts.NQO1.BA_per_Q = height(NQO1_hits)/height(compound_IDs)*100;
counts.NQO1.BA_per_R = height(NQO1_Rhits)/1044*100;
counts.SET.BA_per_Q = height(SET_hits)/height(compound_IDs)*100;
counts.SET.BA_per_R = height(SET_Rhits)/1044*100;
counts.S100A11.BA_per_Q = height(S100A11_hits)/height(compound_IDs)*100;
counts.S100A11.BA_per_R = height(S100A11_Rhits)/1044*100;

% calculate union between 4lines
unionXN = union(XRCC5_hits.drug_names,NQO1_hits.drug_names);
unionSS = union(SET_hits.drug_names,S100A11_hits.drug_names);
unionqueryhits = cell2table(union(unionXN,unionSS));
clear overallhits2 overallhits

% calculate intersection between 4lines
interXN = intersect(XRCC5_hits.drug_names,NQO1_hits.drug_names);
interXSE=intersect(XRCC5_hits.drug_names,SET_hits.drug_names);
interXS1=intersect(XRCC5_hits.drug_names,S100A11_hits.drug_names);
interSES1=intersect(SET_hits.drug_names,S100A11_hits.drug_names);

% calculate union betwen XRCC5 and individual lines
counts.XRCC5.union = cell2table(union(XRCC5_hits.drug_names,XRCC5_hits.drug_names));
counts.NQO1.union = cell2table(union(XRCC5_hits.drug_names,NQO1_hits.drug_names));
counts.SET.union = cell2table(union(XRCC5_hits.drug_names,SET_hits.drug_names));
counts.S100A11.union = cell2table(union(XRCC5_hits.drug_names,S100A11_hits.drug_names));
counts.all.union = unionqueryhits;

% calculate union numbers
counts.XRCC5.UnXhits= height(counts.XRCC5.union);
counts.SET.UnXhits=height(counts.SET.union);
counts.S100A11.UnXhits=height(counts.S100A11.union);
counts.NQO1.UnXhits=height(counts.NQO1.union);
counts.all.UnXhits = height(counts.all.union);

% calculate union percents as percent of XRCC5 hits found 
counts.XRCC5.UperXhits= height(counts.XRCC5.union)/height(XRCC5_hits)*100;
counts.SET.UperXhits=height(counts.SET.union)/height(XRCC5_hits)*100;
counts.S100A11.UperXhits=height(counts.S100A11.union)/height(XRCC5_hits)*100;
counts.NQO1.UperXhits=height(counts.NQO1.union)/height(XRCC5_hits)*100;
counts.all.UperXhits = height(counts.all.union)/height(XRCC5_hits)*100;

% calculate union numbers
counts.XRCC5.nQhits=height(counts.XRCC5.union);
counts.SET.nQhits=height(counts.SET.union);
counts.S100A11.nQhits=height(counts.S100A11.union);
counts.NQO1.nQhits = height(counts.NQO1.union);
counts.all.nQhits = height(counts.all.union);

% calculate union percents as percent of Query compound library
counts.XRCC5.UperQhits=height(counts.XRCC5.union)/height(compound_IDs)*100;
counts.SET.UperQhits=height(counts.SET.union)/height(compound_IDs)*100;
counts.S100A11.UperQhits=height(counts.S100A11.union)/height(compound_IDs)*100;
counts.NQO1.UperQhits = height(counts.NQO1.union)/height(compound_IDs)*100;
counts.all.UperQhits = height(counts.all.union)/height(compound_IDs)*100;

% clear unions
counts.XRCC5.union = 0;
counts.SET.union = 0;
counts.S100A11.union = 0;
counts.NQO1.union = 0;
counts.all.union = 0;
% make table to print
tbXRCC5 = struct2table(counts.XRCC5);
tbSET = struct2table(counts.SET);
tbS100A11 = struct2table(counts.S100A11);
tbNQO1 = struct2table(counts.NQO1);
tball = struct2table(counts.all);

vertcat(tbXRCC5,tbSET,tbS100A11,tbNQO1)
%% reshaping and get all into one table
compound_IDs.Properties.VariableNames{1} = 'drug_names';
bioactives = outerjoin(compound_IDs, SET_b_query(:,[4,7]), 'Type','Left');
bioactives.Properties.VariableNames{1} = 'drug_names';
bioactives = outerjoin(bioactives, S100A11_b_query(:,[4,7]),'Type','Left');
bioactives.Properties.VariableNames{1} = 'drug_names';
bioactives = outerjoin(bioactives, NQO1_b_query(:,[4,7]),'Type','Left');
bioactives.Properties.VariableNames{1} = 'drug_names';
bioactives = outerjoin(bioactives, XRCC5_b_query(:,[4,7]),'Type','Left');
bioactives.Properties.VariableNames{1} = 'drug_names';
bioactives(:,{'drug_names_right','drug_names_right_1','drug_names_right_2','drug_names_right_3'}) = [];

%% Line Plot P Values for Query Compounds 2K Selleck 
bioactives = sortrows(bioactives,'S100A11_bioactive_pVals','descend');
figure()
plot(bioactives.S100A11_bioactive_pVals,'LineWidth',1)
hold on
bioactives = sortrows(bioactives,'NQO1_bioactive_pVals','descend');
plot(bioactives.NQO1_bioactive_pVals,'LineWidth',1)
bioactives = sortrows(bioactives,'XRCC5_bioactive_pVals','descend');
plot(bioactives.XRCC5_bioactive_pVals,'LineWidth',1)
bioactives = sortrows(bioactives,'SET_bioactive_pVals','descend');
plot(bioactives.SET_bioactive_pVals,'LineWidth',1)
legend('S100A11','NQO1','XRCC5','SET');
xlabel('P Values All Bioactive');
ylabel('P Value');
set(gca,'YScale','log')
xlim([1 height(compound_IDs)])
ylim([10E-19 1])
%saveas(gca,fullfile(results_folder,'PVALS_QUERY'))
%print('-f1',fullfile(results_folder,'PVALS_QUERY'),'-dpng','-r500');
close all

%% Plot P Values for all Compounds across the dataset
figure()
hold on

S100A11_bioactive = sortrows(S100A11_bioactive,'S100A11_bioactive_pVals','descend');
plot(S100A11_bioactive.S100A11_bioactive_pVals,'LineWidth',1)

NQO1_bioactive = sortrows(NQO1_bioactive,'NQO1_bioactive_pVals','descend');
plot(NQO1_bioactive.NQO1_bioactive_pVals,'LineWidth',1)

XRCC5_bioactive = sortrows(XRCC5_bioactive,'XRCC5_bioactive_pVals','descend');
plot(XRCC5_bioactive.XRCC5_bioactive_pVals,'LineWidth',1)

SET_bioactive = sortrows(SET_bioactive,'SET_bioactive_pVals','descend');
plot(SET_bioactive.SET_bioactive_pVals,'LineWidth',1)

legend('S100A11','NQO1','XRCC5','SET');
xlabel('P Values All Bioactive');
ylabel('P Value');
set(gca,'YScale','log')
xlim([1 max(endrows)])
ylim([10E-19 1])
%saveas(gca,fullfile(results_folder,'PVALS_ALL'))
%print('-f1',fullfile(results_folder,'PVALS_ALL'),'-dpng','-r500');
close all

%% Plot P Values per Plate
close all
% Select pORACL
pORACL = 'A549_NQO1';
plateIDs = currexp(currexp.CellLine == pORACL & currexp.batch ~= '4L2K_7' & currexp.batch ~= '4L2K_8',:);
batched = horzcat(plateIDs.expt_plate,plateIDs.batch,plateIDs.plate_type);
for a = 1:height(plateIDs);
    NQO1_bioactive.plateIDs = categorical(NQO1_bioactive.plateIDs);
    pVals = NQO1_bioactive(NQO1_bioactive.plateIDs == plateIDs.expt_plate{a,1},:);
    pVals = sortrows(pVals,'NQO1_bioactive_pVals','ascend');
    plot(pVals{:,7})
    hold on
end
legend('Ref B2','Screen B2','Screen B2','Ref B2','Ref B4','Screen B4','Screen B4','Ref B4','Ref B6','Screen B6','Screen B6','Ref B6')
title('P Value by Plate NQO1')
xlabel('Compound')
ylabel('P Value')
hold off
    
pORACL = 'A549_SET'
plateIDs = currexp(currexp.CellLine == pORACL & currexp.batch ~= '4L2K_7' & currexp.batch ~= '4L2K_8',:)
batched = horzcat(plateIDs.expt_plate,plateIDs.batch,plateIDs.plate_type)
display(batched)
for a = 1:height(plateIDs)
    SET_bioactive.plateIDs = categorical(SET_bioactive.plateIDs);
    pVals = SET_bioactive(SET_bioactive.plateIDs == plateIDs.expt_plate{a,1},:);
    pVals = sortrows(pVals,'SET_bioactive_pVals','ascend');
    plot(pVals{:,7})
    hold on
end
legend('Ref B2','Screen B2','Screen B2','Ref B2','Ref B4','Screen B4','Screen B4','Ref B4','Ref B6','Screen B6','Screen B6','Ref B6')
title('P Value by Plate SET')
xlabel('Compound')
ylabel('P Value')
hold off

pORACL = 'A549_S100A11'
plateIDs = currexp(currexp.CellLine == pORACL & currexp.batch ~= '4L2K_7' & currexp.batch ~= '4L2K_8',:)
batched = horzcat(plateIDs.expt_plate,plateIDs.batch,plateIDs.plate_type)
display(batched)
for a = 1:height(plateIDs)
    S100A11_bioactive.plateIDs = categorical(S100A11_bioactive.plateIDs);
    pVals = S100A11_bioactive(S100A11_bioactive.plateIDs == plateIDs.expt_plate{a,1},:);
    pVals = sortrows(pVals,'S100A11_bioactive_pVals','ascend');
    plot(pVals{:,7})
    hold on
end
legend('Ref B1','Screen B1','Screen B1','Ref B1','Ref B3','Screen B3','Screen B3','Ref B3','Ref B5','Screen B5','Screen B5','Ref B5')
title('P Value by Plate S100A11')
xlabel('Compound')
ylabel('P Value')
hold off

pORACL = 'A549_XRCC5'
plateIDs = currexp(currexp.CellLine == pORACL & currexp.batch ~= '4L2K_7' & currexp.batch ~= '4L2K_8',:)
batched = horzcat(plateIDs.expt_plate,plateIDs.batch,plateIDs.plate_type)
display(batched)
for a = 1:height(plateIDs)
    XRCC5_bioactive.plateIDs = categorical(XRCC5_bioactive.plateIDs);
    pVals = XRCC5_bioactive(XRCC5_bioactive.plateIDs == plateIDs.expt_plate{a,1},:);
    pVals = sortrows(pVals,'XRCC5_bioactive_pVals','ascend');
    plot(pVals{:,7})
    hold on
end
legend('Ref B1','Screen B1','Screen B1','Ref B1','Ref B3','Screen B3','Screen B3','Ref B3','Ref B5','Screen B5','Screen B5','Ref B5')
title('P Value by Plate XRCC5')
xlabel('Compound')
ylabel('P Value')
hold off
    
%% Plot Log P by Pathway etc.
close all
bioactives = sortrows(bioactives,'drug_names','ascend'); % reset for consistent figures

bioactives = sortrows(bioactives,{'Pathway','XRCC5_bioactive_pVals','NQO1_bioactive_pVals'},{'ascend','ascend','ascend'});
bioactives1 = bioactives{:,4:7};
bioactives1(isnan(bioactives1))=1; %A(isnan(A)) = 0
%LE_heatmap_pvalues(bioactives1,'Sort by Pathway, Target','P Value');
%saveas(gca,fullfile(results_folder,'PValues_PathwayTarget'));
%print('-f1',fullfile(results_folder,'PValues_PathwayTarget'),'-dpng','-r500');
%close all

bioactives2 = log(bioactives1);
createfigure(bioactives2);
%LE_heatmap_pvalues(bioactives2,'Sort by Pathway,Target','Log P Values');
%saveas(gca,fullfile(results_folder,'LogPValues_PathwayTarget'));
%print('-f1',fullfile(results_folder,'LogPValues_Pathway_X_N'),'-dpng','-r500');
%close all


bioactives = sortrows(bioactives,'drug_names','ascend'); % reset for consistent figures
bioactives = sortrows(bioactives,{'Pathway','XRCC5_bioactive_pVals','NQO1_bioactive_pVals','SET_bioactive_pVals',...
                                    'S100A11_bioactive_pVals'},...
                                    {'ascend','ascend','ascend','ascend','ascend'});
bioactives3 = bioactives{:,4:7};
bioactives3(isnan(bioactives3))=1;
bioactives4 = log(bioactives3);
%LE_heatmap_pvalues(bioactives3,'Sort by Pathway, XRCC5','P Values')
%saveas(gca,fullfile(results_folder,'PValues_PathwayXRCC5'))
%print('-f1',fullfile(results_folder,'PValues_PathwayXRCC5'),'-dpng','-r500');
%close all
LE_heatmap_pvalues(bioactives4,'Sort by Pathway, XRCC5','Log P Values');
%saveas(gca,fullfile(results_folder,'LogPValues_PathwayXRCC5'))
%print('-f1',fullfile(results_folder,'LogPValues_PathwayXRCC5'),'-dpng','-r500');
%close all

bioactives = sortrows(bioactives,'drug_names','ascend'); % reset for consistent figures
bioactives = sortrows(bioactives,{'XRCC5_bioactive_pVals','Pathway','NQO1_bioactive_pVals','SET_bioactive_pVals',...
    'S100A11_bioactive_pVals'},...
    {'ascend','ascend','ascend','ascend','ascend'});
bioactives5 = bioactives{:,4:7};
bioactives5(isnan(bioactives5))=1;
bioactives6 = log(bioactives5);
%LE_heatmap_pvalues(bioactives5,'Sort by Pathway, XRCC5','P Values')
%saveas(gca,fullfile(results_folder,'PValues_PathwayXRCC5'))
%print('-f1',fullfile(results_folder,'PValues_PathwayXRCC5'),'-dpng','-r500');
%close all
LE_heatmap_pvalues(bioactives6,'Sort by XRCC5, Pathway, pORACLs','Log P Values');
saveas(gca,fullfile(results_folder,'LogPValues_XRCC5Pathway'))
print('-f1',fullfile(results_folder,'LogPValues_XRCC5Pathway'),'-dpng','-r500');
close all
