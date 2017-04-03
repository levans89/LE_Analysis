% set categorical
SET_bioactive.cpd_usage = categorical(SET_bioactive.cpd_usage);
NQO1_bioactive.cpd_usage = categorical(NQO1_bioactive.cpd_usage);
S100A11_bioactive.cpd_usage = categorical(S100A11_bioactive.cpd_usage);
XRCC5_bioactive.cpd_usage = categorical(XRCC5_bioactive.cpd_usage);

% name P Values by cell line
XRCC5_bioactive.Properties.VariableNames{7} = 'XRCC5_bioactive_pVals'; 
SET_bioactive.Properties.VariableNames{7} = 'SET_bioactive_pVals';
S100A11_bioactive.Properties.VariableNames{7} = 'S100A11_bioactive_pVals';
NQO1_bioactive.Properties.VariableNames{7} = 'NQO1_bioactive_pVals';

% define subsets of compounds
SET_b_query = SET_bioactive(SET_bioactive.cpd_usage =='query_cpd',:);
SET_b_reference = SET_bioactive(SET_bioactive.cpd_usage =='reference_cpd',:);
SET_b_remainders = SET_bioactive(SET_bioactive.cpd_usage ~='reference_cpd' & SET_bioactive.cpd_usage ~='query_cpd',:);
XRCC5_b_query = XRCC5_bioactive(XRCC5_bioactive.cpd_usage =='query_cpd',:);
XRCC5_b_reference = XRCC5_bioactive(XRCC5_bioactive.cpd_usage =='reference_cpd',:);
XRCC5_b_remainders = XRCC5_bioactive(XRCC5_bioactive.cpd_usage ~='reference_cpd' & XRCC5_bioactive.cpd_usage ~='query_cpd',:);
S100A11_b_query = S100A11_bioactive(S100A11_bioactive.cpd_usage =='query_cpd',:);
S100A11_b_reference = S100A11_bioactive(S100A11_bioactive.cpd_usage =='reference_cpd',:);
S100A11_b_remainders = S100A11_bioactive(S100A11_bioactive.cpd_usage ~='reference_cpd' & S100A11_bioactive.cpd_usage ~='query_cpd',:);
NQO1_b_query = NQO1_bioactive(NQO1_bioactive.cpd_usage =='query_cpd',:);
NQO1_b_reference = NQO1_bioactive(NQO1_bioactive.cpd_usage =='reference_cpd',:);
NQO1_b_remainders = NQO1_bioactive(NQO1_bioactive.cpd_usage ~='reference_cpd' & NQO1_bioactive.cpd_usage ~='query_cpd',:);

% Calculate Bioactive Percentages
percents.bioactive.NQO1 = sum(NQO1_b_query.NQO1_bioactive_pVals <= pVal_thr)/height(NQO1_b_query)*100;
percents.bioactive.SET = sum(SET_b_query.SET_bioactive_pVals <= pVal_thr)/height(SET_b_query)*100;
percents.bioactive.XRCC5 = sum(XRCC5_b_query.XRCC5_bioactive_pVals <= pVal_thr)/height(XRCC5_b_query)*100;
percents.bioactive.S100A11 = sum(S100A11_b_query.S100A11_bioactive_pVals <= pVal_thr)/height(S100A11_b_query)*100;
percents.ref.NQO1 = sum(NQO1_b_reference.NQO1_bioactive_pVals <= pVal_thr)/height(NQO1_b_reference)*100;
percents.ref.SET = sum(SET_b_reference.SET_bioactive_pVals <= pVal_thr)/height(SET_b_reference)*100;
percents.ref.XRCC5 = sum(XRCC5_b_reference.XRCC5_bioactive_pVals <= pVal_thr)/height(XRCC5_b_reference)*100;
percents.ref.S100A11 = sum(S100A11_b_reference.S100A11_bioactive_pVals <= pVal_thr)/height(S100A11_b_reference)*100;
a = struct2table(percents.ref);
b = struct2table(percents.bioactive);
percents = vertcat(a,b);
percents.Properties.RowNames{1}='Reference';
percents.Properties.RowNames{2}='Selleck 2K';
display(percents)
%writetable(percents,fullfile(results_folder,'percents.xlsx'))

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

% Line Plot P Values for Query Compounds 2K Selleck 
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
xlim([1 1835])
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
    
%% plot log p value by pathways
close all
bioactives = sortrows(bioactives,'drug_names','ascend'); % reset for consistent figures

bioactives = sortrows(bioactives,{'Pathway','Target'},{'ascend','ascend'});
bioactives1 = bioactives{:,4:7};
bioactives1(isnan(bioactives1))=1; %A(isnan(A)) = 0
%LE_heatmap_pvalues(bioactives1,'Sort by Pathway, Target','P Value');
%saveas(gca,fullfile(results_folder,'PValues_PathwayTarget'));
%print('-f1',fullfile(results_folder,'PValues_PathwayTarget'),'-dpng','-r500');
%close all

bioactives2 = log(bioactives1);
LE_heatmap_pvalues(bioactives2,'Sort by Pathway,Target','Log P Values');
%saveas(gca,fullfile(results_folder,'LogPValues_PathwayTarget'));
%print('-f1',fullfile(results_folder,'LogPValues_PathwayTarget'),'-dpng','-r500');
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




clear bioactives1 bioactives2 bioactives3 bioactives4 bioactives5 bioactives6