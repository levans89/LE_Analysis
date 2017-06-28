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

%% Calculate Bioactive Percentages
percents.bioactive.NQO1 = sum(NQO1_b_query.NQO1_bioactive_pVals <= pVal_thr)/height(compound_IDs)*100; 
percents.bioactive.SET = sum(SET_b_query.SET_bioactive_pVals <= pVal_thr)/height(compound_IDs)*100; 
percents.bioactive.XRCC5 = sum(XRCC5_b_query.XRCC5_bioactive_pVals <= pVal_thr)/height(compound_IDs)*100; 
percents.bioactive.S100A11 = sum(S100A11_b_query.S100A11_bioactive_pVals <= pVal_thr)/height(compound_IDs)*100;

percents.ref.NQO1 = sum(NQO1_b_reference.NQO1_bioactive_pVals <= pVal_thr)/height(NQO1_b_reference)*100;
percents.ref.SET = sum(SET_b_reference.SET_bioactive_pVals <= pVal_thr)/height(SET_b_reference)*100;
percents.ref.XRCC5 = sum(XRCC5_b_reference.XRCC5_bioactive_pVals <= pVal_thr)/height(XRCC5_b_reference)*100;
percents.ref.S100A11 = sum(S100A11_b_reference.S100A11_bioactive_pVals <= pVal_thr)/height(S100A11_b_reference)*100;

percents_summary = vertcat(struct2table(percents.ref),struct2table(percents.bioactive));
percents_summary.Properties.RowNames{1}='Reference';
percents_summary.Properties.RowNames{2}='Selleck 2K';

display(percents_summary)
writetable(percents_summary,fullfile(results_folder,'percents_bwDMSO_NBT.xlsx'))

%% calculate intersects and unions of query hits
% get hit compounds
SET_Qhits = SET_b_query(SET_b_query.SET_bioactive_pVals<=pVal_thr,:);
S100A11_Qhits = S100A11_b_query(S100A11_b_query.S100A11_bioactive_pVals<=pVal_thr,:);
NQO1_Qhits = NQO1_b_query(NQO1_b_query.NQO1_bioactive_pVals<=pVal_thr,:);
XRCC5_Qhits = XRCC5_b_query(XRCC5_b_query.XRCC5_bioactive_pVals<=pVal_thr,:);

% get hit reference compounds
SET_Rhits = SET_b_reference(SET_b_reference.SET_bioactive_pVals<=pVal_thr,:);
S100A11_Rhits = S100A11_b_reference(S100A11_b_reference.S100A11_bioactive_pVals<=pVal_thr,:);
NQO1_Rhits = NQO1_b_reference(NQO1_b_reference.NQO1_bioactive_pVals<=pVal_thr,:);
XRCC5_Rhits = XRCC5_b_reference(XRCC5_b_reference.XRCC5_bioactive_pVals<=pVal_thr,:);

% calculate bioactive numbers
counts = table;
counts{1,1} = height(XRCC5_Qhits);
counts{2,1} = height(XRCC5_Rhits);
counts.Properties.VariableNames{1} = 'XRCC5';
counts{1,2} = height(NQO1_Qhits);
counts{2,2} = height(NQO1_Rhits);
counts.Properties.VariableNames{2} = 'NQO1';
counts{1,3} = height(SET_Qhits);
counts{2,3} = height(SET_Rhits);
counts.Properties.VariableNames{3} = 'SET';
counts{1,4} = height(S100A11_Qhits);
counts{2,4} = height(S100A11_Rhits);
counts.Properties.VariableNames{4} = 'S100A11';
counts.Properties.RowNames{1} = '# Query 2K Hits / 1835';
counts.Properties.RowNames{2} = '# Reference Hits / 1050';
display (counts)
writetable(counts,fullfile(results_folder,'counts_bwDMSO_NBT.xlsx'));

% calculate intersection between 4lines
percents.inter.XRCC5_NQO1 = height(cell2table(intersect(XRCC5_Qhits.drug_names,NQO1_Qhits.drug_names)))/height(compound_IDs)*100;
percents.inter.XRCC5_SET=height(cell2table(intersect(XRCC5_Qhits.drug_names,SET_Qhits.drug_names)))/height(compound_IDs)*100;
percents.inter.XRCC5_S100A11=height(cell2table(intersect(XRCC5_Qhits.drug_names,S100A11_Qhits.drug_names)))/height(compound_IDs)*100;
percents.inter.S100A11_SET=height(cell2table(intersect(SET_Qhits.drug_names,S100A11_Qhits.drug_names)))/height(compound_IDs)*100;
percents.inter.SET_NQO1=height(cell2table(intersect(SET_Qhits.drug_names,NQO1_Qhits.drug_names)))/height(compound_IDs)*100;
percents.inter.S100A11_NQO1=height(cell2table(intersect(NQO1_Qhits.drug_names,S100A11_Qhits.drug_names)))/height(compound_IDs)*100;
percents.inter.All4 = height(cell2table(intersect(intersect(XRCC5_Qhits.drug_names,NQO1_Qhits.drug_names),intersect(SET_Qhits.drug_names,S100A11_Qhits.drug_names))))/height(compound_IDs)*100;

% calculate union between 4lines
percents.union.XRCC5_NQO1 = height(cell2table(union(XRCC5_Qhits.drug_names,NQO1_Qhits.drug_names)))/height(compound_IDs)*100;
percents.union.XRCC5_S100A11 = height(cell2table(union(XRCC5_Qhits.drug_names,S100A11_Qhits.drug_names)))/height(compound_IDs)*100;
percents.union.XRCC5_SET = height(cell2table(union(XRCC5_Qhits.drug_names,SET_Qhits.drug_names)))/height(compound_IDs)*100;
percents.union.S100A11_SET = height(cell2table(union(S100A11_Qhits.drug_names,SET_Qhits.drug_names)))/height(compound_IDs)*100;
percents.union.S100A11_NQO1 = height(cell2table(union(S100A11_Qhits.drug_names,NQO1_Qhits.drug_names)))/height(compound_IDs)*100;
percents.union.SET_NQO1 = height(cell2table(union(SET_Qhits.drug_names,NQO1_Qhits.drug_names)))/height(compound_IDs)*100;
percents.union.All4=height(cell2table(union(union(S100A11_Qhits.drug_names,SET_Qhits.drug_names),union(XRCC5_Qhits.drug_names,NQO1_Qhits.drug_names))))/height(compound_IDs)*100;

percents = vertcat((struct2table(percents.union)),(struct2table(percents.inter)));
percents.Properties.RowNames{1}= 'Union';
percents.Properties.RowNames{2}='Intersect';
display(percents)
writetable(percents,fullfile(results_folder,'Q_bwDMSO_NBT_unions_intersects.xlsx'));

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
saveas(gca,fullfile(results_folder,'PVALS_QUERY'))
print('-f1',fullfile(results_folder,'PVALS_QUERY'),'-dpng','-r500');
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
xlim([1 4000])
ylim([10E-19 1])
saveas(gca,fullfile(results_folder,'PVALS_ALL'))
print('-f1',fullfile(results_folder,'PVALS_ALL'),'-dpng','-r500');
close all

%% Plot P Values per Plate
close all
% Select pORACL
figure;
pORACL = 'A549_NQO1';
cellplates = plates(plates.CellLine == pORACL,:);
for a = 1:height(cellplates)
    NQO1_bioactive.plateIDs = categorical(NQO1_bioactive.plateIDs);
    pVals = NQO1_bioactive(NQO1_bioactive.plateIDs == cellplates.expt_plate{a,1},:);
    pVals = sortrows(pVals,'NQO1_bioactive_pVals','ascend');
    plot(pVals{:,7})
    hold on
end
legend('Ref B2','Screen B2','Screen B2','Ref B2','Ref B4','Screen B4','Screen B4','Ref B4','Ref B6','Screen B6','Screen B6','Ref B6')
title('P Value by Plate NQO1')
xlabel('Compound')
ylabel('P Value')
hold off
saveas(gca,fullfile(results_folder,'PValues_PathwayTargetNQO1'));

figure;    
pORACL = 'A549_SET';
cellplates = plates(plates.CellLine == pORACL,:);
for a = 1:height(cellplates)
    SET_bioactive.plateIDs = categorical(SET_bioactive.plateIDs);
    pVals = SET_bioactive(SET_bioactive.plateIDs == cellplates.expt_plate{a,1},:);
    pVals = sortrows(pVals,'SET_bioactive_pVals','ascend');
    plot(pVals{:,7})
    hold on
end
legend('Ref B2','Screen B2','Screen B2','Ref B2','Ref B4','Screen B4','Screen B4','Ref B4','Ref B6','Screen B6','Screen B6','Ref B6')
title('P Value by Plate SET')
xlabel('Compound')
ylabel('P Value')
hold off
saveas(gca,fullfile(results_folder,'PValues_PathwayTargetSET'));

figure;
pORACL = 'A549_S100A11';
cellplates = plates(plates.CellLine == pORACL,:);
for a = 1:height(cellplates)
    S100A11_bioactive.plateIDs = categorical(S100A11_bioactive.plateIDs);
    pVals = S100A11_bioactive(S100A11_bioactive.plateIDs == cellplates.expt_plate{a,1},:);
    pVals = sortrows(pVals,'S100A11_bioactive_pVals','ascend');
    plot(pVals{:,7})
    hold on
end
legend('Ref B1','Screen B1','Screen B1','Ref B1','Ref B3','Screen B3','Screen B3','Ref B3','Ref B5','Screen B5','Screen B5','Ref B5')
title('P Value by Plate S100A11')
xlabel('Compound')
ylabel('P Value')
hold off
saveas(gca,fullfile(results_folder,'PValues_PathwayTargetS100A11'));

figure;
pORACL = 'A549_XRCC5';
cellplates = plates(plates.CellLine == pORACL,:);
for a = 1:height(cellplates)
    XRCC5_bioactive.plateIDs = categorical(XRCC5_bioactive.plateIDs);
    pVals = XRCC5_bioactive(XRCC5_bioactive.plateIDs == cellplates.expt_plate{a,1},:);
    pVals = sortrows(pVals,'XRCC5_bioactive_pVals','ascend');
    plot(pVals{:,7})
    hold on
end
legend('Ref B1','Screen B1','Screen B1','Ref B1','Ref B3','Screen B3','Screen B3','Ref B3','Ref B5','Screen B5','Screen B5','Ref B5')
title('P Value by Plate XRCC5')
xlabel('Compound')
ylabel('P Value')
hold off
saveas(gca,fullfile(results_folder,'PValues_PathwayTargetXRCC5'));

%% Plot Log P by Pathway etc.
close all
bioactives = sortrows(bioactives,'drug_names','ascend'); % reset for consistent figures
bioactives = sortrows(bioactives,{'Pathway',...
                                  'XRCC5_bioactive_pVals'}...
                                ,{'ascend',...
                                  'ascend'});
bioactives1 = bioactives{:,4:7};
bioactives1(isnan(bioactives1))=1; 
bioactives2 = log(bioactives1);
LE_heatmap_pvalues(bioactives2,'Sort by Pathway,XRCC5','Log P Values');
saveas(gca,fullfile(results_folder,'LogPValues_Pathway_XRCC5'));
print('-f1',fullfile(results_folder,'LogPValues_Pathway_XRCC5'),'-dpng','-r500');

close all
bioactives = sortrows(bioactives,'drug_names','ascend'); % reset for consistent figures
bioactives = sortrows(bioactives,{'Pathway',...
                                  'Target',...
                                  'XRCC5_bioactive_pVals'}...
                                ,{'ascend',...
                                  'ascend',...
                                  'ascend'});
bioactives1 = bioactives{:,4:7};
bioactives1(isnan(bioactives1))=1; 
bioactives3 = log(bioactives1);
LE_heatmap_pvalues(bioactives3,'Sort by Pathway,Target,XRCC5','Log P Values');
saveas(gca,fullfile(results_folder,'LogPValues_PathwayTargetXRCC5'));
print('-f1',fullfile(results_folder,'LogPValues_Pathway_X_N'),'-dpng','-r500');

close all
bioactives = sortrows(bioactives,'drug_names','ascend'); % reset for consistent figures
bioactives = sortrows(bioactives,{'XRCC5_bioactive_pVals',...
                                  'Pathway',...
                                  'NQO1_bioactive_pVals',...
                                  'SET_bioactive_pVals',...
                                  'S100A11_bioactive_pVals'}...
                                ,{'ascend',...
                                  'ascend',...
                                  'ascend',...
                                  'ascend',...
                                  'ascend'});
bioactives1 = bioactives{:,4:7};
bioactives1(isnan(bioactives1))=1;
bioactives4 = log(bioactives1);
LE_heatmap_pvalues(bioactives4,'Sort by XRCC5, Pathway, pORACLs','Log P Values');
saveas(gca,fullfile(results_folder,'LogPValues_XRCC5PathwaypORALCS'))
print('-f1',fullfile(results_folder,'LogPValues_XRCC5PathwaypORACLS'),'-dpng','-r500');

close all
bioactives = sortrows(bioactives,'drug_names','ascend'); % reset for consistent figures
bioactives = sortrows(bioactives,{'NQO1_bioactive_pVals',...
                                  'XRCC5_bioactive_pVals',...
                                  'Pathway'}...
                                ,{'ascend',...
                                  'ascend',...
                                  'ascend'});
bioactives1 = bioactives{:,4:7};
bioactives1(isnan(bioactives1))=1;
bioactives4 = log(bioactives1);
LE_heatmap_pvalues(bioactives4,'Sort by NQO1, Pathway, pORACLs','Log P Values');
saveas(gca,fullfile(results_folder,'LogPValues_NQO1Pathway'))
print('-f1',fullfile(results_folder,'LogPValues_NQO1Pathway'),'-dpng','-r500');

%% Call Out Different Pathways of Interest
close all
category = 'Neuronal Signaling';
bioactives.Pathway = categorical(bioactives.Pathway);
Neuronal = bioactives(bioactives.Pathway==category,:);
Neuronal = Neuronal((~(Neuronal.NQO1_bioactive_pVals>pVal_thr...
                    & Neuronal.SET_bioactive_pVals>pVal_thr...
                    & Neuronal.S100A11_bioactive_pVals>pVal_thr...
                    & Neuronal.XRCC5_bioactive_pVals>pVal_thr)),:);
Neuronal = sortrows(Neuronal,{'XRCC5_bioactive_pVals',...
                          'NQO1_bioactive_pVals',...
                          'SET_bioactive_pVals',...
                          'S100A11_bioactive_pVals'}...
                        ,{'ascend',...
                          'descend',...
                          'descend',...
                          'descend'});                   
Others1 = Neuronal{:,4:7};
Others1(isnan(Others1))=1;
Others2 = log(Others1);
LE_heatmap_pvalues(Others2,'Sort by XRCC5,NQO1,SET,S100A11','Log P Values');
saveas(gca,fullfile(results_folder,'LogPValues_BioactiveNeuronal'))
print('-f1',fullfile(results_folder,'LogPValues_BioactiveNeuronal'),'-dpng','-r500');

