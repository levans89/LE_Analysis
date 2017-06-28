%% Set P Val Thr
pVal_thr            = 10^-6;

%% set categorical
bioactive.cpd_usage = categorical(bioactive.cpd_usage);
clc
query = bioactive(bioactive.cpd_usage =='query_cpd',:);
reference =bioactive(bioactive.cpd_usage =='reference_cpd',:);
remainders = bioactive(bioactive.cpd_usage ~='reference_cpd' & bioactive.cpd_usage ~='query_cpd',:);

%% Calculate Bioactive Percentages
percents.bioactive = sum(query.bioactive_pVals <= pVal_thr)/height(compound_IDs)*100; 
percents.ref = sum(reference.bioactive_pVals <= pVal_thr)/height(reference)*100;

percents_summary = vertcat(struct2table(percents.ref),struct2table(percents.bioactive));
percents_summary.Properties.RowNames{1}='Reference';
percents_summary.Properties.RowNames{2}='ChemBridgeV2';

display(percents_summary)
writetable(percents_summary,fullfile(results_folder,'percents.xlsx'))

%% calculate intersects and unions of query hits
% get hit compounds
Qhits = query(query.bioactive_pVals<=pVal_thr,:);

% get hit reference compounds
Rhits = reference(reference.bioactive_pVals<=pVal_thr,:);

% calculate bioactive numbers
counts = table;
counts{1,1} = height(Qhits);
counts{2,1} = height(Rhits);
counts.Properties.VariableNames{1} = 'XRCC5';
counts.Properties.RowNames{1} = '# ChemBridge Hits';
counts.Properties.RowNames{2} = '# Reference Hits';
display (counts)
writetable(counts,fullfile(results_folder,'counts.xlsx'));
%% Plot P Values for all Compounds across the dataset
figure()
hold on
bioactive = sortrows(bioactive,'bioactive_pVals','descend');
plot(bioactive.bioactive_pVals,'LineWidth',1)
xlabel('P Values All Bioactive');
ylabel('P Value');
set(gca,'YScale','log')
saveas(gca,fullfile(results_folder,'PVALS_ALL'))
print('-f1',fullfile(results_folder,'PVALS_ALL'),'-dpng','-r500');
close all

%% Plot P Values per Plate
plateIDs = bigscreen
batched = horzcat(plateIDs.expt_plate,plateIDs.batch,plateIDs.plate_type)
display(batched)
for a = 1:height(plateIDs)
    bioactive.plateIDs = categorical(bioactive.plateIDs);
    pVals = bioactive(bioactive.plateIDs == plateIDs.expt_plate{a,1},:);
    pVals = sortrows(pVals,'bioactive_pVals','ascend');
    plot(pVals{:,7})
    hold on
end
title('P Value by Plate')
xlabel('Compound')
ylabel('P Value')
hold off
    
