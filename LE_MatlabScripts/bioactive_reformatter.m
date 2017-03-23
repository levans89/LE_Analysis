%% master table
compound_IDs.Properties.VariableNames{1} = 'drug_names';
targets.Properties.VariableNames{1} = 'targets';
pathways.Properties.VariableNames{1} = 'pathways';
compound_IDs = horzcat(compound_IDs, targets, pathways);

SET_bioactive.cpd_usage = categorical(SET_bioactive.cpd_usage);
NQO1_bioactive.cpd_usage = categorical(NQO1_bioactive.cpd_usage);
S100A11_bioactive.cpd_usage = categorical(S100A11_bioactive.cpd_usage);
XRCC5_bioactive.cpd_usage = categorical(XRCC5_bioactive.cpd_usage);

XRCC5_bioactive.Properties.VariableNames{7} = 'XRCC_bioactive_pVals';
SET_bioactive.Properties.VariableNames{7} = 'SET_bioactive_pVals';
S100A11_bioactive.Properties.VariableNames{7} = 'S100A11_bioactive_pVals';
NQO1_bioactive.Properties.VariableNames{7} = 'NQO1_bioactive_pVals';

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

% reshaping and get all into one table
bioactives = outerjoin(compound_IDs, SET_b_query(:,[4,7]), 'Type','Left');
bioactives.Properties.VariableNames{1} = 'drug_names';
bioactives = outerjoin(bioactives, S100A11_b_query(:,[4,7]),'Type','Left');
bioactives.Properties.VariableNames{1} = 'drug_names';
bioactives = outerjoin(bioactives, NQO1_b_query(:,[4,7]),'Type','Left');
bioactives.Properties.VariableNames{1} = 'drug_names';
bioactives = outerjoin(bioactives, XRCC5_b_query(:,[4,7]),'Type','Left');
bioactives.Properties.VariableNames{1} = 'drug_names';
bioactives(:,{'drug_names_right','drug_names_right_1','drug_names_right_2','drug_names_right_3'}) = [];

% Plot P Values
bioactives = sortrows(bioactives,'S100A11_bioactive_pVals','ascend');
plot(bioactives.S100A11_bioactive_pVals)
hold on
bioactives = sortrows(bioactives,'NQO1_bioactive_pVals','ascend');
plot(bioactives.NQO1_bioactive_pVals)
bioactives = sortrows(bioactives,'XRCC_bioactive_pVals','ascend');
plot(bioactives.XRCC_bioactive_pVals)
bioactives = sortrows(bioactives,'SET_bioactive_pVals','ascend');
plot(bioactives.SET_bioactive_pVals)

% Plot Overlap
