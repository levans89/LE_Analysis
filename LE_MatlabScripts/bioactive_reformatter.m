%% master table
compound_IDs.Properties.VariableNames{1} = 'drug_names';
targets.Properties.VariableNames{1} = 'targets';
pathways.Properties.VariableNames{1} = 'pathways';
compound_IDs = horzcat(compound_IDs, targets, pathways);

%%import bioactive data


SETbioactive.cpd_usage = categorical(SETbioactive.cpd_usage);
NQO1bioactive.cpd_usage = categorical(NQO1bioactive.cpd_usage);
S100A11bioactive.cpd_usage = categorical(S100A11bioactive.cpd_usage);
XRCC5bioactive.cpd_usage = categorical(XRCC5bioactive.cpd_usage);

XRCC5bioactive.Properties.VariableNames{7} = 'XRCC_bioactive_pVals';
SETbioactive.Properties.VariableNames{7} = 'SET_bioactive_pVals';
S100A11bioactive.Properties.VariableNames{7} = 'S100A11_bioactive_pVals';
NQO1bioactive.Properties.VariableNames{7} = 'NQO1_bioactive_pVals';

SET_b_query = SETbioactive(SETbioactive.cpd_usage =='query_cpd',:);
SET_b_reference = SETbioactive(SETbioactive.cpd_usage =='reference_cpd',:);
XRCC5_b_query = XRCC5bioactive(XRCC5bioactive.cpd_usage =='query_cpd',:);
XRCC5_b_reference = XRCC5bioactive(XRCC5bioactive.cpd_usage =='reference_cpd',:);
S100A11_b_query = S100A11bioactive(S100A11bioactive.cpd_usage =='query_cpd',:);
S100A11_b_reference = S100A11bioactive(S100A11bioactive.cpd_usage =='reference_cpd',:);
NQO1_b_query = NQO1bioactive(NQO1bioactive.cpd_usage =='query_cpd',:);
NQO1_b_reference = NQO1bioactive(NQO1bioactive.cpd_usage =='reference_cpd',:);

XRCC5_b_query(:,'cpd_usage') = [];
SET_b_query(:,'cpd_usage') = [];
S100A11_b_query(:,'cpd_usage') = [];
NQO1_b_query(:,'cpd_usage') = [];


% reshaping and get all into one table
bioactives = outerjoin(compound_IDs, SET_b_query, 'Type','Left');
bioactives2 = outerjoin(compound_IDs, S100A11_b_query,'Type','Left');
bioactives3 = outerjoin(compound_IDs, NQO1_b_query,'Type','Left');
bioactives4 = outerjoin(compound_IDs, XRCC5_b_query,'Type','Left');
bioactives5 = outerjoin(bioactives,bioactives2,'Type','Left');
bioactives6 = outerjoin(bioactives3,bioactives4,'Type','Left');
bioactives6.Properties.VariableNames{1} = 'drug_names_compound_IDs_bioactives';
bioactives7 = outerjoin(bioactives5,bioactives6,'Type','Left');
bioactives = bioactives7;
bioactives.Properties.VariableNames{1} = 'drug_names';
bioactives(:,'drug_names_SET_b_query') = [];
bioactives(:,'drug_names_compound_IDs_bioactives2') = [];
bioactives(:,'drug_names_S100A11_b_query') = [];
bioactives(:,{'drug_names_compound_IDs_bioactives_bioactives6','drug_names_NQO1_b_query'}) = [];
bioactives(:,{'drug_names_compound_IDs_bioactives4','drug_names_XRCC5_b_query'}) = [];
clear bioactives2 bioactives3 bioactives4 bioactives5 bioactives6 bioactives7
bioactives(:,{'targets_bioactives2','pathways_bioactives2'}) = [];
bioactives(:,{'targets_bioactives3','pathways_bioactives3'}) = [];
bioactives(:,{'targets_bioactives4','pathways_bioactives4'}) = [];
