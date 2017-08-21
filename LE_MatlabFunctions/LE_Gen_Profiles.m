function [] = LE_Gen_Profiles(plates,cellline,profiles_folder)

paras = Set_Paras(); % set parameters (path,bioactive,screen,cluster)
theseplates = ismember(plates.CellLine,cellline);
plates = plates(theseplates,:);
plate_maps_folder = paras.path.plate_maps; % where to find plate maps
plate_maps = cell(height(plates),1); % preallocate

for m = 1:height(plates) % for all plates in this cell line
    plate_maps{m} = [num2str(plates.expt_plate{m}),'.xlsx']; % find the platemap name
end

for p = 2:height(plates) % for all plates
    
    % 1. Get all inputs for Get_KS_Profiles
    plateID = plates.expt_plate{p};% ID
    plate_type = plates.plate_type{p}; % e.g. reference_plate
    plate_map_file = [plate_maps_folder filesep plate_maps{p}]; % expt plate file
    plate_layout = Read_Plate_Annotation(plate_map_file); % read in expt plate file
    
    % 2. Gather information for profiling
    well_names = {plate_layout.wellData.wellName}';
    drug_names = {plate_layout.wellData.Compound_ID}';
    drug_categories = {plate_layout.wellData.Compound_Category}';
    clone_names = {plate_layout.wellData.Cell_Line}';
    cpd_usage = {plate_layout.wellData.Compound_Usage}';
    dose = {plate_layout.wellData.Dose}';
    cas = {plate_layout.wellData.CAS_Number}';
    targets = {plate_layout.wellData.Target}';
    pathways = {plate_layout.wellData.Pathway}';
    descriptions = {plate_layout.wellData.Description}';
    concentrations = str2double({plate_layout.wellData.Dose_Category})';
    concentrations(isnan(concentrations))=0;
    
    % 3. Define KS control wells
    is_ks_ctrl = ismember(well_names(:),{'A2','A23',...
                                         'B2','B23',...
                                         'C2','C23',...
                                         'D2','D23',...
                                         'E2','E23',...
                                         'F2','F23',...
                                         'G2','G23',...
                                         'H2','H23',...
                                         'I2','I23',...
                                         'J2','J23',...
                                         'K2','K23',...
                                         'L2','L23',...
                                         'M2','M23',...
                                         'N2','N23',...
                                         'O2','O23',...
                                         'P2','P23'});
    
    % 4. Calculate profiles
    Get_KS_Profiles(plateID, plate_type, well_names, drug_names, drug_categories, concentrations,...
        dose, clone_names, cas, targets, pathways, descriptions, ...
        is_ks_ctrl, cpd_usage, paras,...
        'save_destination',profiles_folder);
    
end
end

