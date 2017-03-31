function Show_PVALUE( plate2 )
% Generate heat maps of PVALUE.
% % Show_PVALUE( plate2 )
% % Louise Heinrich 2017.03.29, based on Show_Cell_Number by Chien-Hsiang Hsu, 2016.04.26

%% Commonly used constants
plate_IDs = unique(plate2.plateIDs);
Nplates = length(plate_IDs);
Nt = size(plate2.profiles,2);
max_PVAL = 0 % max(plate2.nCells(:));
Rows = cellstr(('A':'P')');

%% Loop through each plate and plot cell numbers
figure,
sp_idx = 1; % index of subplot

for p = 1:Nplates
    % 1. Get the data associated with this plate
    p_tmp = Select_Rows(plate2,strcmp(plate2.plateIDs,plate_IDs{p}));
    
    % 2. Get the mapping to the physical plate format
    p_physical = Plate2PhysicalPlate(p_tmp);
    p_physical.logbioactive_pVals(isnan(p_physical.logbioactive_pVals))=1;
    % 3. Plot cell number
    for t = 1:Nt
        ax(sp_idx) = subplot(Nplates,Nt,sp_idx);
        imagesc(p_physical.logbioactive_pVals(:,:,t))
        colormap(hot),colorbar
        caxis([-37,max_PVAL])
        title([plate_IDs{p}, ' Log P Value'],'interpreter','none')
        
        ax(sp_idx).YTick = 1:1:size(p_physical.logbioactive_pVals,1);
        ax(sp_idx).YTickLabel = Rows(ax(sp_idx).YTick);
        
        sp_idx = sp_idx + 1;
    end
end

end





