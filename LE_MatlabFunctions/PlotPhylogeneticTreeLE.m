function [H,phyTree] = PlotPhylogeneticTreeLE(cluster_info,tree,clusters_to_color)
% Visualize clusters with phylogenetic tree.
% 
% [H,phyTree] = PlotPhylogeneticTree(cluster_info,tree,clusters_to_color)
% 
% Inputs
%     - cluster_info: output of Cluster_Compounds.
%     - tree: corresponding tree. Make sure that node IDs in tree correspond to row IDs of cluster_info.
%     - clusters_to_color: cluster IDs to color, usually are significant clusters.
% 
% Chien-Hsiang Hsu, 2016.04.12


%% Set up 
nNodes = size(cluster_info,1);
drugClass = cluster_info.drug_categories; 
refFlag = ismember(cluster_info.cpd_usage,{'reference_cpd','negative_ctrl'});
unkown_flag = ismember(cluster_info.cpd_usage,{'query_cpd'});
% refFlag = ~ismember(cluster_info.cpd_usage,{'Unkown'});
% drugClass(refFlag) = cluster_info.drug_categories(refFlag);
drugClass(unkown_flag) = {'Unknown'};
[G,GN] = grp2idx(drugClass);
colorCode = Drug_Category_ColorLE(GN);


%% Get the tree of dendrogram
phyTree = phytree(tree);
H = plot(phyTree,'Type','radial'); % 'radial' | 'equaldaylight'


%% Get xyData to adjust the appearance of the phylogenetic tree
set(H.BranchDots,'Marker','none')

% 1. Get leaf node XY coordinates
xyData = [H.LeafDots.XData',H.LeafDots.YData'];

% 2. Get leaf nodeID (row number in data matrix)
nodeID = str2num(char(cellfun(@(x) x{2},regexpi(get(H.leafNodeLabels,'String'),' ','split'),...
    'UniformOutput',false)));

% 3. Sort xyData by nodeID (then it will be the same order as cluster_info)
[xyData,idx] = sortrows([xyData,nodeID],3);
H.BranchLines(1:nNodes) = H.BranchLines(idx); % sort H.BranchLines(leaf) according


%% Color leaf nodes by colorCode
hold on

% plot reference drgus first
for g = 1:length(GN)
    if ~strcmpi(GN{g},'Unknown')
       plot(xyData(G==g,1),xyData(G==g,2),'o','MarkerFaceColor',colorCode(g,:),...
        'MarkerEdgeColor',colorCode(g,:),'MarkerSize',5,'Linewidth',0.5)       
    end 
end

g = find(strcmpi(GN,'Unknown'));
if ~isempty(g)
    plot(xyData(G==g,1),xyData(G==g,2),'o','MarkerFaceColor',colorCode(g,:),...
            'MarkerEdgeColor',colorCode(g,:),'MarkerSize',3,'Linewidth',0.5)
end


%% Color edges by clusters 
set(H.BranchLines,'LineStyle',':','LineWidth',0.5,'Color',0.3*[1 1 1])
fullCluster = TreeExpanser(tree);

clusterOfInterest = []; % [DHFR, Na+/K+, Glucocoid]
% clusterOfInterest = [27,45,11]; % [DHFR, Na+/K+, Glucocoid]

rng(20150127)
for c = 1:numel(clusters_to_color)
    % Identify member leaf nodes
    memberNodes = find(cluster_info.Cluster==clusters_to_color(c));
    
    % Identify all branches
    branchNodes = find(cellfun(@(x) all(ismember(x,memberNodes)),fullCluster)) + nNodes;
    
    % Set to the right color
    if sum(refFlag(memberNodes)) > 1 % reference cluster
        drugLabel = sortrows(tabulate(drugClass(memberNodes)),-2);
        drugLabel(strcmpi(drugLabel(:,1),'Unknown'),:) = [];
        lineColor = colorCode(strcmpi(GN,drugLabel{1}),:); % use the most frequent label
        lineWidth = 1;
%         clusterLabel = drugLabel{1};
    else
%         lineColor = colorCode(strcmpi(GN,'Unknown'),:); % use the most frequent label
        lineColor = rand(1,3);
        lineWidth = 0.5;
%         clusterLabel = num2str(clusters_to_color(c));
    end
    
    if ismember(clusters_to_color(c),clusterOfInterest)
        switch clusters_to_color(c)
            case 27 % DHFR                
                lineColor = [1,0.4,0];
                
            case 45 % Na+/K+
                lineColor = [0.5696,0.5641,0.1099];
                
            case 11 % Glucocoid 
                lineColor = [0.3199,0.1777,0.6859];
        end
    end
    set(H.BranchLines([memberNodes;branchNodes]),'LineStyle','-','LineWidth',lineWidth,'Color',lineColor)
    
    % Label cluster
%     text(xyData(memberNodes(1),1),xyData(memberNodes(1),2),clusterLabel,...
%         'FontSize',12, 'Color', [1 1 1],'BackgroundColor', [0 0 0],...
%         'HorizontalAlignment', 'Center');
end


%% Associate datacursor with cluster info
Datacursor_Annotation( H.axes.Parent, xyData(:,1:2), cluster_info ); % xyData(:,3) is node ID


%% Final adjustment
set(H.terminalNodeLabels,'String','')
box off
axis(H.axes,'off')
H.axes.XTick = [];
H.axes.YTick = [];

H.axes.Parent.Color = [1 1 1];
H.axes.Parent.Units = 'Normalized';
H.axes.Parent.Position = [0.25,0.13,0.5,0.75];

end


