function pc = Visualize_Plate( plate, varargin)
% Create PCA visualization using the plate struct.
%
% Visualize_Plate( plate, 'Name',Value)
%
% 'Name':
%     'use_mode'
%         - 'time_trace': Plot time traces.
%         - 'point': Plot points (Default).
%         - 'con_trace': Plot concentration traces.
% 
%     'dim': dimension of MDS plot (Default = 3).
% 
%     'time_to_use': Default is using all time points.
% 
%     'cate_to_show': A cell vector indicating categories to be shown. Default is showing all drug categories.
% 
%     'drug_to_show': A cell vector indicating drugs to be shown. Default is showing all drugs.
% 
%     'rm_pos_ctrl': true(Default)|false, show or hide positive control drugs.
% 
%     'legend_mode': 'on'(Default), or 'off'.
% 
%     'text_label': A cell string used for labelling. Default is empty.
% 
%     'color_label': A group indices used for the marker shapes. Default is drug categories.
% 
%     'shape_label': A group indices used for the marker shapes. Default is plateIDs.
% 
%     'feature_for_biplot': A cell of feature names. Default is empty.
% 
% 
% 
% Change Log:
%     2013.11.13: Eliminate 'text_mode', and add 'text_label'.
%     2013.11.15: Add plotting time trace.
%     2014.2.25: Add 'feature_for_biplot'.
%     2014.5.21: Add 'rm_pos_ctrl'.
%     2016.05.16: Adapted to inputParser and reorganized the code.
% 
% 
% Chien-Hsiang Hsu
% 
% Copyright: Altschuler and Wu laboratories
% Author: Chien-Hsiang Hsu at the Altschuler and Wu Laboratories
% For latest updates, check: <http://www.altschulerwulab.org/>.
%
% All rights reserved.


%% Commonly used constants
drug = plate.drug_names;
con = plate.concentrations;
ndrug = length(drug);
dc = strtrim(plate.drug_categories);
[~,GN] = grp2idx(dc);
plateIDs = plate.plateIDs;


%% Parse inputs
P = inputParser;

default_use_mode            = 'point';
default_dim                 = 3;
default_time_to_use         = 1:size(plate.profiles,2);

default_cate_to_show        = GN;
default_drug_to_show        = unique(drug);
default_rm_pos_ctrl         = true;

default_feature_for_biplot  = {};
default_text_label          = {};
default_color_label          = dc;
default_shape_label         = plateIDs;
default_legend_mode         = 'on';

P.addParameter('use_mode',default_use_mode,@(x) ismember(x,{'time_trace','point','con_trace'}));
P.addParameter('dim',default_dim,@isscalar);
P.addParameter('time_to_use',default_time_to_use,@isnumeric);

P.addParameter('cate_to_show',default_cate_to_show,@iscellstr);
P.addParameter('drug_to_show',default_drug_to_show,@iscellstr);
P.addParameter('rm_pos_ctrl',default_rm_pos_ctrl,@islogical);

P.addParameter('legend_mode',default_legend_mode,@ischar);
P.addParameter('text_label',default_text_label,@iscellstr);
P.addParameter('color_label',default_color_label,@(x) length(x)==size(plate.profiles,1));
P.addParameter('shape_label',default_shape_label,@(x) length(x)==size(plate.profiles,1));
P.addParameter('feature_for_biplot',default_feature_for_biplot,@iscellstr);

P.parse(varargin{:});

use_mode            = P.Results.use_mode;
dim                 = P.Results.dim;
time_to_use         = P.Results.time_to_use;

cate_to_show        = P.Results.cate_to_show;
% change to the correct order according to GN
[~,loc] = ismember(GN,cate_to_show);
cate_to_show = cate_to_show(loc(loc~=0));

drug_to_show        = P.Results.drug_to_show;
rm_pos_ctrl         = P.Results.rm_pos_ctrl;

legend_mode         = P.Results.legend_mode;
text_label          = P.Results.text_label;
color_label         = P.Results.color_label;
shape_label         = P.Results.shape_label;
feature_for_biplot  = P.Results.feature_for_biplot;

ntp = length(time_to_use);


%% Prepare data for MDS 
data = plate.profiles;
show_only = ismember(dc,cate_to_show);


switch use_mode
    case 'time_trace'
        dataMds = cat(1,data{:});
        
        dlbl = repmat(drug,ntp,1); % label for drugs        
        tlbl = kron(time_to_use',ones(ndrug,1)); % label for time
        clbl = repmat(con,ntp,1); % label for concentrations   
        dclbl = strcat(dlbl,'_',strtrim(cellstr(num2str(clbl)))); 
        
        color_label = repmat(color_label,ntp,1); % drug category
        shape_label = repmat(shape_label,ntp,1);
        text_label = repmat(text_label,ntp,1);
        
        show_only = repmat(show_only,ntp,1) &... % category_to_show
            (ismember(dlbl,drug_to_show) | ismember(dclbl,drug_to_show)); % drug_to_show
        if rm_pos_ctrl % remove positive control drugs
            show_only = show_only & ~(clbl==0 & ~ismember(color_label,{'DMSO','Unknown'}));
        end
        
        feature_names = plate.feaInfo.name(:);
        
        
    case 'point'
        show_only = ismember(dc,cate_to_show) &...
            (ismember(drug,drug_to_show) | ismember(plate.drug_names,drug_to_show));
        if rm_pos_ctrl % remove positive control drugs
            show_only = show_only & ~strcmp(dc,'Positive Control');
        end
        
        dataMds = cell2mat(data(:,time_to_use)); % input for MDS
        feature_names = cell(ntp,1);
        for t = 1:ntp
            feature_names{t} = strcat(plate.feaInfo.name(:),['_t' num2str(t)]);
        end
        feature_names = cat(1,feature_names{:});
        
    case 'con_trace'
        tmp = data(:,time_to_use);
        dataMds = cell2mat(tmp(:));
        
        % prepare labels for the expanded matrix
        dlbl = repmat(drug,ntp,1); % label for drugs        
        tlbl = kron(time_to_use',ones(ndrug,1)); % label for time
        dtlbl = strcat(dlbl,'_t',strtrim(cellstr(num2str(tlbl))));
        
        clbl = repmat(con,ntp,1); % label for concentrations
        show_only = repmat(show_only,ntp,1) &... % category_to_show
            (ismember(dlbl,drug_to_show) | ismember(dtlbl,drug_to_show)); % drug_to_show            
        color_label = repmat(color_label,ntp,1);
        shape_label = repmat(shape_label,ntp,1);
        text_label = repmat(text_label,ntp,1);
        
        dtslbl = strcat(dlbl,'_t',strtrim(cellstr(num2str(tlbl))),'_',shape_label); 
        
        feature_names = plate.feaInfo.name(:);

    otherwise
            disp('Wrong use_mode specidied!')
            return
end


%% PCA
disp('PCAing...')
[pc,~,~,~,explained] = pca(dataMds);
disp('Done');
disp(['Variation explained: ' num2str(sum(explained(1:dim))) '%'])
for i = 1:dim
    disp(['PC' num2str(i) ': ' num2str(explained(i)) '%'])
end

ytemp = dataMds*pc;
ytemp = ytemp(:,1:dim);


%% Plot 
switch use_mode
    case 'point'
        fh = figure;
        plot_with_color_label( ytemp, color_label, text_label, shape_label, show_only )
        labels = Construct_Labels(plate);
        Datacursor_Annotation( fh, ytemp, labels )
        
    case 'con_trace'
        fh = figure;
        plot_con_traces(ytemp,dtslbl,clbl,color_label,shape_label,show_only,text_label)
        labels = repmat(Construct_Labels(plate),length(time_to_use),1);
        labels.time_point = tlbl;
%         labels = labels(srt,:);
        Datacursor_Annotation( fh, ytemp, labels )

    case 'time_trace'
        fh = figure;
        plot_time_trace(ytemp,tlbl,color_label,shape_label,show_only,text_label)
        labels = repmat(Construct_Labels(plate),length(time_to_use),1);
        labels.time_point = tlbl;
        Datacursor_Annotation( fh, ytemp, labels )
end

if strcmpi(legend_mode,'off')
    legend('off')    
end


if ~isempty(feature_for_biplot)
    hold on
    
%     fea_to_show = ismember(feature_names,feature_for_biplot);
%     fea_to_show = ~cellfun(@isempty,regexp(feature_names,feature_for_biplot));
    fea_to_show = IsSubstr(feature_names,feature_for_biplot);
    
    fea_names = feature_names(fea_to_show);
    retained_var = diag(pc(fea_to_show,1:dim)*diag(lat(1:dim))*pc(fea_to_show,1:dim)'); % retained variances of each feature after projection.
    retained_var_percent = retained_var./diag(cov(dataMds(:,fea_to_show)));
    fea_dir = pc(fea_to_show,1:dim); % feature vectors with unequal length in the projected space.
    fea_dir = diag(sqrt(retained_var))*bsxfun(@rdivide,fea_dir,sqrt(diag(fea_dir*fea_dir'))); % let feature vectors have length proportioal to their std.
    scale_fac = max(abs(ytemp(:)))/max(sqrt(diag(fea_dir*fea_dir')));
    fea_dir = fea_dir*scale_fac;
    

    for f = 1:size(fea_dir)
        switch dim
            case 3
                plot3([0,fea_dir(f,1)],[0,fea_dir(f,2)],[0,fea_dir(f,3)],'k','linewidth',2)
                text(fea_dir(f,1),fea_dir(f,2),fea_dir(f,3),...
                    [fea_names{f},...
                    '(' num2str(retained_var(f),2) ', ',...
                    num2str(retained_var_percent(f),2) ')'],...
                    'interpreter','none','fontsize',10,'fontweight','bold')
            case 2
                plot([0,fea_dir(f,1)],[0,fea_dir(f,2)],'k','linewidth',2)
                text(fea_dir(f,1),fea_dir(f,2),...
                    [fea_names{f},...
                    '(' num2str(retained_var(f),2) ', ',...
                    num2str(retained_var_percent(f),2) ')'],...
                    'interpreter','none','fontsize',10,'fontweight','bold')
        end
        
    end
end



end


%------------------------%
function plot_con_traces(ytemp,dtslbl,clbl,color_label,shape_label,show_only,text_label)

dim = size(ytemp,2);
if ~ismember(dim,[2,3])
    error('Wrong dimension.')
end

% [G_c,GN_c] = grp2idx(clbl(clbl~=0));
% drug_num = size(ytemp(clbl~=0),1)/length(GN_c); % number of drugs in a time points. 

[G_drug,GN_drug] = grp2idx(dtslbl); % group of drug_time_plate
[G_color,GN_color] = grp2idx(color_label);

ls = {'.-','.--','.:','.-.','.-','.--','.:','.-.'}; % for line
ls2 = {'o','s','p','d','*','h','o','s','p','d','h','*','o','s','p','d','h','*'}; % for points
color_code = Drug_Category_Color(GN_color);

tmp = shape_label(show_only); % degenerate shape of those not being shown
shape_label(~show_only) = tmp(1);
[G_shape,GN_shape] = grp2idx(shape_label);


%%% Plot drugs without concentration %%%
plot_idx = clbl==0 & show_only;
switch dim
    case 3        
        scatter3(ytemp(plot_idx,1),ytemp(plot_idx,2),ytemp(plot_idx,3),10,...
            color_code(G_color(plot_idx),:),'fill')
       
    case 2
        scatter(ytemp(plot_idx,1),ytemp(plot_idx,2),10,...
            color_code(G_color(plot_idx),:),'fill')

end
hold on


%%% Plot concentration traces %%%
for d = 1:length(GN_drug)
    plot_idx = (G_drug==d) & (clbl~=0) & (show_only);
    
    if sum(plot_idx)==0
        continue
    end
    
    plot_data = ytemp(plot_idx,:);
    if ~isempty(text_label)
        plot_text = text_label(plot_idx);
    else
        plot_text = {};
    end
    
    % sort by concentration
    [~,srt] = sort(clbl(plot_idx,:),'ascend');
    plot_data = plot_data(srt,:);
    
    d_first = find(plot_idx,1); % just need the first instance to set color, shape...
    
    switch color_label{d_first}
        case {'DMSO','DMSO + TGF'}
%             opts = {ls{G_shape(d_first)},'color',color_code(G_color(d_first),:),'linewidth',1};
            opts = {'o','markerfacecolor',color_code(G_color(d_first),:),...
                'markeredgecolor',color_code(G_color(d_first),:),'markersize',4};
            opts2 = {'o','markerfacecolor',color_code(G_color(d_first),:),...
                'markeredgecolor',color_code(G_color(d_first),:),'markersize',4};
        
        case 'Unknown'
            opts = {ls{G_shape(d_first)},'color',color_code(G_color(d_first),:),'linewidth',1,...
                'markerfacecolor',color_code(G_color(d),:),'markersize',20};
            opts2 = {ls2{G_shape(d_first)},'markerfacecolor',color_code(G_color(d_first),:),...
                'markeredgecolor','k','markersize',7};
        
        case 'Nonbioactive'
            opts = {ls{G_shape(d_first)},'color',[0.8,0.8,0.8],'linewidth',1};
            opts2 = {'o','markerfacecolor',[0.8,0.8,0.8],...
                'markeredgecolor',[0.8,0.8,0.8],'markersize',2};
            
        case 'DMSO_no_TGF'
            opts = {'o','markerfacecolor',color_code(G_color(d_first),:),...
                'markeredgecolor',color_code(G_color(d_first),:),'markersize',4};
            opts2 = {'o','markerfacecolor',color_code(G_color(d_first),:),...
                'markeredgecolor',color_code(G_color(d_first),:),'markersize',4};
        
        otherwise
            opts = {ls{G_shape(d_first)},'color',color_code(G_color(d_first),:),'linewidth',1,...
                'markerfacecolor',color_code(G_color(d_first),:),'markersize',12};
            opts2 = {ls2{G_shape(d_first)},'markerfacecolor',color_code(G_color(d_first),:),...
                'markeredgecolor','k','markersize',7};
    end
    
    % Plot the end point (highest concentration) again
    switch dim
        case 3
            h(G_shape(d_first),G_color(d_first)) = plot3(plot_data(:,1),plot_data(:,2),plot_data(:,3),opts{:});
            hold on
            plot3(plot_data(1,1),plot_data(1,2),plot_data(1,3),opts2{:});
            if ~isempty(plot_text)
                text(plot_data(1,1),plot_data(1,2),plot_data(1,3),plot_text(1),...
                     'interpreter','none')
            end

        case 2
            h(G_shape(d_first),G_color(d_first)) = plot(plot_data(:,1),plot_data(:,2),opts{:});
            plot2(plot_data(1,1),plot_data(1,2),opts2{:});
            if ~isempty(plot_text)
                text(plot_data(1,1),plot_data(1,2),plot_text(1),'interpreter','none')
            end
    end
end

rg = max(abs(ytemp(:)));
%xlim([-rg,rg])
%ylim([-rg,rg])
%zlim([-rg,rg])
grid on

if length(GN_shape)>1
    for r = 1:length(GN_shape)
        switch dim
            case 3
                f(r) = plot3(0,0,0,ls{r},'color',[0 0 0],'linewidth',1,'markerfacecolor',[0 0 0],'markersize',9);
            case 2
                f(r) = plot(0,0,ls{r},'color',[0 0 0],'linewidth',1,'markerfacecolor',[0 0 0],'markersize',9);
        end        
        set(f(r),'visible','off')   
    end
    legend([h(1,ishghandle(h(1,:))),f],{GN_color{any(ishghandle(h))},GN_shape{:}},'interpreter','none');
else
    legend(h(1,ishghandle(h(1,:))),GN_color(ishghandle(h(1,:))),'interpreter','none');
end    

end



function plot_time_trace(ytemp,tlbl,color_label,shape_label,show_only,text_label)

dim = size(ytemp,2);
if ~ismember(dim,[2,3])
    error('Wrong dimension.')
end
ls = {'.-','.--','.:','.-.','.-','.--','.:','.-.'}; % for line
ls2 = {'o','s','p','d','*','h','o','s','p','d','h','*','o','s','p','d','h','*'}; % for points

[G_t,GN_t] = grp2idx(tlbl);
drug_num = size(ytemp,1)/length(GN_t); % number of drugs in a time points. 

[G_color,GN_color] = grp2idx(color_label);
color_code = Drug_Category_Color(GN_color); 

tmp = shape_label(show_only); % degenerate shape of those not being shown
shape_label(~show_only) = tmp(1);
[G_shape,GN_shape] = grp2idx(shape_label);

dd = max(abs(ytemp(:)))*0.01; % small distance between text and data point


for d = 1:drug_num
    if ~show_only(d)
        continue
    end
    
    if strcmpi(color_label{d},'DMSO')
        opts = {ls{G_shape(d)},'color',color_code(G_color(d),:),'linewidth',1};
        opts2 = {'o','markerfacecolor',color_code(G_color(d),:),...
            'markeredgecolor',color_code(G_color(d),:),'markersize',2};
        
    elseif strcmpi(color_label{d},'Unknown')
        opts = {ls{G_shape(d)},'color',color_code(G_color(d),:),'linewidth',1,...
            'markerfacecolor',color_code(G_color(d),:),'markersize',12};
        opts2 = {ls2{G_shape(d)},'markerfacecolor',color_code(G_color(d),:),...
            'markeredgecolor','k','markersize',7};
        
    else
        opts = {ls{G_shape(d)},'color',color_code(G_color(d),:),'linewidth',1,...
            'markerfacecolor',color_code(G_color(d),:),'markersize',12};
        opts2 = {ls2{G_shape(d)},'markerfacecolor',color_code(G_color(d),:),...
            'markeredgecolor','k','markersize',7};
    end
    
    plot_data = ytemp(d:drug_num:end,:);
    switch dim
        case 2
            h(G_shape(d),G_color(d)) = plot(plot_data(:,1),plot_data(:,2),opts{:});
            if ~isempty(text_label)
                text(plot_data(end,1)+dd,plot_data(end,2)+dd,text_label{d},...
                    'color',[0,0,0],'fontweight','bold','interpreter','none');
            end            
            hold on
        case 3
            h(G_shape(d),G_color(d)) = plot3(plot_data(:,1),plot_data(:,2),plot_data(:,3),opts{:});
            if ~isempty(text_label)
                text(plot_data(end,1)+dd,plot_data(end,2)+dd,plot_data(end,3)+dd,text_label{d},...
                    'color',[0,0,0],'fontweight','bold','interpreter','none');
            end            
            hold on
            plot3(plot_data(end,1),plot_data(end,2),plot_data(end,3),opts2{:});
        otherwise
            disp('Wrong dimension.')
    end
end

rg = max(abs(ytemp(:)));
%xlim([-rg,rg])
%ylim([-rg,rg])
%zlim([-rg,rg])
grid on

if length(GN_shape)>1
    for r = 1:length(GN_shape)
        switch dim
            case 3
                f(r) = plot3(0,0,0,ls{r},'color',[0,0,0],'linewidth',1,'markerfacecolor',[0,0,0],'markersize',9);
            case 2
                f(r) = plot(0,0,ls{r},'color',[0,0,0],'linewidth',1,'markerfacecolor',[0,0,0],'markersize',9);
        end        
        set(f(r),'visible','off')       
    end
    legend([h(1,ishghandle(h(1,:))),f],{GN_color{any(ishghandle(h))},GN_shape{:}},'interpreter','none')
else
    legend(h(1,ishghandle(h(1,:))),GN_color(ishghandle(h(1,:))),'interpreter','none')
end

end


function labels = Construct_Labels(plate)
label_struct = rmfield(plate,{'profiles','feaInfo'}); % ,'nCells'
labels = struct2table(label_struct);
end


function is_substr = IsSubstr(str_to_check,str_pattern)
is_substr = false(length(str_to_check),length(str_pattern));

for s = 1:length(str_pattern)
    is_substr(:,s) = ~cellfun(@isempty,regexp(str_to_check,str_pattern{s}));
end
is_substr = any(is_substr,2);
end












