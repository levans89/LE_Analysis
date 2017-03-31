function LE_heatmap_pvalues(cdata1,y,t)
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'YTick',[1 600 1200 1835],...
    'XTickLabel',{'SET','S100A11','NQO1','XRCC5'},...
    'XTick',[1 2 3 4],...
    'Layer','top',...
    'YDir','reverse');
box(axes1,'on');
hold(axes1,'on');

% Create image
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
colormap('gray')
% Create xlabel
xlabel('pORACL');

% Create ylabel
ylabel(y);

% Create title
title(t);

% Create colorbar
colorbar('peer',axes1);
ylim([0 1835])

