
ref_cluster.categories2 = cell2table(ref_cluster.categories)
ref_cluster.counts2 = cell2table(tabulate(cluster_info{1}.drug_categories))
ref_cluster.drugcounts = join(ref_cluster.categories2,ref_cluster.counts2)
ref_cluster.clustercounts=tabulate(cluster_info{1}.Cluster)



drugcounts = table2array(ref_cluster.drugcounts(:,2))
drugcounts = repmat(drugcounts,1,139)
x = ref_cluster.stat./drugcounts

clustercounts = ref_cluster.clustercounts(:,2)';
clustercounts = repmat(clustercounts,32,1);
y = ref_cluster.stat./clustercounts



x_ax = x(:)
y_ax = y(:)

scatter(x_ax,y_ax)
i =1;
fh = figure;
ax(i) = subplot(1,1,i,'Parent',fh)
imagesc(x)
colormap(ax(i),'hot')
colorbar(ax(i))
ax(i).YTick = 1:length(ref_cluster(i).categories);
ax(i).YTickLabel = ref_cluster(i).categories;
ax(i).TickLabelInterpreter = 'none';

xlabel('Cluster #')
ylabel('Compound Category (Pathway)')
title('Percentage Drugs in a Compound Category Appearing in each Cluster')    
    
    
    close all
    