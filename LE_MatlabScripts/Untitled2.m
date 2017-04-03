pathways.Pathway = categorical(pathways.Pathway);
pathwaysize = tabulate(pathways.Pathway)

for i = 1:length(cluster_info)
  query_cluster(i).categories = pathwaysize(:,1)
  query_cluster(i).cluster_ID = unique(cluster_info{i}.Cluster);
  query_cluster(i).stat = zeros(length(query_cluster(i).categories),length(query_cluster(i).cluster_ID));
  for j = 1:size(query_cluster(i).stat,1)
        is_class = strcmpi(cluster_info{i}.drug_categories,query_cluster(i).categories{j});
        tb = tabulate(cluster_info{i}.Cluster(is_class));
        
        [is_in,loc] = ismember(tb(:,1),query_cluster(i).cluster_ID);
        query_cluster(i).stat(j,loc(is_in)) = tb(:,2);
  end
      