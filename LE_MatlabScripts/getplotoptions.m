%getplotoptions

[targets,pathways, compound_IDs]=getselleckpathways();
uniquepathways_S = table2cell(unique(pathways));
uniquetargets_S = table2cell(unique(targets));
plot_opt_query = {'cate_to_show',horzcat(uniquepathways_S',{'DMSO'})};
plot_opt_reference = {'cate_to_show',{'DMSO','Actin','AuroraB','DNA','ER','HDAC','HSP90','MT','PLK','Proteasome','mTOR'}}; %plot_opt_reference = {'cate_to_show',{'DMSO','Proteasome','mTOR'}};
plot_opt_ALL = {'cate_to_show',horzcat(plot_opt_query{1,2},plot_opt_reference{1,2})};
