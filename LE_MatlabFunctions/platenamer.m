% Louise Heinrich in Altschuler & Wu Labs
% September 8 2017
% Lookup full platename from experimental db using barcode (for use in
% plate printing)

function [input_table]= platenamer(input_table, plates_folder,results_folder)
platelist = dir(plates_folder);
platelist = table({platelist.name}.', 'VariableNames', {'name'});
for i = 1:height(input_table)
    barcode = cellstr(table2cell(input_table(i,'plateIDs')));
    for s = 1:height(platelist)
        pos = strfind(platelist.name{s},barcode);
        if pos>0;
            %break
            %display(s)
            input_table.Var3{i,1}=platelist.name{s};
        else
        end
    end
end
input_table.Properties.VariableNames{3} = 'plateNames';
cd(results_folder)
writetable(input_table, 't.csv');
end
