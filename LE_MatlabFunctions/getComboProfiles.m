% Leanna Morinishi 2017.4.01
% Concatenate_Profiles from NBT set across combinations of pORACL lines. 

function catStruct = getComboProfiles(line1, line2, line3, line4)

switch nargin
    case 1
        allProfiles = line1.profiles;
        allDrugs = line1.drug_names;
        allDose = line1.dose;
        allCats = line1.drug_categories;
        allPlates = line1.plate_types;
        allTargets = line1.targets;
        allPathways = line1.pathways;
        j = size(line1.profiles,1)+1;
        
    case 2
        line1_drugInfo = table(line1.drug_names,...
            line1.dose, line1.drug_categories,...
            line1.plate_types, line1.targets,...
            line1.pathways);
        line1_un = unique(line1_drugInfo, 'rows');

        line2_drugInfo = table(line2.drug_names,...
            line2.dose, line2.drug_categories,...
            line2.plate_types, line2.targets,...
            line2.pathways);
        line2_un = unique(line2_drugInfo, 'rows');
        
        a = getRows(line1_un, line2_un);
        drugs_in_common = line1_un(a,:);
        allProfiles = zeros(size(line1.profiles{1},1), size(line1.profiles{1},2)*2);
        allDrugs = repmat({''}, size(line1.profiles,1),1); 
        allDose = repmat({''}, size(line1.profiles,1),1);
        allCats = repmat({''}, size(line1.profiles,1),1);
        allPlates = repmat({''}, size(line1.profiles,1),1);
        allTargets = repmat({''}, size(line1.profiles,1),1);
        allPathways = repmat({''}, size(line1.profiles,1),1);

        j=1;
        for i = 1:size(drugs_in_common)
            drugi = drugs_in_common{i,1};
            dosei = drugs_in_common{i,2};
            cati = drugs_in_common{i,3};
            line1_rows = (strcmp(line1.drug_names, drugi) & strcmp(line1.dose, dosei));
            line2_rows = (strcmp(line2.drug_names, drugi) & strcmp(line2.dose, dosei));
            mini = min([sum(line1_rows), sum(line2_rows)]);

            line1_profs = cell2mat(line1.profiles(line1_rows,:));
            line2_profs = cell2mat(line2.profiles(line2_rows,:));

            allProfiles(j:(j+mini-1),:) = [line1_profs(1:mini,:) line2_profs(1:mini,:)];
            allDrugs(j:(j+mini-1)) = drugi(1);
            allDose(j:(j+mini-1)) = dosei(1);
            allCats(j:(j+mini-1)) = cati(1);
            allPlates(j:(j+mini-1)) = drugs_in_common{i,4};
            allTargets(j:(j+mini-1)) = drugs_in_common{i,5};
            allPathways(j:(j+mini-1)) = drugs_in_common{i,6};
            j = j+mini;
        end
        
    case 3
        line1_drugInfo = table(line1.drug_names,...
            line1.dose, line1.drug_categories,...
            line1.plate_types, line1.targets,...
            line1.pathways);
        line1_un = unique(line1_drugInfo, 'rows');

        line2_drugInfo = table(line2.drug_names,...
            line2.dose, line2.drug_categories,...
            line2.plate_types, line2.targets,...
            line2.pathways);
        line2_un = unique(line2_drugInfo, 'rows');
        
        line3_drugInfo = table(line3.drug_names,...
            line3.dose, line3.drug_categories,...
            line3.plate_types, line3.targets,...
            line3.pathways);
        line3_un = unique(line3_drugInfo, 'rows');
        
        a = getRows(line1_un, line2_un, line3_un);
        
        drugs_in_common = line1_un(a,:);
        allProfiles = zeros(size(line1.profiles{1},1), size(line1.profiles{1},2)*3);
        allDrugs = repmat({''}, size(line1.profiles,1),1); 
        allDose = repmat({''}, size(line1.profiles,1),1);
        allCats = repmat({''}, size(line1.profiles,1),1);
        allPlates = repmat({''}, size(line1.profiles,1),1);
        allTargets = repmat({''}, size(line1.profiles,1),1);
        allPathways = repmat({''}, size(line1.profiles,1),1);
        
        j=1;
        for i = 1:size(drugs_in_common)
            drugi = drugs_in_common{i,1};
            dosei = drugs_in_common{i,2};
            cati = drugs_in_common{i,3};
            line1_rows = (strcmp(line1.drug_names, drugi) & strcmp(line1.dose, dosei));
            line2_rows = (strcmp(line2.drug_names, drugi) & strcmp(line2.dose, dosei));
            line3_rows = (strcmp(line3.drug_names, drugi) & strcmp(line3.dose, dosei));
            mini = min([sum(line1_rows), sum(line2_rows), sum(line3_rows)]);

            line1_profs = cell2mat(line1.profiles(line1_rows,:));
            line2_profs = cell2mat(line2.profiles(line2_rows,:));
            line3_profs = cell2mat(line3.profiles(line3_rows,:));

            allProfiles(j:(j+mini-1),:) = [line1_profs(1:mini,:) line2_profs(1:mini,:) line3_profs(1:mini,:)];
            allDrugs(j:(j+mini-1)) = drugi(1);
            allDose(j:(j+mini-1)) = dosei(1);
            allCats(j:(j+mini-1)) = cati(1);
            allPlates(j:(j+mini-1)) = drugs_in_common{i,4};
            allTargets(j:(j+mini-1)) = drugs_in_common{i,5};
            allPathways(j:(j+mini-1)) = drugs_in_common{i,6};
            j = j+mini;
        end
        
    case 4
        line1_drugInfo = table(line1.drug_names,...
            line1.dose, line1.drug_categories,...
            line1.plate_types, line1.targets,...
            line1.pathways);
        line1_un = unique(line1_drugInfo, 'rows');

        line2_drugInfo = table(line2.drug_names,...
            line2.dose, line2.drug_categories,...
            line2.plate_types, line2.targets,...
            line2.pathways);
        line2_un = unique(line2_drugInfo, 'rows');
        
        line3_drugInfo = table(line3.drug_names,...
            line3.dose, line3.drug_categories,...
            line3.plate_types, line3.targets,...
            line3.pathways);
        line3_un = unique(line3_drugInfo, 'rows');
        
        line4_drugInfo = table(line4.drug_names,...
            line4.dose, line4.drug_categories,...
            line4.plate_types, line4.targets,...
            line4.pathways);
        line4_un = unique(line4_drugInfo, 'rows');
        a = getRows(line1_un, line2_un, line3_un, line4_un);
        
        drugs_in_common = line1_un(a,:);
        allProfiles = zeros(size(line1.profiles{1},1), size(line1.profiles{1},2)*4);
        allDrugs = repmat({''}, size(line1.profiles,1),1); 
        allDose = repmat({''}, size(line1.profiles,1),1);
        allCats = repmat({''}, size(line1.profiles,1),1);
        allPlates = repmat({''}, size(line1.profiles,1),1);
        allTargets = repmat({''}, size(line1.profiles,1),1);
        allPathways = repmat({''}, size(line1.profiles,1),1);
        
        j=1;
        for i = 1:size(drugs_in_common)
            drugi = drugs_in_common{i,1};
            dosei = drugs_in_common{i,2};
            cati = drugs_in_common{i,3};
            line1_rows = (strcmp(line1.drug_names, drugi) & strcmp(line1.dose, dosei));
            line2_rows = (strcmp(line2.drug_names, drugi) & strcmp(line2.dose, dosei));
            line3_rows = (strcmp(line3.drug_names, drugi) & strcmp(line3.dose, dosei));
            line4_rows = (strcmp(line4.drug_names, drugi) & strcmp(line4.dose, dosei));
            mini = min([sum(line1_rows), sum(line2_rows), sum(line3_rows), sum(line4_rows)]);

            line1_profs = cell2mat(line1.profiles(line1_rows,:));
            line2_profs = cell2mat(line2.profiles(line2_rows,:));
            line3_profs = cell2mat(line3.profiles(line3_rows,:));
            line4_profs = cell2mat(line4.profiles(line4_rows,:));
            
            allProfiles(j:(j+mini-1),:) = [line1_profs(1:mini,:) line2_profs(1:mini,:) line3_profs(1:mini,:) line4_profs(1:mini,:)];
            allDrugs(j:(j+mini-1)) = drugi(1);
            allDose(j:(j+mini-1)) = dosei(1);
            allCats(j:(j+mini-1)) = cati(1);
            allPlates(j:(j+mini-1)) = drugs_in_common{i,4};
            allTargets(j:(j+mini-1)) = drugs_in_common{i,5};
            allPathways(j:(j+mini-1)) = drugs_in_common{i,6};
            j = j+mini;
        end
        
end

allProfiles = allProfiles(1:j-1, :);
allDrugs = allDrugs(1:j-1, :);
allDose = allDose(1:j-1, :);
allCats = allCats(1:j-1, :);
allPlates = allPlates(1:j-1, :);
allTargets = allTargets(1:j-1, :);
allPathways = allPathways(1:j-1, :);
catStruct = struct('profiles', {num2cell(allProfiles,2)},...
    'drug_names', {allDrugs},...
    'dose', {allDose},...
    'drug_categories', {allCats},...
    'plate_types', {allPlates},...
    'targets', {allTargets},...
    'pathways', {allPathways});

end