for i=1:6
a = vertcat(classification_accuracy.XRCC5.cv_report{1,i}{1,1},classification_accuracy.XRCC5.cv_report{1,i}{2,1},classification_accuracy.XRCC5.cv_report{1,i}{3,1},classification_accuracy.XRCC5.cv_report{1,i}{4,1},classification_accuracy.XRCC5.cv_report{1,i}{4,1},classification_accuracy.XRCC5.cv_report{1,i}{5,1},classification_accuracy.XRCC5.cv_report{1,i}{6,1},classification_accuracy.XRCC5.cv_report{1,i}{7,1},classification_accuracy.XRCC5.cv_report{1,i}{8,1},classification_accuracy.XRCC5.cv_report{1,i}{9,1},classification_accuracy.XRCC5.cv_report{1,i}{10,1},classification_accuracy.XRCC5.cv_report{1,i}{11,1});
xrcc5_aw = cell2table(a);
xrcc5_aw.Properties.VariableNames{1} = 'Truth';
xrcc5_aw.Properties.VariableNames{1} = 'Label';
xrcc5_aw.Properties.VariableNames{2} = 'Truth';
xrcc5_aw.Properties.VariableNames{3} = 'Prediction';
xrcc5_aw.Properties.VariableNames{4} = 'Iteration';
xrcc5_aw.log = strcmpi(xrcc5_aw.Prediction,xrcc5_aw.Truth);
xrcc5_aw = xrcc5_aw(xrcc5_aw.log == 0,:);
accuracy.XRCC5(i).offenders = cell2table(tabulate(xrcc5_aw.Label));
accuracy.XRCC5(i).offenders = accuracy.XRCC5(i).offenders(accuracy.XRCC5(i).offenders.Var2 > 1,:);
accuracy.XRCC5(i).offenders = sortrows(accuracy.XRCC5(i).offenders,'Var1','ascend');
accuracy.XRCC5(i).offenders(:,'Var3') = [];
end

for i=1:6
a = vertcat(classification_accuracy.SET.cv_report{1,i}{1,1},classification_accuracy.SET.cv_report{1,i}{2,1},classification_accuracy.SET.cv_report{1,i}{3,1},classification_accuracy.SET.cv_report{1,i}{4,1},classification_accuracy.SET.cv_report{1,i}{4,1},classification_accuracy.SET.cv_report{1,i}{5,1},classification_accuracy.SET.cv_report{1,i}{6,1},classification_accuracy.SET.cv_report{1,i}{7,1},classification_accuracy.SET.cv_report{1,i}{8,1},classification_accuracy.SET.cv_report{1,i}{9,1},classification_accuracy.SET.cv_report{1,i}{10,1},classification_accuracy.SET.cv_report{1,i}{11,1});
SET_aw = cell2table(a);
SET_aw.Properties.VariableNames{1} = 'Truth';
SET_aw.Properties.VariableNames{1} = 'Label';
SET_aw.Properties.VariableNames{2} = 'Truth';
SET_aw.Properties.VariableNames{3} = 'Prediction';
S100A11_aw.Properties.VariableNames{4} = 'Iteration';
SET_aw.log = strcmpi(SET_aw.Prediction,SET_aw.Truth);
SET_aw = SET_aw(SET_aw.log == 0,:);
accuracy.SET(i).offenders = cell2table(tabulate(SET_aw.Label));
accuracy.SET(i).offenders = accuracy.SET(i).offenders(accuracy.SET(i).offenders.Var2 > 1,:);
accuracy.SET(i).offenders = sortrows(accuracy.SET(i).offenders,'Var1','ascend');
accuracy.SET(i).offenders(:,'Var3') = [];
end

for i=1:6
a = vertcat(classification_accuracy.S100A11.cv_report{1,i}{1,1},classification_accuracy.S100A11.cv_report{1,i}{2,1},classification_accuracy.S100A11.cv_report{1,i}{3,1},classification_accuracy.S100A11.cv_report{1,i}{4,1},classification_accuracy.S100A11.cv_report{1,i}{4,1},classification_accuracy.S100A11.cv_report{1,i}{5,1},classification_accuracy.S100A11.cv_report{1,i}{6,1},classification_accuracy.S100A11.cv_report{1,i}{7,1},classification_accuracy.S100A11.cv_report{1,i}{8,1},classification_accuracy.S100A11.cv_report{1,i}{9,1},classification_accuracy.S100A11.cv_report{1,i}{10,1},classification_accuracy.S100A11.cv_report{1,i}{11,1});
S100A11_aw = cell2table(a);
S100A11_aw.Properties.VariableNames{1} = 'Truth';
S100A11_aw.Properties.VariableNames{1} = 'Label';
S100A11_aw.Properties.VariableNames{2} = 'Truth';
S100A11_aw.Properties.VariableNames{3} = 'Prediction';
S100A11_aw.Properties.VariableNames{4} = 'Iteration';
S100A11_aw.log = strcmpi(S100A11_aw.Prediction,S100A11_aw.Truth);
S100A11_aw = S100A11_aw(S100A11_aw.log == 0,:);
accuracy.S100A11(i).offenders = cell2table(tabulate(S100A11_aw.Label));
accuracy.S100A11(i).offenders = accuracy.S100A11(i).offenders(accuracy.S100A11(i).offenders.Var2 > 1,:);
accuracy.S100A11(i).offenders = sortrows(accuracy.S100A11(i).offenders,'Var1','ascend');
accuracy.S100A11(i).offenders(:,'Var3') = [];
end

for i=1:6
a = vertcat(classification_accuracy.NQO1.cv_report{1,i}{1,1},classification_accuracy.NQO1.cv_report{1,i}{2,1},classification_accuracy.NQO1.cv_report{1,i}{3,1},classification_accuracy.NQO1.cv_report{1,i}{4,1},classification_accuracy.NQO1.cv_report{1,i}{4,1},classification_accuracy.NQO1.cv_report{1,i}{5,1},classification_accuracy.NQO1.cv_report{1,i}{6,1},classification_accuracy.NQO1.cv_report{1,i}{7,1},classification_accuracy.NQO1.cv_report{1,i}{8,1},classification_accuracy.NQO1.cv_report{1,i}{9,1},classification_accuracy.NQO1.cv_report{1,i}{10,1},classification_accuracy.NQO1.cv_report{1,i}{11,1});
NQO1_aw = cell2table(a);
NQO1_aw.Properties.VariableNames{1} = 'Truth';
NQO1_aw.Properties.VariableNames{1} = 'Label';
NQO1_aw.Properties.VariableNames{2} = 'Truth';
NQO1_aw.Properties.VariableNames{3} = 'Prediction';
NQO1_aw.Properties.VariableNames{4} = 'Iteration';
NQO1_aw.log = strcmpi(NQO1_aw.Prediction,NQO1_aw.Truth);
NQO1_aw = NQO1_aw(NQO1_aw.log == 0,:);
accuracy.NQO1(i).offenders = cell2table(tabulate(NQO1_aw.Label));
accuracy.NQO1(i).offenders = accuracy.NQO1(i).offenders(accuracy.NQO1(i).offenders.Var2 > 1,:);
accuracy.NQO1(i).offenders = sortrows(accuracy.NQO1(i).offenders,'Var1','ascend');
accuracy.NQO1(i).offenders(:,'Var3') = [];
end

load cpd_list.mat
xrcc5 = cell(6,166);
nqo1 = cell(6,166)
set = cell(6,166)
s100a11 = cell(6,166)
for i=1:6
    for j = 1:size(cpd_list,1)
        xrcc5{i,j} = (ismember(accuracy.XRCC5(i).offenders.Var1,cpd_list{j}))';
        xrcc5{i,j} = xrcc5{i,j}((ismember(accuracy.XRCC5(i).offenders.Var1,cpd_list{j}))>0)';
        set{i,j} = (ismember(accuracy.SET(i).offenders.Var1,cpd_list{j}))';
        set{i,j} = set{i,j}((ismember(accuracy.SET(i).offenders.Var1,cpd_list{j}))>0)';
        s100a11{i,j} = (ismember(accuracy.S100A11(i).offenders.Var1,cpd_list{j}))';
        s100a11{i,j} = s100a11{i,j}((ismember(accuracy.S100A11(i).offenders.Var1,cpd_list{j}))>0)';
        nqo1{i,j} = (ismember(accuracy.NQO1(i).offenders.Var1,cpd_list{j}))';
        nqo1{i,j} = nqo1{i,j}((ismember(accuracy.NQO1(i).offenders.Var1,cpd_list{j}))>0)';
        
        if isempty(xrcc5{i,j})== 1
            xrcc5_m(i,j)=0
        else
            xrcc5_m(i,j)=1
        end
        
        if isempty(set{i,j})== 1
            set_m(i,j)=0
        else
            set_m(i,j)=1
        end
        
        if isempty(s100a11{i,j})== 1
            s100a11_m(i,j)=0
        else
            s100a11_m(i,j)=1
        end
        
        if isempty(nqo1{i,j})== 1
            nqo1_m(i,j)=0
        else
            nqo1_m(i,j)=1
        end
    end
end
% figure
% imagesc(xrcc5_m)
% title('XRCC5')
% figure
% imagesc(set_m)
% title('SET')
% figure
% imagesc(s100a11_m)
% title('S100A11')
% figure
% imagesc(nqo1_m)
% title('NQO1')
allplates = vertcat(xrcc5_m,set_m,nqo1_m,s100a11_m)
imagesc(allplates)