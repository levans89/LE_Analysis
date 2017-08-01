for s = 1:length(st)
st{s,1}=num2str(st{s,1})
end
extract = Exp_DB(ismember(Exp_DB.expt_plate,st),:)
masterplates = plates;
plates = extract;