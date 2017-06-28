load('V:\projects\2010\2010_06_CD_Tag\code\NBT_ORACL_code_and_data_All\data\plate_2014002036.mat')

qt6 = plate.data(:,:,6);
qt5 = plate.data(:,:,5);
qt4 = plate.data(:,:,4);
qt3 = plate.data(:,:,3);
qt2 = plate.data(:,:,2);
%qt1 = plate.data(:,:,1);

%NSC26113_1(1,:) = qt1{117}(1,:);
%NSC33570_1(1,:) = qt1{167}(1,:);
NSC26113_2(1,:) = qt2{117}(1,:);
NSC33570_2(1,:) = qt2{167}(1,:);
NSC26113_3(1,:) = qt3{117}(1,:);
NSC33570_3(1,:) = qt3{167}(1,:);
NSC26113_4(1,:) = qt4{117}(1,:);
NSC33570_4(1,:) = qt4{167}(1,:);
NSC26113_5(1,:) = qt5{117}(1,:);
NSC33570_5(1,:) = qt5{167}(1,:);
NSC26113_6(1,:) = qt6{117}(1,:);
NSC33570_6(1,:) = qt6{167}(1,:);

%NSC26113 = horzcat(NSC26113_1,NSC26113_2,NSC26113_3,NSC26113_4,NSC26113_5,NSC26113_6);
%NSC33570 = horzcat(NSC33570_1,NSC33570_2,NSC33570_3,NSC33570_4,NSC33570_5,NSC33570_6);
NSC26113 = horzcat(NSC26113_2,NSC26113_3,NSC26113_4,NSC26113_5,NSC26113_6);
NSC33570 = horzcat(NSC33570_2,NSC33570_3,NSC33570_4,NSC33570_5,NSC33570_6);

for a = 1:22
%DMSO_1(a,:) = qt1{a};
DMSO_2(a,:) = qt2{a};
DMSO_3(a,:) = qt3{a};
DMSO_4(a,:) = qt4{a};
DMSO_5(a,:) = qt5{a};
DMSO_6(a,:) = qt6{a};
end
%DMSO = horzcat(DMSO_1,DMSO_2,DMSO_3,DMSO_4,DMSO_5,DMSO_6);
DMSO = horzcat(DMSO_2,DMSO_3,DMSO_4,DMSO_5,DMSO_6);

load('V:\projects\2010\2010_06_CD_Tag\code\NBT_ORACL_code_and_data_All\data\plate_2014002034.mat')
rt6 = plate.data(:,:,6);
rt5 = plate.data(:,:,5);
rt4 = plate.data(:,:,4);
rt3 = plate.data(:,:,3);
rt2 = plate.data(:,:,2);
%rt1 = plate.data(:,:,1);

for a = 16:62
%DMSO_r1(a-15,:) = rt1{a};
DMSO_r2(a-15,:) = rt2{a};
DMSO_r3(a-15,:) = rt3{a};
DMSO_r4(a-15,:) = rt4{a};
DMSO_r5(a-15,:) = rt5{a};
DMSO_r6(a-15,:) = rt6{a};
end
%DMSO = horzcat(DMSO_r1,DMSO_r2,DMSO_r3,DMSO_r4,DMSO_r5,DMSO_r6);
DMSO = horzcat(DMSO_r2,DMSO_r3,DMSO_r4,DMSO_r5,DMSO_r6);


for a = 194:213
%proteasome_r1(a,:) = rt1{a};
proteasome_r2(a,:) = rt2{a};
proteasome_r3(a,:) = rt3{a};
proteasome_r4(a,:) = rt4{a};
proteasome_r5(a,:) = rt5{a};
proteasome_r6(a,:) = rt6{a};
end

for a = 230:234
%proteasome_r1(a,:) = rt1{a};
proteasome_r2(a,:) = rt2{a};
proteasome_r3(a,:) = rt3{a};
proteasome_r4(a,:) = rt4{a};
proteasome_r5(a,:) = rt5{a};
proteasome_r6(a,:) = rt6{a};
end

%proteasome_r = horzcat(proteasome_r1,proteasome_r2,proteasome_r3,proteasome_r4,proteasome_r5,proteasome_r6);
proteasome_r = horzcat(proteasome_r2,proteasome_r3,proteasome_r4,proteasome_r5,proteasome_r6);

for a = 1:15
%ref_cpds1(a,:) = rt1{a};
ref_cpds2(a,:) = rt2{a};
ref_cpds3(a,:) = rt3{a};
ref_cpds4(a,:) = rt4{a};
ref_cpds5(a,:) = rt5{a};
ref_cpds6(a,:) = rt6{a};
end

for a = 63:77
%ref_cpds1(a,:) = rt1{a};
ref_cpds2(a,:) = rt2{a};
ref_cpds3(a,:) = rt3{a};
ref_cpds4(a,:) = rt4{a};
ref_cpds5(a,:) = rt5{a};
ref_cpds6(a,:) = rt6{a};
end

for a = 94:193
%ref_cpds1(a,:) = rt1{a};
ref_cpds2(a,:) = rt2{a};
ref_cpds3(a,:) = rt3{a};
ref_cpds4(a,:) = rt4{a};
ref_cpds5(a,:) = rt5{a};
ref_cpds6(a,:) = rt6{a};
end

for a = 240:264
%ref_cpds1(a,:) = rt1{a};
ref_cpds2(a,:) = rt2{a};
ref_cpds3(a,:) = rt3{a};
ref_cpds4(a,:) = rt4{a};
ref_cpds5(a,:) = rt5{a};
ref_cpds6(a,:) = rt6{a};
end

%ref_cpds = horzcat(ref_cpds1,ref_cpds2,ref_cpds3,ref_cpds4,ref_cpds5,ref_cpds6);
ref_cpds = horzcat(ref_cpds2,ref_cpds3,ref_cpds4,ref_cpds5,ref_cpds6);

all = vertcat(ref_cpds,proteasome_r,NSC26113,NSC33570);
imagesc(all)
colormap(jet)
