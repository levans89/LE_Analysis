function plate = Select_Rows( plate, chosen_flag )
% Select rows of data in a given plate struct
% 
% plate = Select_Rows( plate, chosen_flag )
% 
%  - chosen_row is a logical vector with the length = size(plate.profiles,1)
% 
% Copyright: Altschuler and Wu laboratories
% Author: Chien-Hsiang Hsu at the Altschuler and Wu Laboratories
% For latest updates, check: <http://www.altschulerwulab.org/>.
%
% All rights reserved.


if ~islogical(chosen_flag)
    error('Chosen_row should be logical.')
end

if length(chosen_flag)~=size(plate.profiles,1)
    error('Chosen_row should have the same length as number of profiles.')
end

f_names = fieldnames(plate);

for f = 1:length(f_names)
    if size(plate.(f_names{f}),1)==length(chosen_flag) && ~strcmp(f_names{f},'feaInfo')
        plate.(f_names{f}) = plate.(f_names{f})(chosen_flag,:,:);
    end
end

end

