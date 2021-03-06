function colors = Drug_Category_ColorLE(GN)
% Consistent colors for labeling drugs
% y = my_color(GN)
%
% Copyright: Altschuler and Wu laboratories
% Author: Chien-Hsiang Hsu at the Altschuler and Wu Laboratories
% For latest updates, check: <http://www.altschulerwulab.org/>.
%
% All rights reserved.

GN(strcmpi(GN,'Proteosome')) = {'Proteasome'};
GN(strcmpi(GN,'Hsp90')) = {'HSP90'};
GN(strcmpi(GN,'unclassified')) = {'Unclassified'};

colors = zeros(length(GN),3);

% Specified color
    color_spec = [0.8       0.8       0.8       % DMSO
                  0         0         1         % DNA              
                  0         0         1         % DNA Damage
                  1.0       0         0         % Proteasome
                  1.0       0         0         % Proteases
                  0         1.0       0         % HDAC
                  0         1.0       0         % Epigenetics
                  0.5961    0.3059    0.6392    % Hsp90 (purple)
                  1.0000    0.8276    0         % MT
                  0    0.3448         0         % mTOR
                  0    0.3448         0         % PI3K/Akt/mTOR
                  0    0.3448         0         % AKT
                  0.7       0.7       0.7       % Unclassified
                  0.7       0.7       0.7       % Cytoskeletal Signaling
                  0.7       0.7       0.7       % Others       
                  1.0000    0.8276    0         % Microbiology
                  0.6207    0.3103    0.2759    % Actin
                  0    1.0000    0.7586         % AuroraB
                  0    1.0000    0.7586         % Cell Cycle
                  0    1.0000    0.7586         % PLK
                  0    0.5172    0.5862         % ER
                  0.5862    0.8276    0.3103    % Endocrinology & Hormones
                  0.3       0.3       0.3       % Positive Control
                  0.9922    0.7059    0.3843    % Apoptosis
                  0.8471    0.7020    0.3961    % Neuronal Signaling
                  1.0       0         0         % Ubiquitin
                  0.5529    0.8275    0.7804    % Angiogenesis
                  0.3294	0.1882	0.0196      % Stem Cells &  Wnt',...
                  0.5490	0.3176	0.0392      % Protein Tyrosine Kinase',...
                  0.7490	0.5059	0.1765      % JAK/STAT',...
                  0.8745	0.7608	0.4902      % MAPK',...
                  0.9647	0.9098	0.7647      % Metabolism',...
                  0.7804	0.9176	0.8980      % NF-?B',...
                  0.5020	0.8039	0.7569      % TGF-beta/Smad',...
                  0.2078	0.5922	0.5608      % GPCR & G Protein',...
                  0.0039	0.4000	0.3686      % Transmembrane Transporters'
                  ];

color_order =  {'DMSO',...
                'DNA',...
                'DNA Damage',...
                'Proteasome',...
                'Proteases',...                
                'HDAC',...
                'Epigenetics',...
                'HSP90',...
                'MT',...
                'mTOR',...
                'PI3K/Akt/mTOR',...
                'AKT',...
                'Unclassified',...
                'Cytoskeletal Signaling',...
                'Others',...
                'Microbiology',...
                'Actin',...
                'AuroraB',...
                'Cell Cycle',...
                'PLK',...
                'ER',...
                'Endocrinology & Hormones',...
                'Positive Control',...
                'Apoptosis',...
                'Neuronal Signaling',...
                'Ubiquitin',...
                'Angiogenesis',...
                'Stem Cells &  Wnt',...
                'Protein Tyrosine Kinase',...
                'JAK/STAT',...
                'MAPK',...
                'Metabolism',...
                'NF-?B',...
                'TGF-beta/Smad',...
                'GPCR & G Protein',...
                'Transmembrane Transporters'};            
 
            
[is_spec,loc] = ismember(GN,color_order);
colors(is_spec,:) = color_spec(loc(is_spec),:);


% Assign random color to the categories without specified color
Nmiss = sum(~is_spec);
if Nmiss < 10
    colors(~is_spec,:) = Distinct_Colors(Nmiss,colors(is_spec,:));
else
    colors(~is_spec,:) = hsv(Nmiss);
end
end


%== Helper functions ==%
function new_colors = Distinct_Colors(n,used_colors)

new_colors = zeros(n,3);

% Generate candidate colors
n_grid = 100;  % number of grid divisions along each axis in RGB space
x = linspace(0,1,n_grid);
[R,G,B] = ndgrid(x,x,x);
candidates = [R(:) G(:) B(:)];

for i = 1:n
    min_D = min(pdist2(candidates,used_colors),[],2);
    [~,idx] = max(min_D);
    if isempty(idx)
        idx = 1;
    end
    new_colors(i,:) = candidates(idx,:);
    used_colors = [used_colors;new_colors(i,:)];
end

end

