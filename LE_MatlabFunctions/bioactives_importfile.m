function tableout = bioactives_importfile(workbookFile,sheetName,startRow,endRow)
%IMPORTFILE Import data from a spreadsheet
%   DATA = IMPORTFILE(FILE) reads data from the first worksheet in the
%   Microsoft Excel spreadsheet file named FILE and returns the data as a
%   table.
%
%   DATA = IMPORTFILE(FILE,SHEET) reads from the specified worksheet.
%
%   DATA = IMPORTFILE(FILE,SHEET,STARTROW,ENDROW) reads from the specified
%   worksheet for the specified row interval(s). Specify STARTROW and
%   ENDROW as a pair of scalars or vectors of matching size for
%   dis-contiguous row intervals. To read to the end of the file specify an
%   ENDROW of inf.
%
%	Non-numeric cells are replaced with: NaN
%
% Example:
%   XRCC5bioactive = importfile('XRCC5_bioactive.xlsx','Sheet1',1,3528);


%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 3
    startRow = 2;
    endRow = 4000;
end

%% Import the data
[~, ~, raw] = xlsread(workbookFile, sheetName, sprintf('A%d:M%d',startRow(1),endRow(1)));
for block=2:length(startRow)
    [~, ~, tmpRawBlock] = xlsread(workbookFile, sheetName, sprintf('A%d:M%d',startRow(block),endRow(block)));
    raw = [raw;tmpRawBlock]; %#ok<AGROW>
end
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[2,3,4,8,9,10,11,12]);
raw = raw(:,[1,5,6,7,13]);
profiles = raw(:,[14:end]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
I = cellfun(@(x) ischar(x), raw);
raw(I) = {NaN};
data = reshape([raw{:}],size(raw));

%% Create table
tableout = table;

%% Allocate imported array to column variable names
tableout.plateIDs = data(:,1);
tableout.well_names = cellVectors(:,1);
tableout.drug_categories = cellVectors(:,2);
tableout.drug_names = cellVectors(:,3);
tableout.dose = data(:,2);
tableout.concentrations = data(:,3);
tableout.bioactive_pVals = data(:,4);
tableout.targets = cellVectors(:,4);
tableout.pathways = cellVectors(:,5);
tableout.descriptions = cellVectors(:,6);
tableout.cas = cellVectors(:,7);
tableout.cpd_usage = cellVectors(:,8);
tableout.nCells = data(:,5);


