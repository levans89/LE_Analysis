function data = importfile(workbookFile, sheetName, range)
%IMPORTFILE Import data from a spreadsheet
%   DATA = IMPORTFILE(FILE) reads all numeric data from the first worksheet
%   in the Microsoft Excel spreadsheet file named FILE and returns the
%   numeric data.
%
%   DATA = IMPORTFILE(FILE,SHEET) reads from the specified worksheet.
%
%   DATA = IMPORTFILE(FILE,SHEET,RANGE) reads from the specified worksheet
%   and from the specified RANGE. Specify RANGE using the syntax
%   'C1:C2',where C1 and C2 are opposing corners of the region.
%
%	Non-numeric cells are replaced with: NaN
%
% Example:
%   Untitled = importfile('2016018104.xlsx','BADWELLS','A1:Y17');
%
%   See also XLSREAD.

% Auto-generated by MATLAB on 2017/03/13 11:30:58

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If no range is specified, read all data
if nargin <= 2 || isempty(range)
    range = '';
end

%% Import the data
[~, ~, raw] = xlsread(workbookFile, sheetName, range);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {0}; % Replace non-numeric cells
data = reshape(raw,[384,1]);
clear raw R

