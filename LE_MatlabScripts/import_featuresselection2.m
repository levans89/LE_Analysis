function [name,category_id,channel_id,use,use_2ch,use_all] = import_featuresselection2(filename, startRow, endRow)
%IMPORTFILE1 Import numeric data from a text file as column vectors.
%   [NAME1,CATEGORY_ID1,CHANNEL_ID1,USE1,USE_2CH,USE_ALL] =
%   IMPORTFILE1(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [NAME1,CATEGORY_ID1,CHANNEL_ID1,USE1,USE_2CH,USE_ALL] =
%   IMPORTFILE1(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [name1,category_id1,channel_id1,use1,use_2ch,use_all] = importfile1('feature_select_NBT.txt',2, 660);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/07/11 13:43:32

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: text (%q)
%	column2: double (%f)
%   column3: text (%q)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%f%q%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
name = dataArray{:, 1};
category_id = dataArray{:, 2};
channel_id = dataArray{:, 3};
use = dataArray{:, 4};
use_2ch = dataArray{:, 5};
use_all = dataArray{:, 6};


