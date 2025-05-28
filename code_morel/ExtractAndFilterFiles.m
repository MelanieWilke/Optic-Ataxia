function [fullfilenames] = ExtractAndFilterFiles(directory, keyword)
% Extracts all files from a given directory
% and returns a sorted list of the full pathname of all files
% filtered for keyword
if nargin < 1
    directory = pwd;
end
if nargin < 2 
    keyword = '';
end
% extract filenames and directorynames
files = dir(directory);
files = files(~[files.isdir] & contains({files.name},keyword));
fullfilenames = arrayfun(@(x) [x.folder filesep x.name], files, 'UniformOutput', false);
end
