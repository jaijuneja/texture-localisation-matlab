function [index, images, path] = build_index(imgFolder, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
% 
% BUILD INDEX
% [index, images, path] = build_index(imgFolder, 'imgFormats',
% valImgFormats, 'numWords', valNumWords)
%
% Generates index of images using all the images in a given folder. Can 
% specify the number of words using the 'numWords' property.
%
% Inputs:
%   - imgFolder:    String containing path to image folder
%   
%   Optional Properties:
%       - 'imgFormats': Cell array containing strings of valid image
%                       formats - e.g {'.jpg', '.png', '.tiff'}. By default
%                       accepts JPEG and PNG images
%       - 'numWords':   Number of visual words to use in inverted index
%
% Outputs:
%   - index:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%   - images:   For n images, a 1xn cell array where each element contains
%               the path to each file

if ~ischar(imgFolder) || ~isdir(imgFolder)
    error('Input "imgFolder" must be a string containing a path')
elseif ~isequal(imgFolder(end), '/')
    imgFolder = strcat(imgFolder, '/'); % Path needs to end with /
end

opts.numWords = 1e5;
opts.imgFormats = {'.jpg', '.png'};
opts = vl_argparse(opts, varargin);

path.images = imgFolder;
path.etc = strcat(path.images, 'etc/');
path.index = strcat(path.etc, 'index.mat');
path.cor = strcat(path.etc, 'cor.mat');
path.world = strcat(path.etc, 'world.mat');

% Obtain listing of all files in image path
image_files = dir(path.images);

% Initiate blank cell array of image paths
images = cell(1, length(image_files));
not_img = [];

% Populate cell array with paths to images
for i = 1:length(image_files)
    if ~image_files(i).isdir && ...
            sum(strcmp(image_files(i).name(end-3:end), opts.imgFormats))
        % If file is not a directory and ends with one of the extensions
        % listed in opts.imgFormats then add its path to images(i)
        images(i) = { strcat(path.images, image_files(i).name) };    
    else
        % If the listed file is not an image, record its index
        not_img = [not_img i];
    end
end
images(not_img) = []; % Remove all entries that are directories (not images)

% Iinitialize the index structure and estimate visual word vocabulary
% Save the generated index in the etc path
if ~exist(path.etc, 'dir')
    mkdir(path.etc);
end
% If index doesn't already exist, generate it
if ~exist(path.index, 'file')
    index = visualindex_build(images, 'numWords', opts.numWords);
    % Add images to index
    index = visualindex_add(index, images, 1:length(images));
    % Save index
    save(path.index, 'index');
else
    index = load(path.index);
    index = index.index;
end

end