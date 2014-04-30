function resize_images(imgFolder, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
%
% RESIZE_IMAGES
% resize_images(imgFolder, 'maxDimension', valMaxDimension, 'fileSuffix',
% valFileSuffix, 'saveFolder', valSaveFolder)
%
% Resizes all images in a folder to a maximum specified dimension. Images
% are saved to a new folder by default. To overwrite the old images set the
% property 'saveFolder' to the same value as imgFolder.
%
% Inputs:
%   - imgFolder:    String containing path to a folder with images to be
%                   resized
%
%   Optional Properties:
%       - 'maxDimension':	Sets maximum number of pixels along major axis
%                           of image. Set to 1024 by default
%       - 'fileSuffix':     String to be added to the end of the original
%                           file names
%       - 'saveFolder':     String containing the path to the folder where
%                           images are to be saved. Set equal to imgFolder
%                           to save in original folder or to overwrite
%                           files. By default valSavePath = imgFolder +
%                           '_resized'

% Some initial error checking of input arguments
if ~ischar(imgFolder) || ~isdir(imgFolder)
    error('Input "imgFolder" must be a string containing a path')
elseif ~isequal(imgFolder(end), filesep)
    imgFolder(end+1) = filesep; % Path needs to end with /
end

% Default property values
opts.saveFolder = strcat(imgFolder(1:end-1), '_resized/');
opts.maxDimension = 1024;
opts.fileSuffix = '';
% Parse optional input arguments
opts = vl_argparse(opts, varargin);

% Error checking of optional properties
if ~ischar(opts.fileSuffix)
    error('Property "fileSuffix" must be a string')
end
if ~ischar(opts.saveFolder)
    error('Property "saveFolder" must be a string')
elseif ~isequal(opts.saveFolder(end), filesep)
    opts.saveFolder(end+1) = filesep; % Path needs to end with /
end

% Obtain listing of all files in image folder
image_files = dir(imgFolder);

% Populate cell array with paths to images
for i = 1:length(image_files)
    if ~image_files(i).isdir
        impath = strcat(imgFolder, image_files(i).name);
        try
            im = imread(impath);
        catch
            continue;
        end
        [height, width, ~] = size(im);
        isLandscape = width > height;
        if isLandscape
            im = imresize(im, [NaN, opts.maxDimension]);
        else
            im = imresize(im, [opts.maxDimension, NaN]);
        end
        % Save resized image
        if ~exist(opts.saveFolder, 'dir')
            mkdir(opts.saveFolder);
        end
        imwrite(im, strcat(opts.saveFolder, image_files(i).name, ...
            opts.fileSuffix));
    end
end

end