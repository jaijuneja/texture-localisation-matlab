% Resizes all images in a folder to the maximum dimension specified below
% Images can either be overwritten or saved to a different path

% First get paths to all images in folder
opts.image_path = 'test_images/outside_tiling/';
opts.save_path = 'test_images/outside_tiling_lowres/';
opts.max_dimension = 1024; % Maximum size is 1024 pixels height or width
opts.file_suffix = ''; % If left as empty string, new images will have same names

% Obtain listing of all files in image path
image_files = dir(opts.image_path);

% Initiate blank cell array of image paths
images = cell(1, length(image_files));

% Populate cell array with paths to images
for i = 1:length(image_files)
    if ~image_files(i).isdir
        impath = strcat(opts.image_path, image_files(i).name);
        im = imread(impath);
        [height, width, ~] = size(im);
        isLandscape = width > height;
        if isLandscape
            im = imresize(im, [NaN, opts.max_dimension]);
        else
            im = imresize(im, [opts.max_dimension, NaN]);
        end
        % Save resized image
        if ~exist(opts.save_path, 'dir')
            mkdir(opts.save_path);
        end
        imwrite(im, strcat(opts.save_path, image_files(i).name, opts.file_suffix));
    end
end