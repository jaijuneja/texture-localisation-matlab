%% Initialise: get paths to all images in database
% Look at demo_matchFeatures - basically implement that code in a loop
opts.image_path = 'test_images/set1_lowres/';
opts.peak_thresh = 0;
opts.edge_thresh = 0;

% Obtain listing of all files in image path
image_files = dir(opts.image_path);

% Initiate blank cell array of image paths
images = cell(1, length(image_files));
features = cell(3, length(image_files));
isdir = [];
IDcounter = 0;
% Obtain features for all images
for i = 1:length(image_files)
    if ~image_files(i).isdir
        IDcounter = IDcounter + 1;
        impath = strcat(opts.image_path, image_files(i).name);
        try
            im = imread(impath);
        catch
            continue; % If not an image, jump to next file
        end
        images{i} = im;
        % make single
        im = im2single(im) ;
        if size(im,3) > 1, im = rgb2gray(im); end
        [features{2, i}, features{3, i}] = vl_sift(im);
        features{1, i} = IDcounter;
    else
        isdir = [isdir i];
    end
end
features(:, isdir) = [];
images(isdir) = [];

% Match one image with all others
% Image to query:
queryImgID = 1;
queryImgNdx = cell2mat(features(1,:))==queryImgID;
queryImg = features(:, queryImgNdx);
frame_a = queryImg{2}; descrip_a = queryImg{3};

querySearchImgs = features(:, ~queryImgNdx);
matches_scores = cell(size(querySearchImgs)); % Row 2 = matches, row 3 = scores

for i = 1:length(querySearchImgs)
    descrip_b = querySearchImgs{3, i};
    [matches_scores{2,i}, matches_scores{3, i}] = vl_ubcmatch(...
        descrip_a, descrip_b);
    matches_scores{1,i} = querySearchImgs{1,i}; % ID of image being matched
end

% Get top match
maxcell = cellfun(@(x) max(x(:)), matches_scores(3,:));
[best_match, ndx] = max(maxcell);
bestMatchID = matches_scores{1, ndx};
bestMatchNdx = find(cell2mat(features(1,:))==bestMatchID);
queryImgNdx = find(queryImgNdx);
sift_mosaic(images{queryImgNdx}, images{bestMatchNdx});