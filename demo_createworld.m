%% Initialise: get paths to all images in database

opts.image_path = 'test_images/adjacency_demo/'; % Make sure this ends with /
opts.etc_path = strcat(opts.image_path, 'etc/');
opts.index_path = strcat(opts.etc_path, 'index.mat');
opts.cor_path = strcat(opts.etc_path, 'cor.mat');
opts.world_path = strcat(opts.etc_path, 'world.mat');
opts.numWords = round(5e4 * (6/16));

% Obtain listing of all files in image path
image_files = dir(opts.image_path);

% Initiate blank cell array of image paths
images = cell(1, length(image_files));
isdir = [];

% Populate cell array with paths to images
for i = 1:length(image_files)
    if ~image_files(i).isdir
        images(i) = { strcat(opts.image_path, image_files(i).name) };
    else
        isdir = [isdir i]; % If the listed file is a directory, record its index 
    end
end
images(isdir) = []; % Remove all entries that are directories (not images)

%% Iinitialize the index structure to estimate visual word vocabulary

% Save the generated index in the specified path
if ~exist(opts.etc_path, 'dir')
    mkdir(opts.etc_path);
end
% If index doesn't already exist, generate it
if ~exist(opts.index_path, 'file')
    index = visualindex_build(images, 'numWords', opts.numWords);
    % Add images to index
    index = visualindex_add(index, images, 1:length(images));
    % Save index
    save(opts.index_path, 'index');
else
    load(opts.index_path);
end

%% Find correspondences between images

if ~exist(opts.cor_path, 'file')
    % Build correlation structure
    cor = get_correspondence(index);
    % Save it
    save(opts.cor_path, 'cor');
else
    load(opts.cor_path);
end

%% Build global map of features

if ~exist(opts.world_path, 'file')
    % Build correlation structure
    world = build_world(index, cor);
    % Save it
    save(opts.world_path, 'world');
else
    load(opts.world_path);
end

plot_features(world);

%% Plot mosaic of corresponding images
% Take first image (ID = 1) as an example
% find images which can be mapped to global frame (excluding reference img)
im_ref_id = 1;
ims_mappable = find(cellfun(@(x)(~isempty(x) & ~isequal(x, eye(3))), cor.H_to_ref));
im_ref = imread(index.index.names{im_ref_id});

num_pieces = length(ims_mappable);

im_mosaic = cell(1, num_pieces);
mosaic.pieces = cell(1, num_pieces);
mosaic.origins = cell(1, num_pieces);
% Need co-ordinates of bottom left

H_current = eye(3); % Homography to reference frame is initially identity
figure;
for match_ndx = 1:num_pieces
    im_new_ndx = ims_mappable(match_ndx);
    im_new = imread(index.index.names{im_new_ndx});
    H = cor.H_to_ref{im_new_ndx};
    [im_mosaic{match_ndx}, ...
        mosaic.pieces{match_ndx}, ...
        mosaic.origins{match_ndx}] ...
        = get_mosaic_piece(im_ref, im_new, H);
    
    subplot(1, num_pieces, match_ndx), subimage(im_mosaic{match_ndx})
end

%% Superimpose all images on a global map
[image_map, origin] = build_mosaic(mosaic, im_ref);
figure;
imagesc(image_map); axis off, axis equal