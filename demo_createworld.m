%% Get paths to all images in database and generate index struct

[index, images, path] = build_index('test_images/adjacency_demo/');

%% Find correspondences between images

if ~exist(path.cor, 'file')
    % Build correlation structure
    cor = build_correspondence(index, 'refImg', 1);
    % Save it
    save(path.cor, 'cor');
else
    load(path.cor);
end

view(cor.graph);

%% Build global map of features

if ~exist(path.world, 'file')
    % Build correlation structure
    world = build_world(index, cor);
    % Save it
    save(path.world, 'world');
else
    load(path.world);
end

figure; plot_features(world);

%% Plot mosaic of corresponding images
% Need to put this all in get_mosaic_piece with the function
% get_mosaic_piece at the bottom of the file, as it is only called within
% the file itself

mosaic = get_mosaic_pieces(index, cor);

%% Superimpose all images on a global map
im_ref = imread(index.index.names{cor.ref_img});
[image_map, origin] = build_mosaic(mosaic, im_ref);
figure;
imagesc(image_map); axis off, axis equal