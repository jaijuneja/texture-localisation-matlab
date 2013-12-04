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

mosaic = get_mosaic_pieces(index, cor);

%% Superimpose all images on a global map
[image_map, origin] = build_mosaic(index, mosaic, cor);
figure; imagesc(image_map); axis off, axis equal