%% Get paths to all images in database and generate index struct

[index, images, path] = build_index('test_images/concrete_resized/', 'numWords', 5e4);

%% Find correspondences between images

if ~exist(path.cor, 'file')
    % Build correlation structure
    cor = build_correspondence(index, 'percentThresh', 0.4, 'numThresh', 20);
    % Save it
    save(path.cor, 'cor');
else
    load(path.cor);
end

view(cor.graph);

%% Browse correspondence structure
plot_images(index, cor);

%% Change reference image
cor = set_refimg(cor, 9);
save(path.cor, 'cor');
%% Get global camera poses
cor = get_poses(index, cor);
figure; plot3d_poses(index, cor);
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

%% Superimpose image mosaic with feature plot
figure;
plot_everything(index, world, cor, 'showFeatures', true, ...
    'matchesOnly', true, 'showImgBorders', false, 'showMosaic', true)

%% Correct image for perspective distortion
cor = fix_persp_distortion(index, cor);