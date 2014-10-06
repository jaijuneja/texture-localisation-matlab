function plot3d_poses(model, cor, world, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 30/01/2014
% -------------------------------------------------------------------------
% PLOT3D_POSES
% plot3d_poses(model, cor, world, varargin)
%
% Plots the estimated 3D positions of cameras relative to the world frame.
% Note that cor.H_to_world needs to be constructed using get_poses before
% this function can be called.
%
% Inputs:
%   - model
%   - cor
%   - world
%
%   Optional Properties:
%       - TBD

opts.showMosaic = false;
opts.showFeatures = false;
opts.matchesOnly = true;
opts.globalFeatsOnly = false;
opts.camSize = 0.1;
opts.arrowLength = 0.2;
opts.scaleFactor = 1;
opts.queryCameras = [];
opts.onlyDrawQueryCams = false;
opts.truePosition = [];
opts.computeOverlap = false;
opts = vl_argparse(opts, varargin);

ims_mappable = find(cellfun(@(x)(~isempty(x)), cor.H_to_world));

truePos = opts.truePosition;

offset = [0; 0];
if ~isequal(opts.scaleFactor, 1) && ~opts.showMosaic
    cor = transform_world(cor, opts.scaleFactor);
elseif opts.showMosaic
    if isequal(opts.scaleFactor, 1), opts.scaleFactor = 1000; end
    % Transform world to mosaic scale
    cor = transform_world(cor, opts.scaleFactor);
    
    % Build mosaic
    cor.H_to_ref = cor.H_to_world;
    mosaic = get_mosaic_pieces(model, cor);
    image_map = build_mosaic(model, mosaic, cor);
    
    % Determine mosaic offset
    offset = plot_transformations(model, cor, 'dontPlot', true);
    
    % Plot mosaic
    mos = warp(sparse(size(image_map,1), size(image_map,2)),image_map);
    alpha = double(~isnan(image_map(:,:,1)) | ~isnan(image_map(:,:,2)) ...
        | ~isnan(image_map(:,:,3)));
    set(mos,'FaceAlpha',  'texturemap', 'AlphaDataMapping', 'none', ...
        'AlphaData', alpha);
    hold on
end

if ~isempty(opts.truePosition)
    if ~isequal(opts.scaleFactor, 1)
        T_scale = [1 0 0; 0 1 0; 0 0 1/opts.scaleFactor];
        truePos = transform_points(truePos, T_scale);
    end
    truePos(1,:) = truePos(1,:) + offset(1);
    truePos(2,:) = truePos(2,:) + offset(1);
    h = fill(truePos(1,:), truePos(2,:), 'black');
    set(h, 'FaceAlpha', 0.5);
    hold on
end
    
% Plot features
if opts.showFeatures && opts.matchesOnly
    plot_feature_matches(world, cor, 'globalFeatsOnly', opts.globalFeatsOnly, ...
        'xOffset', offset(1), 'yOffset', offset(2), 'scaleFactor', ...
        opts.scaleFactor, 'showLegend', false);
    hold on
elseif opts.showFeatures
    plot_features(world, 'xOffset', offset(1), 'yOffset', offset(2), ...
        'scaleFactor', opts.scaleFactor, 'showLegend', false);
    hold on
end

% scales = zeros(1, length(ims_mappable));
% for i = 1:length(ims_mappable)
%     img = ims_mappable(i);
%     H_wj = inv(cor.H_to_world{img});
%     H_wj = H_wj / H_wj(3,3);
%     gammaG = cor.intrinsics \ H_wj * cor.intrinsics;
%     gamma = median(svd(gammaG));
%     scales(i) = gamma;
% end
% meanscale = mean(scales);
% stdscale = std(scales);
% cor = transform_world(cor, meanscale);
% opts.scaleFactor = meanscale;

arrow_len = opts.arrowLength * opts.scaleFactor;
arrowhead_size = 0.9 * opts.scaleFactor;
camsize = opts.camSize * opts.scaleFactor;

w.orig = [0 0 0]' + [offset; 0] * opts.showMosaic;
w.x = [1 0 0]' * arrow_len;
w.y = [0 1 0]' * arrow_len;
w.z = [0 0 1]' * arrow_len;
quiver3(repmat(w.orig(1),3,1), repmat(w.orig(2),3,1), ...
    repmat(w.orig(3),3,1), w.x, w.y, w.z, ...
    'black', 'LineWidth', 2, 'MaxHeadSize', arrowhead_size);
hold on

% Create matrix of colours
ColOrd = get(gca,'ColorOrder');
[m,~] = size(ColOrd);

for i = 1:length(ims_mappable)
    ColRow = rem(i,m);
    if ColRow == 0
      ColRow = m;
    end
    % Get the color
    Col = ColOrd(ColRow,:);

    img = ims_mappable(i);
    im_info = imfinfo(model.index.names{img});
    imsize = [im_info.Width; im_info.Height];

    H_wj = inv(cor.H_to_world{img});
    H_wj = H_wj / H_wj(3,3);
    
    [R, t] = decompose_homog(H_wj, cor.intrinsics);
    
    c = draw_camera(R, t, w, camsize, imsize);

    vertices = transform_rect(cor.H_to_world{img}, imsize);

    if opts.showMosaic
        % Offset plot to align with mosaic
        [vertices, c] = apply_offset(vertices, c, offset);
    end

    % Plot image frame projected onto world
    plot(vertices(1,:), vertices(2,:), 'LineWidth', 1.5, 'Color', Col);

    if ~opts.onlyDrawQueryCams
        % Plot camera
        h = fill3(c.plotsquare(1,:),c.plotsquare(2,:),c.plotsquare(3,:), Col);
        set(h,'FaceAlpha',0.5)
        plot3(c.plotcam(1,:), c.plotcam(2,:), c.plotcam(3,:), 'black');

        quiver3(repmat(c.orig(1),3,1), repmat(c.orig(2),3,1), ...
            repmat(c.orig(3),3,1), c.plotarrows(:,1), c.plotarrows(:,2), ...
            c.plotarrows(:,3), 'Color', Col, 'LineWidth', 2, ...
            'MaxHeadSize', arrowhead_size);
    end
end

% If a set of query cameras are input under the option 'queryCameras' then
% plot them
if isstruct(opts.queryCameras)
    for i = 1:length(opts.queryCameras.R)
        qry_c = draw_camera(opts.queryCameras.R{i}, ...
            opts.queryCameras.t{i}, w, camsize, opts.queryCameras.imsize);

        % If you want to use display the reprojected query image using the
        % decomposed values of R and t rather than the homography extracted
        % from geometric verification then uncomment below:
        homog = cor.intrinsics * [opts.queryCameras.R{i}(:,1:2) opts.queryCameras.t{i}];
        homog = inv(homog); homog = homog/homog(3,3);
        qry_vertices = transform_rect(homog, ...
            opts.queryCameras.imsize);
        
%         qry_vertices = transform_rect(opts.queryCameras.H_to_world{i}, ...
%             opts.queryCameras.imsize);
        if opts.showMosaic
            % Offset plot to align with mosaic
            [qry_vertices, qry_c] = apply_offset(qry_vertices, qry_c, offset);
        end
        
        plot(qry_vertices(1,:), qry_vertices(2,:), 'LineWidth', 1.5, ...
            'Color', 'black');

        h = fill3(qry_c.plotsquare(1,:), qry_c.plotsquare(2,:), ...
            qry_c.plotsquare(3,:), 'black');
        set(h,'FaceAlpha',0.5)
        plot3(qry_c.plotcam(1,:), qry_c.plotcam(2,:), qry_c.plotcam(3,:), ...
            'black');

        quiver3(repmat(qry_c.orig(1),3,1), repmat(qry_c.orig(2),3,1), ...
            repmat(qry_c.orig(3),3,1), qry_c.plotarrows(:,1), ...
            qry_c.plotarrows(:,2), qry_c.plotarrows(:,3), 'Color', 'black', ...
            'LineWidth', 2, 'MaxHeadSize', arrowhead_size);
        
        if opts.computeOverlap
            overlap = get_overlap(truePos, qry_vertices(:,1:end-1));
            fprintf(['Query Img #' num2str(i) ' Overlap: ' ...
                num2str(overlap) '\n']);
        end
    end
end

% Reverse y-axis so that plot is aligned with image mosaic
grid on, box off
set(gca, 'YDir', 'reverse', 'ZDir', 'reverse')
title('')
axis equal, axis tight, hold off

% -----------------------------------------------
function new_vertices = transform_rect(H, imsize)
% -----------------------------------------------
% Transforms the corners of a rectangle with dimensions imsize(2) x
% imsize(1) by the transformation matrix H. Returns the new vertices.

H_top = H(1:2,:);
H_bot = H(3,:);

% Vertices of image in homogeneous co-ordinates (note that 5th vertex is
% same as first, so that when plotting them they form a closed polygon)
im_vertices =   [  1    imsize(1)   imsize(1)           1   1
                   1            1   imsize(2)   imsize(2)   1
                   1            1           1           1   1   ];
          
% Vertices of transformed image
new_vertices = (H_top * im_vertices) ./ ...
    ([H_bot * im_vertices; H_bot * im_vertices]);

% ----------------------------------------------------
function cam = draw_camera(R, t, wld, camsize, imsize)
% ----------------------------------------------------
wld.x = wld.x * camsize / max(wld.x);
wld.y = wld.y * camsize / max(wld.y);
wld.z = wld.z * camsize / max(wld.z);

cam.x = R' * wld.x - R'*t;
cam.y = R' * wld.y - R'*t;
cam.z = R' * wld.z - R'*t;
cam.orig = - R' * t;
cam.plot = [cam.x cam.orig cam.y cam.orig cam.z];

xvec = cam.x - cam.orig;
yvec = cam.y - cam.orig;
zvec = cam.z - cam.orig;

cam.plotarrows = [xvec yvec zvec]';

imsize = imsize / max(imsize);
xvec = xvec * imsize(1);
yvec = yvec * imsize(2);

cam.plotcam = [cam.orig, cam.z+xvec/2+yvec/2, cam.orig, cam.z-xvec/2+yvec/2, ...
    cam.orig, cam.z+xvec/2-yvec/2, cam.orig, cam.z-xvec/2-yvec/2, ...
    cam.z+xvec/2-yvec/2, cam.z+xvec/2+yvec/2, cam.z-xvec/2+yvec/2, ...
    cam.z-xvec/2-yvec/2];
cam.plotsquare = [cam.z-xvec/2-yvec/2, cam.z+xvec/2-yvec/2, ...
    cam.z+xvec/2+yvec/2, cam.z-xvec/2+yvec/2];

% ------------------------------------------------------------
function [vertices, cam] = apply_offset(vertices, cam, offset)
% ------------------------------------------------------------
vertices(1,:) = vertices(1,:) + offset(1);
vertices(2,:) = vertices(2,:) + offset(2);
cam.plotsquare(1,:) = cam.plotsquare(1,:) + offset(1);
cam.plotsquare(2,:) = cam.plotsquare(2,:) + offset(2);
cam.plotcam(1,:) = cam.plotcam(1,:) + offset(1);
cam.plotcam(2,:) = cam.plotcam(2,:) + offset(2);
cam.orig(1) = cam.orig(1) + offset(1);
cam.orig(2) = cam.orig(2) + offset(2);