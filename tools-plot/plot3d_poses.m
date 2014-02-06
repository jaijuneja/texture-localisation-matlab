function plot3d_poses(model, cor, world, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 08/12/2013
% -------------------------------------------------------------------------
% PLOT3D_POSES
%

opts.showMosaic = false;
opts.showFeatures = false;
opts.matchesOnly = true;
opts.globalFeatsOnly = false;
opts.camSize = 0.1;
opts.arrowLength = 0.2;
opts.scaleFactor = 1;
opts = vl_argparse(opts, varargin);

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

ims_mappable = find(cellfun(@(x)(~isempty(x)), cor.H_to_world));

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
    
    r1r2t = cor.intrinsics \ H_wj;

    r1 = r1r2t(:,1);
    scale1 = norm(r1);
    r2 = r1r2t(:,2);
    scale2 = norm(r2);
    % r2 = r2/norm(r2);
    r3 = cross(r1, r2);
    % r3 = r3/norm(r3);
    % r1 = cross(r2, r3);

    R = [r1 r2 r3];
    % Improvement to commented out method above - see planar-struct3.pdf
    [U, ~, V] = svd(R);
    R = U * V';
    
    % Want to scale t so that it represents r1r2t as closely as possible
    scalet = mean([scale1 scale2]);
    % Alternative scaling
    % scalet2 = norm(r1r2t(:,1:2))/norm(R(:,1:2));
    t = r1r2t(:,3);
    t = t/scalet;
    
    c = draw_camera(R, t, w, camsize, imsize);

    vertices = transform_rect(cor.H_to_world{img}, imsize);
    
    if opts.showMosaic
        % Offset plot to align with mosaic
        [vertices, c] = apply_offset(vertices, c, offset);
    end
    
    plot(vertices(1,:), vertices(2,:), 'LineWidth', 1.5, 'Color', Col);
    
    h = fill3(c.plotsquare(1,:),c.plotsquare(2,:),c.plotsquare(3,:), Col);
    set(h,'FaceAlpha',0.5)
    plot3(c.plotcam(1,:), c.plotcam(2,:), c.plotcam(3,:), 'black');
    
    quiver3(repmat(c.orig(1),3,1), repmat(c.orig(2),3,1), ...
        repmat(c.orig(3),3,1), c.plotarrows(:,1), c.plotarrows(:,2), ...
        c.plotarrows(:,3), 'Color', Col, 'LineWidth', 2, ...
        'MaxHeadSize', arrowhead_size);
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