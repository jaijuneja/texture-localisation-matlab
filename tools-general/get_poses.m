function [cor, H_w0] = get_poses(model, cor)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 17/12/2013
% -------------------------------------------------------------------------
%
% GET_POSES
% [cor, H] = get_poses(model, cor)
%
% Removes perspective distortion from an image mosaic by prompting user for
% four points that should be mapped to rectilinear co-ordinates.
%
% Inputs:
%   - model:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%   - cor:      Correspondence structure containing links between different
%               images (graph representation using an adjacency matrix).
%               Type 'help build_correspondence' for more info
%
% Outputs:
%   - cor:  The cor.H_to_world variable is now populated
%   - H_w0: Planar homography from world co-ordinates to the reference
%           image plane. This value is also assigned to cor.H_world_toref

if isempty(cor.intrinsics)
    error('No calibration matrix detected. Set it using set_intrinsics!')
end

% Plot the original mosaic
mosaic = get_mosaic_pieces(model, cor);
image_map = build_mosaic(model, mosaic, cor);
figure; imagesc(image_map);
hold on, axis equal, axis tight

% Prompt user to input 4 points that should map to rectilinear co-ordinates
numPoints = 4;
x = zeros(1, numPoints);
y = zeros(1, numPoints);

for i = 1:numPoints
    [x(i), y(i)] = ginput(1);
    plot(x(i), y(i), 'r+')
end

plot([x x(1)], [y y(1)], 'r', 'LineWidth', 2)
drawnow;

% Convert rectangle from mosaic co-ordinates to world co-ords (mosaic is in
% positive co-ordinates only, whereas world includes negative values)
offsets_bef = plot_transformations(model, cor, 'dontPlot', true, ...
    'fromframe', 'ref');
x = x - offsets_bef(1);
y = y - offsets_bef(2);

% Compute rectangle that fits around the selected points
xy_rec = to_square(x, y);
xnew = xy_rec(1,:);
ynew = xy_rec(2,:);

% Compute the matrix A such that A * h = 0, where h is the vector of
% unknowns in H
A = zeros(8, 9);
for i = 1:4
    A_tmp = [xnew(i), ynew(i), 1, 0, 0, 0, -x(i)*xnew(i), -x(i)*ynew(i), -x(i); ...
        0, 0, 0, xnew(i), ynew(i), 1, -y(i)*xnew(i), -y(i)*ynew(i), -y(i)];
    A(i*2-1:i*2, :) = A_tmp;
end

% Hence we obtain the null space of A by computing SVD
h = null(A);
H = reshape(h, 3, 3)';
% img_scale = 1000; % Scale in pixels of output mosaic image

% Calculate aspect ratio alpha
KinvH = cor.intrinsics \ H;
alpha = norm(KinvH(:,2))/norm(KinvH(:,1));
% Calculate transformation from world plane to reference plane
H_w0 = H * [1 0 0; 0 1/alpha 0; 0 0 1];
cor.H_world_toref = H_w0;

% Correct global image transformations for perspective distortion
mappable = find(cellfun(@(x)(~isempty(x)), cor.H_to_ref));
cor.H_to_world = cor.H_to_ref;
for i = mappable
    cor.H_to_world{i} = H_w0 \ cor.H_to_ref{i};
    cor.H_to_world{i} = cor.H_to_world{i} / cor.H_to_world{i}(3,3);
end

% Plot the new mosaic
cor_plot = cor;
cor_plot.H_to_ref = cor_plot.H_to_world;
% Scale it up so that it can be visualised properly on the mosaic
% cor_plot = transform_world(cor_plot, img_scale);
figure; plot_transformations(model, cor_plot);

% Plot the mosaic on the world plane
% newmosaic = get_mosaic_pieces(model, cor_plot);
% new_image_map = build_mosaic(model, newmosaic, cor_plot);
% figure; imagesc(new_image_map);
% hold on
%
% Show the rectangle transformed to the world plane
% xy_rec = transform_world(H_w0, img_scale) * [x; y; ones(1, 4)];
xy_rec = H_w0 * [x; y; ones(1, 4)];
xy_rec(1:2, :) = xy_rec(1:2, :) ./ [xy_rec(3, :); xy_rec(3, :)];
xy_rec = xy_rec(1:2, :);
%
% offsets_aft = plot_transformations(model, cor_plot, 'dontPlot', true);
% xy_rec(1, :) = xy_rec(1, :) + offsets_aft(1);
% xy_rec(2, :) = xy_rec(2, :) + offsets_aft(2);
%
% Plot rectangle on new mosaic
plot([xy_rec(1,:) xy_rec(1,1)], [xy_rec(2,:) xy_rec(2,1)], 'r', ...
    'LineWidth', 2)
% test_result(xy_rec);
hold off

function vertices = to_square(x, y)
rectangle = [x; y];
pts = length(rectangle);
vertices = zeros(2, pts);

dist_to_00 = pdist([0 0; rectangle']);
dist_to_00 = dist_to_00(1:pts);
 
[~, idx] = min(dist_to_00);

square_clockwise = [0 1 1 0; 0 0 1 1];

for i = 1:pts
    j = ndx(idx + i - 1, pts);
    vertices(:, j) = square_clockwise(:, i);
end

function index = ndx(index, numPoints)

if index > numPoints
    index = index - numPoints;
elseif index <= 0
    index = numPoints + index;
end

function dot_prod = test_result(xy_rec)
% Check to make sure that the output is in rectilinear co-ordinates (i.e.
% dot product of edges is zero).
pts = length(xy_rec);
line_seg = zeros(2, pts); % Line vectors
for i = 1:pts
    j = ndx(i+1, pts);
    line_seg(:, i) = [xy_rec(1,j) - xy_rec(1,i); xy_rec(2,j) - xy_rec(2,i)];
end

dot_prod = zeros(1, pts);
for i = 1:pts
    j = ndx(i+1, pts);
    dot_prod(i) = dot(line_seg(:, i), line_seg(:,j));
end
dot_sum = sum(abs(dot_prod));

assert(dot_sum < 1e-5);