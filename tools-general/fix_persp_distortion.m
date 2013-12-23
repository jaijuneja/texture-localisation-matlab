function [cor, H] = fix_persp_distortion(model, cor, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 17/12/2013
% -------------------------------------------------------------------------
%
% FIX_PERSP_DISTORTION
% [cor, H] = fix_persp_distortion(model, cor)
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
%   - cor:  Correspondence structure with corrected global transformations
%           cor.H_to_ref
%   - H:    Projective transformation that corrects the original image for
%           perspective distortion (3x3 matrix)

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
offsets_bef = plot_transformations(model, cor, 'dontPlot', true);
x = x - offsets_bef(1);
y = y - offsets_bef(2);

% Compute rectangle that fits around the selected points
xy_rec = rect_fit(x, y);

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

% Correct global image transformations for perspective distortion
mappable = find(cellfun(@(x)(~isempty(x)), cor.H_to_ref));
for i = mappable
    cor.H_to_ref{i} = H \ cor.H_to_ref{i};
    cor.H_to_ref{i} = cor.H_to_ref{i} * cor.H_to_ref{i}(3,3);
end

% Plot the new mosaic
newmosaic = get_mosaic_pieces(model, cor);
new_image_map = build_mosaic(model, newmosaic, cor);
figure; imagesc(new_image_map);
hold on, axis equal, axis tight

xy_rec = H \ [x; y; ones(1, 4)];
xy_rec(1:2, :) = xy_rec(1:2, :) ./ [xy_rec(3, :); xy_rec(3, :)];
xy_rec = xy_rec(1:2, :);

offsets_aft = plot_transformations(model, cor, 'dontPlot', true);
xy_rec(1, :) = xy_rec(1, :) + offsets_aft(1);
xy_rec(2, :) = xy_rec(2, :) + offsets_aft(2);

% Plot rectangle on new mosaic
plot([xy_rec(1,:) xy_rec(1,1)], [xy_rec(2,:) xy_rec(2,1)], 'r', 'LineWidth', 2)
hold off

% test_result(xy_rec);

H = inv(H);

function vertices = rect_fit(x, y)

pts = length(x);

line_seg = zeros(2, pts); % Line vectors
line_len = zeros(1, pts); % Length of each line
for i = 1:pts
    j = ndx(i+1, pts);
    line_seg(:, i) = [x(j) - x(i); y(j) - y(i)];
    line_len(i) = pdist([0 0; line_seg(:, i)']);
end

% Get lenghts of rectangle edges that should be parallel
line_len_1 = line_len([1 3]);
line_len_2 = line_len([2 4]);
% Determine which edge in each pair is shortest
[~, minline1] = min(line_len_1);
[~, minline2] = min(line_len_2);

% Get indices of each of these lines
if isequal(minline1, 2)
    minline1 = 3;
end
minline2 = minline2 * 2;

minline = [minline1 minline2];
[~, linesame] = max(line_len(minline));
maxminline = minline(linesame);

vertices = zeros(2, pts);

vertices(:, maxminline) = [x(maxminline); y(maxminline)];
vertices(:, ndx(maxminline+1, pts)) = ...
    [x(ndx(maxminline+1, pts)); y(ndx(maxminline+1, pts))];

normal = [line_seg(2, maxminline); -line_seg(1, maxminline)];
normal = normal/pdist([0 0; normal']);

i = ndx(maxminline+2, pts);
vertices(:, i) = ...
    vertices(:, ndx(i-1, pts)) + normal * dot(line_seg(:, ndx(i-1, pts)), normal);
i = ndx(maxminline+3, pts);
vertices(:, i) = ...
    vertices(:, ndx(i-1, pts)) - line_seg(:, maxminline);

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