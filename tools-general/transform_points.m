function new_points = transform_points(points, H)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 04/02/2014
% -------------------------------------------------------------------------
%
% TRANSFORM_POINTS
% new_points = transform_points(points, H)
%
% Transforms a set of points using the transformation H.
%
% Inputs:
%   - points:   2xn matrix for n points, of the form:
%               [[x1; y1] [x2; y2] ... [xn; yn]]
%   - H:        3x3 transformation matrix
%
% Outputs:
%   - new_points:   Matrix of transformed points; has the same dimensions
%                   as input argument 'points'

new_points = H * [points; ones(1, size(points, 2))];

new_points = new_points(1:2,:) ./ [new_points(3,:); new_points(3,:)];

end