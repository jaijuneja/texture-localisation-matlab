function overlap = get_overlap(verts1, verts2, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 24/03/2014
% -------------------------------------------------------------------------
%
% GET_OVERLAP
% overlap = get_overlap(verts1, verts2, opts.numSteps, 'plotOverlap',
% valPlotOverlap)
%
% Calculates the overlap (or intersection) of two polygons as a percentage
% of their union. For two polygons A and B, it computes: (A n B)/(A u B)
% Thus, for two polygons that are completely coincident, the overlap is 1.
% For two polygons that do not intersect, the overlap is 0.
%
% Inputs:
%   - verts1, verts2:   Vertices of the two polygons in the form:
%                       [x1 x2 x3 x4 ...; y1 y2 y3 y4 ...];
%
%   Optional Properties:
%       - 'plotOverlap':    Plot the overlap region of the polygons on a
%                           new figure.
%       - 'numSteps':       Number of steps between max and min values of
%                           polygon vertices. A larger number of steps
%                           increases the numerical accuracy. 500 by
%                           default.

opts.plotOverlap = false;
opts.numSteps = 500;
opts = vl_argparse(opts, varargin);

% Create a box bounding the two polygons
xmin = min([verts1(1, :) verts2(1, :)]);
xmax = max([verts1(1, :) verts2(1, :)]);

ymin = min([verts1(2, :) verts2(2, :)]);
ymax = max([verts1(2, :) verts2(2, :)]);

% Determine the step size given the number of steps
step = min((xmax-xmin)/opts.numSteps, (ymax-ymin)/opts.numSteps);

% Refine bounding box given the step size
xmin = floor(xmin/step) * step;
xmax = ceil(xmax/step) * step;
ymin = floor(ymin/step) * step;
ymax = ceil(ymax/step) * step;

% Compute a grid for the bounding box
[xs, ys] = meshgrid(xmin:step:xmax, ymin:step:ymax);
% Binary matrix containing 1 if the centre of the grid is contained within
% the polygon
inPoly1 = inpolygon(xs, ys, verts1(1,:), verts1(2,:));
inPoly2 = inpolygon(xs, ys, verts2(1,:), verts2(2,:));

% Compute overlap percentage
overlap = sum(sum((inPoly1 & inPoly2))) / sum(sum((inPoly1 | inPoly2)));

if opts.plotOverlap
    overlapPlot = (inPoly1 | inPoly2) / 0.5;
    overlapPlot = overlapPlot + (inPoly1 & inPoly2)/0.5;
    figure;
    imagesc(overlapPlot);
    axis equal
    colormap gray
    axis tight
end

end