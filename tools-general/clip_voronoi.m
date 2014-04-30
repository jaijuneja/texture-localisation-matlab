function [DT, numNewCtrs] = clip_voronoi(DT, region)
region = num2cell(region);
[xmin xmax ymin ymax] = region{:};
centers = DT.Points;
newCenters = [];

[V, R] = voronoiDiagram(DT);

toLeft = find(V(:,1) < xmin | V(:,1)==Inf);
for i = toLeft'
    ctrsToLeft = find(cellfun(@(x)any(x==i), R));
    for j = ctrsToLeft
        % Reflect centers across left border
        newCenters = vertcat(newCenters, [2*xmin - centers(j, 1), centers(j, 2)]);
    end
end

toRight = find(V(:,1) > xmax | V(:,1)==Inf);
for i = toRight'
    ctrsToRight = find(cellfun(@(x)any(x==i), R));
    for j = ctrsToRight
        % Reflect centers across left border
        newCenters = vertcat(newCenters, [2*xmax - centers(j, 1), centers(j, 2)]);
    end
end

above = find(V(:,2) > ymax | V(:,2)==Inf);
for i = above'
    ctrsAbove = find(cellfun(@(x)any(x==i), R));
    for j = ctrsAbove
        % Reflect centers across left border
        newCenters = vertcat(newCenters, [centers(j, 1), 2*ymax - centers(j, 2)]);
    end
end

below = find(V(:,2) < ymin | V(:,2)==Inf);
for i = below'
    ctrsBelow = find(cellfun(@(x)any(x==i), R));
    for j = ctrsBelow
        % Reflect centers across left border
        newCenters = vertcat(newCenters, [centers(j, 1), 2*ymin - centers(j, 2)]);
    end
end

% Remove duplicates
newCenters = unique(newCenters, 'rows');
numNewCtrs = size(newCenters, 1);

% New centers for Delaunay triangulation
newCenters = [newCenters; centers];

DT = delaunayTriangulation(newCenters);
end