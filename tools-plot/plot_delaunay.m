function plot_delaunay(DT)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 27/03/2014
% -------------------------------------------------------------------------
%
% PLOT_DELAUNAY
%
% Plot the triangles from the Delaunay triangulation DT. There is an
% in-built MATLAB function called 'triplot' that already achieves this,
% however it only operates on objects of the class 'delaunayTriangulation',
% which are read-only. This function has been created for cases where the
% output DT from the function 'delaunayTriangulation' has been modified and
% needs to be plotted.
%
% Inputs:
%   - DT:   Struct which contains the same fields as a
%           delaunayTriangulation object

for i = 1:size(DT.ConnectivityList, 1)
    pts = DT.Points(DT.ConnectivityList(i,:), :);
    pts(end+1,:) = pts(1,:);
    plot(pts(:,1), pts(:,2), 'LineWidth', 2, 'Color', 'black');
    hold on
end

end