function cor = set_refimg(cor, refImg)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 27/11/2013
% -------------------------------------------------------------------------
%
% SET_REFIMG
% cor = set_refimg(cor, refImg)
%
% Set the reference image and calculate the transformation of all other
% images to the reference frame. Transformations are saved to cor.H_to_ref.
%
% Inputs:
%   - cor:      Correspondence structure containing links between different
%               images (type 'help build_correspondence' for more info)
%   - refImg:   ID of the reference image
%
% Outputs:
%   - cor

cor.ref_img = refImg;

% Delete all old values in cor.H_to_ref
cor.H_to_ref = cell(size(cor.H_to_ref));

% Compute minimum spanning tree to traverse graph
cor.img_order = graphtraverse(sparse(cor.adjacency), refImg, ...
    'Method', 'BFS', 'Directed', false);

cor.H_to_ref{refImg} = eye(3); % First image acts as global co-ordinate frame
previous_node = refImg; % Start at ref img and then traverse minimum spanning tree
last_node_branched_ndx = 1;
for i = 2:length(cor.img_order)
    current_node = cor.img_order(i);
    % If current node is not connected to previous node, then update
    % previous_node
    if isempty(find(cor.img_matches{current_node} == previous_node, 1))
        node_branched_ndx = find(ismember(cor.img_order(last_node_branched_ndx+1:i-1), ...
            cor.img_matches{current_node}));
        last_node_branched_ndx = last_node_branched_ndx + node_branched_ndx(1);
        previous_node = cor.img_order(last_node_branched_ndx);
    end
    H_previous = cor.H_to_ref{previous_node};
    previous_node_ndx = cor.img_matches{current_node} == previous_node;
    H_current_to_previous = cor.H{current_node}{previous_node_ndx};
    H = H_previous / H_current_to_previous;
    cor.H_to_ref{current_node} = H / H(3,3);
end

end