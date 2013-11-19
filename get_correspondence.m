function cor = get_correspondence(model, percent_thresh)
% Maybe instead have the threshold as a percentage of the best match score
% AND greater than a certain value
if nargin == 1, percent_thresh = 0.4; end
abs_thresh = 20;

cor.id = cell(1, length(model.index.ids));
cor.img_matches = cell(1, length(model.index.ids));
cor.feature_matches = cell(1, length(model.index.ids));
cor.scores = cell(1, length(model.index.ids));
cor.H = cell(1, length(model.index.ids));
cor.H_to_ref = cell(1, length(model.index.ids)); % Transformation to reference image

cor.adjacency = zeros(length(model.index.ids));

for i = model.index.ids
    [ids, scores, result] = visualindex_query(model, i);
    score_thresh = scores(2) * percent_thresh; % 60% of best score by default
    goodmatch = find(scores > score_thresh & scores > abs_thresh ...
        & ids ~= result.query);
    cor.id{i} = result.query;
    cor.img_matches{i} = ids(goodmatch);
    cor.scores{i} = scores(goodmatch);
    cor.H{i} = result.H(goodmatch);
    cor.feature_matches{i} = result.matches(goodmatch);
    cor.adjacency(i, ids(goodmatch)) = 1;
end

% Construct cell array of image names
node_names = cell(1, numel(model.index.ids));
node_names(:) = {'Image '};
node_ids = strtrim(cellstr(num2str(model.index.ids'))');
node_names(:) = strcat(node_names(:), node_ids(:));
% Construct biograph of images using adjacency matrix
cor.graph = biograph(cor.adjacency, node_names);
% Use view(graph) to display biograph
% Use imagesc(adjacency), axis equal to display adjacency matrix

% Force symmetry of adjacency matrix
adjacency_symmetric = cor.adjacency & cor.adjacency';
% Compute minimum spanning tree to traverse graph
order = graphtraverse(sparse(adjacency_symmetric), 1, 'Method', 'BFS', 'Directed', false);

cor.H_to_ref{1} = eye(3); % First image acts as global co-ordinate frame
previous_node = 1; % Start at image 1 and then traverse minimum spanning tree
last_node_branched_ndx = 1;
for i = 2:length(order)
    current_node = order(i);
    % If current node is not connected to previous node, then update
    % previous_node
    if isempty(find(cor.img_matches{current_node} == previous_node))
        node_branched_ndx = find(ismember(order(last_node_branched_ndx+1:i-1), ...
            cor.img_matches{current_node}));
        last_node_branched_ndx = last_node_branched_ndx + node_branched_ndx(1);
        previous_node = order(last_node_branched_ndx);
    end
    H_previous = cor.H_to_ref{previous_node};
    previous_node_ndx = find(cor.img_matches{current_node} == previous_node);
    H_current_to_previous = cor.H{current_node}{previous_node_ndx};
    H = H_current_to_previous * H_previous;
    cor.H_to_ref{current_node} = H / H(3,3);
end

end