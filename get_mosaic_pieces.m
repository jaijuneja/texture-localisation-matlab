function mosaic = get_mosaic_pieces(model, cor)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
% 
% GET_MOSAIC_PIECES
% mosaic = get_mosaic_pieces(model, cor)
%
% Do something
%
% Inputs:
%   - model:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%   - cor:      Correspondence structure containing links between different
%               images (graph representation using an adjacency matrix).
%               Type 'help build_correspondence' for more info
%
% Outputs:
%   - mosaic:   Structure containing individual pieces of the mosaic, which
%               are obtained by calling stitch_images
%                   *   mosaic.pieces is a 1xn cell array for n pieces in
%                       the where each element contains an image
%                       transformed relative to the reference image. The
%                       ref image has ID cor.ref_img
%                   *   mosaic.origins is a 1xn cell array where, for each
%                       piece of the mosaic, the origin of the reference
%                       image is given in pixel co-ordinates ([row col])

% Find images which can be mapped to global frame (excluding reference img)
ims_mappable = find(cellfun(@(x)(~isempty(x) & ~isequal(x, eye(3))), ...
    cor.H_to_ref));

num_pieces = length(ims_mappable);

im_ref = imread(model.index.names{cor.ref_img});

% Initialise loop variables
mosaic.pieces = cell(1, num_pieces);
mosaic.origins = cell(1, num_pieces);

% Get all the pieces of the mosaic and put them in the mosaic structure
for match_ndx = 1:num_pieces
    im_new_ndx = ims_mappable(match_ndx);
    im_new = imread(model.index.names{im_new_ndx});
    H = cor.H_to_ref{im_new_ndx};
    [~, mosaic.pieces{match_ndx}, mosaic.origins{match_ndx}] ...
        = stitch_images(im_ref, im_new, H);
end

end