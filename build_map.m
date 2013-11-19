function image_map = build_map(mosaic, im_ref)

image_map = im2double(im_ref); % mosaic.reference{1};
origin = [size(im_ref, 1), 1]; % mosaic.origins{1};
num_matches = numel(mosaic.pieces);

for i = 1:num_matches
    map_size = size(image_map);
    space_top = origin(1);
    space_bottom = map_size(1) - space_top;
    space_left = origin(2);
    space_right = map_size(2) - space_left;

    image_map_tmp = mosaic.pieces{i};
    origin_tmp = mosaic.origins{i};
    
    map_size_tmp = size(image_map_tmp);
    space_top_tmp = origin_tmp(1);
    space_bottom_tmp = map_size_tmp(1) - space_top_tmp;
    space_left_tmp = origin_tmp(2);
    space_right_tmp = map_size_tmp(2) - space_left_tmp;

    % Adjust images so that they have the same origin
    if space_top > space_top_tmp
        toprows = space_top - space_top_tmp;
        image_map_tmp = [nan(toprows, map_size_tmp(2), 3); image_map_tmp];
        map_size_tmp = size(image_map_tmp);
    end
    if space_top_tmp > space_top
        toprows = space_top_tmp - space_top;
        image_map = [nan(toprows, map_size(2), 3); image_map];
        map_size = size(image_map);
        origin(1) = origin(1) + toprows;
    end
    
    if space_bottom > space_bottom_tmp
        bottomrows = space_bottom - space_bottom_tmp;
        image_map_tmp = [image_map_tmp; nan(bottomrows, map_size_tmp(2), 3)];
        map_size_tmp = size(image_map_tmp);
    end
    if space_bottom_tmp > space_bottom
        bottomrows = space_bottom_tmp - space_bottom;
        image_map = [image_map; nan(bottomrows, map_size(2), 3)];
        map_size = size(image_map);
    end 
    
    if space_left > space_left_tmp
        colsleft = space_left - space_left_tmp;
        image_map_tmp = [nan(map_size_tmp(1), colsleft, 3), image_map_tmp];
        map_size_tmp = size(image_map_tmp);
    end
    if space_left_tmp > space_left
        colsleft = space_left_tmp - space_left;
        image_map = [nan(map_size(1), colsleft, 3), image_map];
        map_size = size(image_map);
        origin(2) = origin(2) + colsleft;
    end
        
    if space_right > space_right_tmp
        colsright = space_right - space_right_tmp;
        image_map_tmp = [image_map_tmp, nan(map_size_tmp(1), colsright, 3)];
        map_size_tmp = size(image_map_tmp);
    end
    if space_right_tmp > space_right
        colsright = space_right_tmp - space_right;
        image_map = [image_map, nan(map_size(1), colsright, 3)];
        map_size = size(image_map);
    end
    
    mass = ~isnan(image_map) + ~isnan(image_map_tmp);
    image_map(isnan(image_map)) = 0 ;
    image_map_tmp(isnan(image_map_tmp)) = 0 ;
    image_map = (image_map + image_map_tmp) ./ mass;
    image_map(image_map == 0) = NaN;
end

end