% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
%
% Generates a mosaic of two input images - im_ref and im_new, given a
% transformation matrix H between the two. Returns the image mosaic both
% including and excluding the reference image. Call
% [mosaic, mosaic_noref, origin] = get_mosaic_piece(im_ref, im_new, H)
%
% Inputs:
%   - im_ref:   Reference image
%   - im_new:   Second image
%   - H:        3x3 transformation matrix matching im_new to im_ref
%
% Outputs:
%   - mosaic:       Mosaic of im_ref and im_new
%   - mosaic_noref: Same as mosaic but excludes the reference image, so
%                   that we only have im_new transformed into the
%                   co-ordinate frame of im_ref
%   - origin:       The origin of the output image in pixel co-ordinates,
%                   given by the bottom left corner of the reference image
%
function [mosaic, mosaic_noref, origin] = get_mosaic_piece(im_ref, im_new, H)
    box2 = [1  size(im_new,2) size(im_new,2)  1 ;
            1  1           size(im_new,1)  size(im_new,1) ;
            1  1           1            1 ] ;
    box2_ = H \ box2 ;
    box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
    box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
    ur = min([1 box2_(1,:)]):max([size(im_ref,2) box2_(1,:)]) ;
    vr = min([1 box2_(2,:)]):max([size(im_ref,1) box2_(2,:)]) ;

    [u,v] = meshgrid(ur,vr) ;
    im1_ = vl_imwbackward(im2double(im_ref),u,v) ;

    z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
    u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
    v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
    im2_ = vl_imwbackward(im2double(im_new),u_,v_) ;

    % Find the origin of the reference image
    [i, j, ~] = ind2sub(size(im1_), find(~isnan(im1_))) ;
    origin = [max(i), min(j)] ; % Origin is lower left co-ord of ref image
    mosaic_noref = im2_ ;
    
    mass = ~isnan(im1_) + ~isnan(im2_) ;
    im1_(isnan(im1_)) = 0 ;
    im2_(isnan(im2_)) = 0 ;
    mosaic = (im1_ + im2_) ./ mass ;
end