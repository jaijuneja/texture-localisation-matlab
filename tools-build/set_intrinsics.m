function cor = set_intrinsics(cor, focal_length, principal_pt, skew)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 28/01/2014
% -------------------------------------------------------------------------
%
% SET_INTRINSICS
% cor = set_intrinsics(cor, focal_length, principal_pt, skew)
%
% Generate a correspondence structure using the index of images.
%
% Inputs:
%   - cor:              Connectivity structure (see build_correspondence)
%   - focal_length:     2-vector [fx fy] containing x and y focal lengths 
%                       of camera
%   - principal_pt:     2-vector [px py] containing x and y principal point
%                       offset
%   - skew:             Non-zero when pixels aren't square, zero by default
%
% Outputs:
%   - cor:  Correspondence structure (see build_correspondence), which now
%           includes the intrinsic calibration matrix of camera. Takes the
%           form:
%           [ fx s  px
%             0  fy py
%             0  0  1 ]

if nargin < 4
    skew = 0;
end

cor.intrinsics = [ 	focal_length(1)	skew            principal_pt(1)
                   	0           	focal_length(2) principal_pt(2)
                   	0             	0               1               ];
end