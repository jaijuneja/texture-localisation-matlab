function plot_images(model, cor)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 09/12/2013
% -------------------------------------------------------------------------
%
% PLOT_IMAGES
% plot_images(model, cor)
%
% Plots all of the images in the database. When an image is clicked on, a
% figure pops up showing the selected image and all of its matches, as well
% as their matching scores. A matched image  can then be clicked on to show
% an interactive plot of feature-to-feature matches between the two images.
% This allows the user to conveniently browse the correspondence structure
% ('cor') and visualise the matches between images.
%
% Inputs:
%   - model:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%   - cor:      Correspondence structure containing links between different
%               images (represented by a graph). Type 'help
%               build_correspondence' for more info

clf;
set(gcf,'name','Images in Database','numbertitle','off')

num_ims = length(model.index.ids);

for i = 1:num_ims
    vl_tightsubplot(num_ims, i, 'Margin', 1e-2);
    im = imread(model.index.names{i});
    data.h(i) = imagesc(im); % Return image handle
    title(['Image ' num2str(i)])
    set(data.h(i), 'ButtonDownFcn', @zoomIntoImage);
    set(gca, 'XTick', [], 'YTick', []);
    axis equal, axis tight
end

% Data for interactive plot
data.cor = cor;
data.model = model;
guidata(gcf, data);

function zoomIntoImage(h, event, data)

data = guidata(h);

% Determine which image was selected and find its matches
im_ndx = find(h == data.h) ;
im_matches = data.cor.img_matches{im_ndx};
num_matches = length(im_matches);

clf;
set(gcf,'name',['Image ' num2str(im_ndx) ' and Its Matches'], ...
    'numbertitle','off')

% Plot data
numcols = 3;
main_imspace = numcols^2;
numrows = ceil((main_imspace+num_matches+1)/numcols);

% Plot reference image (that was clicked on)
subplot(numrows, numcols, 1:main_imspace)
im = imread(data.model.index.names{im_ndx}) ;
imagesc(im);
title(['Image ' num2str(im_ndx) ' (click on a match)'])
set(gca, 'XTick', [], 'YTick', []);
axis equal, axis tight

% Initialise match info (for interactive plot)
match.ref_img = im;
match.ref_ndx = im_ndx;
match.ref_frames = data.model.index.frames{im_ndx};
match.H = cell(1, num_matches);
match.frames = cell(1, num_matches);
match.img = cell(1, num_matches);
match.img_ndx = zeros(1, num_matches);

% Plot each matched image and get info about each match
for i = 1:num_matches
    match_ndx = data.cor.img_matches{im_ndx}(i);
    im2 = imread(data.model.index.names{match_ndx});
    
    % Plot matched image
    subplot(numrows, numcols, main_imspace+i);
    match.h(i) = imagesc(im2);
    title(['Match #' num2str(i) ':' char(10) 'Image ' num2str(match_ndx)])
    xlabel(['Score: ' num2str(data.cor.scores{im_ndx}(i))])
    set(match.h(i), 'ButtonDownFcn', @zoomIntoMatch);
    set(gca, 'XTick', [], 'YTick', []);
    axis equal, axis tight
    
    % Get matched image data
    match.H{i} = data.cor.H{im_ndx}{i};
    match.frames{i} = data.model.index.frames{match_ndx};
    match.feature_matches{i} = data.cor.feature_matches{im_ndx}{i};
    match.img{i} = im2;
    match.img_ndx(i) = match_ndx;
end

% One extra subplot for the back arrow
subplot(numrows, numcols, main_imspace+num_matches+1)
title(['Click arrow' char(10) 'to go back'])
arrow = text(0,0.5,'\leftarrow','FontSize',50);
set(arrow, 'ButtonDownFcn', @goBack);
set(gca, 'color', 'none')
axis equal, axis tight, axis off

% Need cor and model in case the back arrow is clicked
match.cor = data.cor;
match.model = data.model;
guidata(gcf, match);


function zoomIntoMatch(h, event, match)

match = guidata(h);

% Determine which match was selected
ndx = find(h == match.h) ;

figure; clf;
set(gcf,'name',['Matches Between Image ' num2str(match.ref_ndx) ...
    ' and Image ' num2str(match.img_ndx(ndx))],'numbertitle','off')

% Plot the matches
plotMatches(match.ref_img, match.img{ndx}, match.ref_frames, ...
    match.frames{ndx}, match.feature_matches{ndx}, 'homography', match.H{ndx});

function goBack(h, event, match)
match = guidata(h);
% Return to main plot
plot_images(match.model, match.cor)