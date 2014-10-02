%% function [f_maps] = getFeatures(im_rgb,scale,which_feature,opts)
% 
% INPUTS
%  im_rgb        - NxMx3 query image; C can be 1 (grayscale) or 3 (rgb color)
%  scale         - how many times should the image be downsampled?
%  which_feature - which feature types to compute? can contain multiple entries
%  opts          - parameter settings (see setEnvironment)
%
% OUTPUTS
%  f_maps        - NxMxF array of F feature maps input im_rgb
% 
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2014 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

%% temporary get features for the test image
matfiledir = 'C:\Users\luong_nguyen\Dropbox\DDIresearch\crisp-boundaries';
source_svs = 'tp10-867-1';
%source_svs = 'tp10-611';
purple_source = load(fullfile(matfiledir,[source_svs 'training_purple.mat']));
training_data_purple_source =purple_source.training_data_purple;
pink_source = load(fullfile(matfiledir,[source_svs 'training_pink.mat']));
training_data_pink_source = pink_source.training_data_pink;

%% get the rotation matrix 
% source image
source_training_data = [training_data_purple_source(:,1:2000) training_data_pink_source(:,1:6000)];
[U,~,~] = svd(source_training_data,0);
source_rotation_matrix = [-U(:,1) U(:,2:3)]';
r = I(:,:,1); g = I(:,:,2); b = I(:,:,3);
rgb_coords = [r(:)'; g(:)'; b(:)'];
rotated_coordinates = source_rotation_matrix*double(rgb_coords);
theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
f_maps = reshape(theta,size(r));