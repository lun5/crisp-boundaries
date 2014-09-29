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

function [f_maps] = getFeatures_oppCol(im_rgb,scale,which_feature,opts)
    
    im = [];
    r = im_rgb(:,:,1); g = im_rgb(:,:,2); b = im_rgb(:,:,3);
    rotated_coordinates = opts.features.rotation_matrix*double([r(:)'; g(:)'; b(:)']);
    
    if strcmp(which_feature,'hue_opp')
        theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
        im = reshape(theta,size(r));
    end
    %%
    if (strcmp(which_feature,'saturation_opp'))
        sat = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
        im = reshape(sat,size(r));
    end
    %%
    if (strcmp(which_feature,'brightness_opp'))
        brightness = rotated_coordinates(1,:);
        im = reshape(brightness,size(r));
    end
    %% Do I need to account for an option to have all of the features?
    
    %% downsample
    im = imresize(im,2^(-(scale-1)));
        
    %%
    f_maps = im;
end