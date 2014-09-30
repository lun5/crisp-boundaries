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
    % hue
    theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
    im_theta = reshape(theta,size(r));
    figure; imagesc(im_theta); axis equal; axis([0 size(im_theta,2) 0 size(im_theta,1)]);
    colormap(hsv); colorbar;
    % saturation
    sat = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
    im_sat = reshape(sat,size(r));
    figure; imagesc(im_sat); axis equal; axis([0 size(im_sat,2) 0 size(im_sat,1)]);
    colormap(jet); colorbar; 
    % brightness
    brightness = rotated_coordinates(1,:);
    im_brightness = reshape(brightness,size(r));
    figure; imagesc(im_brightness); axis equal; axis([0 size(im_brightness,2) 0 size(im_brightness,1)]);
    colormap(jet); colorbar; 
    
    if strcmp(which_feature,'hue opp')
        im = im_theta;
    end
    %%
    if (strcmp(which_feature,'saturation opp'))       
        im = im_sat; 
    end
    %%
    if (strcmp(which_feature,'brightness opp'))
        im = im_brightness;
    end
    %% 
    if (strcmp(which_feature,'hsb oppCol'))
        im = cat(3,im_theta, im_sat, imbrightness);
    end
        
    %% downsample
    im = imresize(im,2^(-(scale-1)));
        
    %%
    f_maps = im;
end