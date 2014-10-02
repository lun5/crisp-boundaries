%% testing the effect of opts.sig on the joint distribution and calculation
% of PMI
% Luong Nguyen 9/25/14

%% 
addpath(genpath(pwd));
%I = imread('test_images/253027.jpg'); % zebra
%I = imread('test_images/syntheticImage.png'); % synthetic image
I = imread('test_images/tp10-611_22528_16384_2048_2048.tif'); % H&E
imshow(I);
rect = getrect;
I = imcrop(I,rect);

opts = setEnvironment('speedy');
I = im2uint8(I);
if (size(I,3)==1)
    I = repmat(I,[1 1 3]);
end
figure; imshow(I);
% calculate features
num_scales = opts.num_scales;
scale_offset = opts.scale_offset;
f_maps = getFeatures(double(I)/255,num_scales+scale_offset,opts.features.which_features,opts);
Nsamples = opts.kde.Nkernels;
%F = sampleF(f_maps,Nsamples,opts);
% vary opts.sig to observe effects on joint distribution and pmi
sigvec = [0.1,0.25,0.5,1, 2];
pairStatsvec = zeros(length(sigvec),5);
% parameters for ploting joint density and pmi
tol = opts.kde.kdtree_tol;
x = 0:0.01:1;
y = 0:0.01:1;
[X,Y] = meshgrid(x,y);
Fim = [X(:),Y(:)];

for i = 1:length(sigvec)
    opts.sig = sigvec(i);
    [F, pairStats] = sampleF_withStats(f_maps,Nsamples,opts);
    pairStatsvec(i,:) = pairStats;
    Fsym = [F; [F(:,2) F(:,1)]]; % symmetric F(A,B) = F(B,A). 
    p = kde(Fsym',0.05,[],'e');
    %% calculate joint density
    pd = evaluate(p,Fim',tol);
    reg = opts.p_reg;
    pJoint = reg + pd;
    pJoint_mesh = reshape(pJoint, size(X));
    figure;contourf(x,y,log(pJoint_mesh),30); axis square; colorbar;
    xlabel('Luminance A'); ylabel('Luminance B');
    set(gca,'XTick',0:0.1:1)
    set(gca,'YTick',0:0.1:1)
    %% evaluate p(A)p(B)
    N = floor(size(Fim,2)/2); assert((round(N)-N)==0);
    p2_1 = marginal(p,1:N);
    p2_2 = marginal(p,N+1:(2*N));
    p2 = joinTrees(p2_1,p2_2,0.5);
    pMarg_x = evaluate(p2,X(:)',tol);
    pMarg_y = evaluate(p2,Y(:)',tol);
    pProd = pMarg_x.*pMarg_y +reg;
    %% calculate pmi
    pmi = log((pJoint.^(opts.joint_exponent))./pProd);
    pmi_mesh = reshape(pmi,size(X));
    figure;contourf(x,y,pmi_mesh,20);
    axis square; colorbar; 
    xlabel('Luminance A'); ylabel('Luminance B');
    set(gca,'XTick',0:0.1:1)
    set(gca,'YTick',0:0.1:1)
end

close all;

% datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
% I = imread(fullfile(datadir,'tp10-867-1_47104_22528_2048_2048.tif'));
% imshow(I);
% 'tp10-611_22528_14336_2048_2048.tif'
% 'tp10-611_53248_12288_2048_2048.tif'
% 'tp10-834-2_36864_12288_2048_2048.tif'
% 'mws09-778a_20480_22528_2048_2048.tif'
% 'tp10-762-1_28672_20480_2048_2048.tif'
% 'tp10-611_59392_26624_2048_2048.tif'
% tp10-611_67584_18432_2048_2048
% tp10-611_22528_16384_2048_2048