%% 
addpath(genpath(pwd));
% I = imread('test_images/253027.jpg'); % zebra
% I = imread('test_images/syntheticImage.png'); % synthetic image
% I = imread('test_images/synthetic3.png');
% I = imread('test_images/tp09-96_20480_10240_2048_2048.tif'); % H&E
datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
if ~ exist(datadir,'dir')
    datadir = '/Users/lun5/Research/color_deconvolution/TissueImages/';
end
I = imread(fullfile(datadir,'tp10-867-1_4096_20480_2048_2048.tif'));
% I = imread(fullfile(datadir,'tp10-867-1_47104_22528_2048_2048.tif'));
% I = imread(fullfile(datadir,'tp10-867-1_26624_24576_2048_2048.tif'));
% I = imread(fullfile(datadir,'tp10-867-1_34816_18432_2048_2048.tif'));
% I = imread(fullfile(datadir,'tp10-611_22528_16384_2048_2048.tif'));
imshow(I);
% rect = getrect;
% I_purple = imcrop(I,round(rect));
% rect = getrect;
% I_pink = imcrop(I,round(rect));
% addpath('../color-deconvolution/');
% training_data = cat(2, raw2rgb(I_purple), raw2rgb(I_pink));
% [U,~,~] = svd(training_data + normrnd(0,1),0);
% rotation_matrix = [-U(:,1) U(:,2:3)]';
opts = setEnvironment('speedy');

% rect = getrect;
rect = [440         746        1178         489];
I = imcrop(I,round(rect));
imshow(I);
I = im2uint8(I);
if (size(I,3)==1)
    I = repmat(I,[1 1 3]);
end

num_scales = opts.num_scales;
scale_offset = opts.scale_offset;
opts.features.which_features = {'hue opp'};
%options are 'hue opp', 'saturation opp','brightness opp','hsb opp';
getRotMat; % calculate the rotation matrix 
opts.features.rotation_matrix = rotation_matrix;
f_maps = getFeatures_oppCol(double(I),num_scales+scale_offset,opts.features.which_features{1},opts);

opts.sig = 5;
% opts.features.which_features = {'luminance'};
% f_maps = getFeatures(double(I), num_scales+scale_offset,opts.features.which_features{1},opts);
Nsamples = opts.kde.Nkernels;
F = sampleF(f_maps,Nsamples,opts);
Fsym = [F; [F(:,2) F(:,1)]]; % symmetric F(A,B) = F(B,A). 
p = kde(Fsym',0.05,[],'e');
%% Interactive selection of red, green, blue circles on the zebra
figure; imshow(I); hold on;

[rowSub,colSub] = ginput;
rowSub = round(rowSub); colSub = round(colSub); 
c_vecs = {'r','r','g','g','w','w'};
% shape inserter for 6 combinations: 
% pink-pink: white circle, purple purple: green circle, white white: red
% circle, pink purple: white square, pink white: red square, purle-white:
% green square

for i =1:floor(length(colSub)/2)
    coord1 = [rowSub(2*(i-1)+1),colSub(2*(i-1)+1)];
    coord2 = [rowSub(2*i),colSub(2*i)];
    r_shape = 10; % radius of the shape drawn
    plot([coord1(1),coord2(1)],[coord1(2),coord2(2)],'o','MarkerSize',4,...
        'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
    if mod(i,2) == 1
        viscircles((coord1 + coord2)/2,10, 'EdgeColor',c_vecs{i},'DrawBackgroundCircle',false);
    else       
        rectangle('Position',[(coord1(1) + coord2(1))/2 - r_shape ...
            (coord1(2) + coord2(2))/2 - r_shape 2*r_shape 2*r_shape],...,
            'LineWidth', 2, 'EdgeColor',c_vecs{i})
    end
end

hold off;

%% create the contour plots for P_AB
% look up help kde/hist
tol = opts.kde.kdtree_tol;
dataRange = max(f_maps(:)) - min(f_maps(:)); 
x = min(f_maps(:)):dataRange/500:max(f_maps(:)); 
y = x;

[X,Y] = meshgrid(x,y);
Fim = [X(:),Y(:)];
pd = evaluate(p,Fim',tol);
pd_mesh =  reshape(pd, size(X));
%[pd,x,y] = hist(p,500,[1,2]);
% 
% figure;mesh(x,y,pd_mesh); axis square; colorbar;
% xlabel('Luminance A'); ylabel('Luminance B');
% set(gca,'XTick',0:0.1:1)
% set(gca,'YTick',0:0.1:1)
% 

% figure;contourf(x,y,pd_mesh,30); axis square; colorbar;
% xlabel('Luminance A'); ylabel('Luminance B');
% set(gca,'XTick',0:0.1:1);set(gca,'YTick',0:0.1:1)
% 
%% 
%linearInd = sub2ind(size(f_maps{1}), rowSub, colSub);
% if later not work, add {1} behind f_maps
linearInd = sub2ind(size(f_maps), colSub, rowSub);
feats = f_maps(linearInd);
colormap(jet);
%% joint probabilities
reg = opts.p_reg;
pJoint = reg + pd;
% in a mesh
pJoint_mesh = reshape(pJoint, size(X));
% figure;mesh(x,y,log(pJoint_mesh)); axis square; colorbar;
% xlabel('Luminance A'); ylabel('Luminance B');
% 

figure;contourf(x,y,log(pJoint_mesh),30); axis square; colorbar;
%xlabel('Luminance A'); ylabel('Luminance B');
xlabel(opts.features.which_features); ylabel(opts.features.which_features);
%set(gca,'XTick',x);set(gca,'YTick',y)

hold on;
for i =1:floor(length(feats)/2)
    coord = [feats(2*(i-1)+1) feats(2*i)]; % coordinate of points picked interactively
    r_shape = 0.3;
    plot(coord(1),coord(2),'o','MarkerSize',4,...
        'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
    if mod(i,2) == 1
        viscircles(coord,r_shape, 'EdgeColor',c_vecs{i},'DrawBackgroundCircle',false);
    else 
        rectangle('Position', [ coord(1) - r_shape coord(2) - r_shape ...
            2*r_shape 2*r_shape], 'LineWidth', 2, 'EdgeColor',c_vecs{i});
    end
end
hold off;
%% useful functions to look up
% help kde/evaluate
% help kde/getPoints
%% evaluate p(A)p(B)
N = floor(size(Fim,2)/2); assert((round(N)-N)==0);
p2_1 = marginal(p,1:N);
p2_2 = marginal(p,N+1:(2*N));
p2 = joinTrees(p2_1,p2_2,0.5);
pMarg_x = evaluate(p2,X(:)',tol);
pMarg_y = evaluate(p2,Y(:)',tol);
pProd = pMarg_x.*pMarg_y +reg;

%% calculate pmi
% rf = learnPMIPredictor(f_maps,p,opts);
% pmi = fastRFreg_predict(Fim,rf);
pmi = log((pJoint.^(opts.joint_exponent))./pProd);
Z_pmi = reshape(pmi,size(X));
% 
% figure;mesh(x,y,Z_pmi);axis square; colorbar;
% xlabel('Luminance A'); ylabel('Luminance B');
% 
% 
figure;[C_pmi,h_pmi]=contourf(x,y,Z_pmi,20);
axis square; colorbar; 
xlabel(opts.features.which_features); ylabel(opts.features.which_features);
%set(gca,'XTick',x);set(gca,'YTick',y)

hold on;
for i =1:floor(length(feats)/2)
    coord = [feats(2*(i-1)+1) feats(2*i)]; % coordinate of points picked interactively
    r_shape = 0.3;
    plot(coord(1),coord(2),'o ','MarkerSize',4,...
        'MarkerEdgeColor',c_vecs{i},'MarkerFaceColor',c_vecs{i});
    if mod(i,2) == 1
        viscircles(coord,r_shape, 'EdgeColor',c_vecs{i},'DrawBackgroundCircle',false);
    else 
        rectangle('Position', [ coord(1) - r_shape coord(2) - r_shape ...
            2*r_shape 2*r_shape], 'LineWidth', 2, 'EdgeColor',c_vecs{i});
    end
end
hold off; 

%%
rf = learnPMIPredictor(f_maps,p,opts);
pmi = fastRFreg_predict(Fim,rf);
