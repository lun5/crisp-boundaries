%# some sample images
close all; clear all;
I = imread('coins.png');
I_transp = imread('peppers.png');

%# create a gaussian mask for transparency
[r,c,~] = size(I_transp);
M = fspecial('gaussian', [r c], mean([r c]./5));
M = (M-min(M(:)))./range(M(:));

%# show overlayed images
figure, imshow(I, 'XData',[1 c], 'YData',[1 r]), hold on
hImg = imshow(I_transp);
set(hImg, 'AlphaData',M);

%# draw a rectangle
rectangle('Position',[255 120 100 100], 'LineWidth',2, 'EdgeColor','b');

rectangle('Position',[0.59,0.35,3.75,1.37],...
          'Curvature',[0.8,0.4],...
         'LineWidth',2,'LineStyle','--')
daspect([1,1,1])

