% Demo of crisp boundaries
%
% -------------------------------------------------------------------------
% Crisp Boundaries Toolbox
% Phillip Isola, 2014 [phillpi@mit.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------


%% setup
% first cd to the directory containing this file, then run:
% compile; % this will check to make sure everything is compiled properly; if it is not, will try to compile it


%% Detect boundaries
% you can control the speed/accuracy tradeoff by setting 'type' to one of the values below
% for more control, feel free to play with the parameters in setEnvironment.m

%type = 'super_speedy'; % use this for fastest results
type = 'speedy'; % use this for slightly slower but crisper results
%type = 'accurate'; % use this for slow, but high quality results
%I = imread('test_images/101027.jpg');
%I = imread('test_images/253027.jpg'); % zebra
I = imread('test_images/134067.jpg'); % leopard
datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
if ~ exist(datadir,'dir')
    datadir = '/Users/lun5/Research/color_deconvolution/TissueImages/';
end
I = imread(fullfile(datadir,'tp10-867-1_47104_22528_2048_2048.tif'));
rect = [911        1324         207         129];
I = imcrop(I,round(rect));
% datadir = 'C:\Users\luong_nguyen\Dropbox\DDIresearch\crisp-boundaries\breastColony';
% I = imread(fullfile(datadir, 'LLV232 D04 10x max proj.tif'));
% I = imadjust(I);

[E,E_oriented] = findBoundaries(I,type);
close all; subplot(121); imshow(I); subplot(122); imshow(1-mat2gray(E));

%% Segment image
% builds an Ultrametric Contour Map from the detected boundaries (E_oriented)
% then segments image based on this map
%
% this part of the code is only supported on Mac and Linux

if (~ispc)
    
    thresh = 0.1; % larger values give fewer segments
    E_ucm = contours2ucm_crisp_boundaries(E_oriented,type);
    S = ucm2colorsegs(E_ucm,I,thresh);

    close all; subplot(121); imshow(I); subplot(122); imshow(S);
end