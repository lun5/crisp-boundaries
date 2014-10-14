Ws_current = Ws{1};
% black_green = [31,162];
% white_green = [44,133];
% green_green = [20,120];
% black_white = [72,152];
% coords = [black_green;white_green;green_green;black_white];

figure; imshow(I); hold on;
%coords = ginput;coords = round(coords);
plot(coords(:,1), coords(:,2),'kx','MarkerSize',10);
ind_coords = sub2ind([size(I,2),size(I,1)],coords(:,1),coords(:,2));
hold off;

aff_im = zeros(size(I(:,:,1)));
for i = 1:length(coords)
    im = reshape(Ws_current(ind_coords(i),:),size(I,2), size(I,1));
    aff_im = aff_im + im';
end

aff_im_discr = aff_im;
aff_im_discr(aff_im > 0 & aff_im < 450) = 1;
aff_im_discr(aff_im > 450) = 3;
figure; imshow(I); hold on; 
cspy(aff_im,'markersize',5, 'colormap', 'jet', 'levels', 7); 
cspy(aff_im_discr,'markersize',5, 'colormap', 'jet', 'levels', 7); 
colorbar; 
hold off;

figure;
h = histogram(aff_im(aff_im > 0),'DisplayStyle','bar' );
h.FaceColor = [0.8 .8 .8];

