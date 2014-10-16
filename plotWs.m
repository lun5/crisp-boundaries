coords = [195,102;120,223;168,257];
%coords = [114,45;71,35;168,87;169,44];
imshow(I);
[rowSub,colSub] = ginput;
coords = round([rowSub,colSub]);
% dist = 3;
% for i = 1:length(coords)
%     rowSub(2*i-1) = coords(i,1) - dist;
%     rowSub(2*i) = coords(i,1) + dist;
%     colSub(2*i-1) = coords(i,2) ;
%     colSub(2*i) = coords(i,2) ;
% end
    
Ws_current = Ws{1};
% black_green = [31,162];
% white_green = [44,133];
% green_green = [20,120];
% black_white = [72,152];
% coords = [black_green;white_green;green_green;black_white];

figure; imshow(I); hold on;
%coords = ginput;coords = round(coords);
plot(coords(:,1), coords(:,2),'wx','MarkerSize',10);
ind_coords = sub2ind([size(I,2),size(I,1)],coords(:,1),coords(:,2));
hold off;

Ws_spots = Ws_current(ind_coords,:)*100;
aff_im = zeros(size(I(:,:,1)));
figure
u = uicontrol('Style','slider','Position',[10 50 20 340],...
    'Min',1,'Max',16,'Value',1);
for i = 1:size(coords,1)
    im = reshape(Ws_spots(i,:),size(I,2), size(I,1));
    %figure;
    h = histogram(nonzeros(im),'DisplayStyle','stairs' );
    %h.FaceColor = [0.8 .8 .8]; h.BinWidth= max(nonzeros(Ws_spots))/50;
    h.Normalization = 'probability'; 
    ax = axis;axis([min(nonzeros(Ws_spots)) max(nonzeros(Ws_spots)) 0 1]);
    u.Value = i;
    M(i) = getframe(gcf);    
    aff_im = aff_im + im';
end
% movie(M,5,3);
% movie2avi(M,fullfile(pwd,'results','opp.avi'), 'compression', 'None','fps',3);
% movie2avi(M,fullfile(pwd,'results','pink_pink.avi'), 'compression', 'None','fps',3);

figure; imshow(I); hold on; 
cspy(aff_im,'markersize',15, 'colormap', 'jet', 'levels', 7); 
colorbar; 
hold off;

figure;
h = histogram(aff_im(aff_im > 0),'DisplayStyle','bar' );
h.FaceColor = [0.8 .8 .8];h.NumBins = 20;
h.Normalization = 'probability'; ax = axis;
axis([0 max(nonzeros(Ws_spots)) 0 1]);
