Ws_current = Ws{1};

imshow(I);hold on;
%coords = ginput; coords = round(coords);
coords = [193,102;122,226;168,257];
plot(coords(:,1), coords(:,2),'ko', 'MarkerSize', 5);
c_vecs = {'b','g','r'};
r_shape = 5;
for i = 1:length(linInd)
    rectangle('Position',[(coords(i,1) - r_shape) ...
    (coords(i,2)- r_shape) 2*r_shape 2*r_shape],...,
    'LineWidth', 2, 'EdgeColor',c_vecs{i})
end

linInd = sub2ind([size(I,2) size(I,1)], coords(:,1), coords(:,2));

aff_im = zeros(size(I(:,:,1)));
for i = 1:length(linInd)
    im = reshape(Ws_current(linInd(i),:),size(I,2),size(I,1));
    aff_im = aff_im + im';
end

figure;
imshow(I); hold on;
cspy(aff_im*100,'markersize',15,'colormap','jet', 'levels',7);
hold off;
colorbar

figure;
h = hist(aff_im(aff_im > 0));
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';