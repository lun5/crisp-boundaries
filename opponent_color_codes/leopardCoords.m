coords = [195,102;120,223;168,257];
%coords = [114,45;71,35;168,87;169,44];
coords_A = [coords(:,1) - 3, coords(:,2) ];
coords_B = [coords(:,1) + 3, coords(:,2) ];
rowSub = [coords_A(:,1); coords_B(:,1)];
colSub = [coords_A(:,2); coords_B(:,2)];

dist = 7;
for i = 1:length(coords)
    rowSub(2*i-1) = coords(i,1) - dist;
    rowSub(2*i) = coords(i,1) + dist;
    colSub(2*i-1) = coords(i,2) ;
    colSub(2*i) = coords(i,2) ;
end
    
