coords = [193,102;122,226;168,257];
coords_A = [coords(:,1) - 3, coords(:,2) - 3];
coords_B = [coords(:,1) + 3, coords(:,2) + 3];
rowSub = [coords_A(:,1); coords_B(:,1)];
colSub = [coords_A(:,2); coords_B(:,2)];

for i = 1:length(coords)
    rowSub(2*i-1) = coords(i,1) - 3;
    rowSub(2*i) = coords(i,1) + 3;
    colSub(2*i-1) = coords(i,2) - 3;
    colSub(2*i) = coords(i,2) + 3;
end
    
