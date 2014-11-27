[colSub,rowSub] = ginput;
rowSub = round(rowSub); colSub = round(colSub); 
linearInd = sub2ind(size(f_maps), rowSub, colSub);