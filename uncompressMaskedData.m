function [ output] = uncompressMaskedData(mask,X_masked,idx_start,idx_end)

if ~exist(idx_start,'var')
    idx_start = 1;
end
if ~exist(idx_end,'var')
    idx_end = size(X_masked,2);
end


output = zeros(size(mask,1),size(mask,2),idx_end-idx_start+1);

[mask_y,mask_x] = find(mask);

for mask_idx = 1:length(mask_y)
    output(mask_y(mask_idx),mask_x(mask_idx),:) = X_masked(mask_idx,:);
end

end

