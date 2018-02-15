function [ X_mask ] = compressDataWithMask(final_mask,data)
    %COMPRESSDATAWITHMASK Summary of this function goes here
    %   Detailed explanation goes here
    [mask_y,mask_x] = find(final_mask);
    X_mask = zeros(length(mask_y),size(data,3));
    ctr = 1;

    for mask_idx = 1:length(mask_y)
        X_mask(mask_idx,:) = data(mask_y(mask_idx),mask_x(mask_idx),:);
    end

end

