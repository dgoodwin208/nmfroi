input_file = '../Arcon/2001_substack.tif';
output_directory = 'output/';
if ~exist(output_directory,'dir')
    mkdir(output_directory);
end
Fs = 300; %Hz

data = load3DTif_uint16(input_file);
output_filename = 'output/sparsenmfnnls_pieces.mat';

%For robustness, we separate the data into overlapping sections and
%calculate the signal mask (1=pixel is signal, 0 is signal is background),
%then combine the different masks. The number of overlapping sections is
%num_pieces-1
num_pieces = 5;


%% restructure the data vector to put each pixel in a row

elements_per_subpiece = size(data,3)/num_pieces;
frequency_resolution = Fs/elements_per_subpiece;
sprintf('Frequency domain resolution %f',frequency_resolution);
BigA = [];
BigY = [];
for n=1:num_pieces-1
    %note that this goes for two spans, effectively 100% overlap
    t_indices = (n-1)*elements_per_subpiece+1:1:(n+1)*elements_per_subpiece;
    t = squeeze(data(1,1,t_indices));
    
    window = hanning(length(t));
    
    X_mags = abs(fft(t));
    num_bins = length(X_mags);
    log(X_mags(1:num_bins/2));
    
    X = zeros(size(data,1)*size(data,2),num_bins/2);
    ctr = 1;
    for y=1:size(data,1)
        for x=1:size(data,2)
            t = squeeze(data(y,x,t_indices) - mean(data(y,x,t_indices)));
            t = t.*window;
            T = abs(fft(t));
            X(ctr,:) = T(1:num_bins/2);
            ctr = ctr +1;
        end
    end
    
    
    k = 2;
    [A,Y,numIter,tElapsed,finalResidual]=sparsenmfnnls(X,k); %removing sparsity
    
    
    %Save all the results
    BigA(n,:,:) = A;
    BigY(n,:,:) = Y;
    
    
    % unpack this (yes there should be a way to do this with reshape!)
    masks = zeros(size(data,1),size(data,1),2);
    ctr = 1;
    for y=1:size(data,1)
        for x=1:size(data,2)
            masks(y,x,1) = A(ctr,1);
            masks(y,x,2) = A(ctr,2);
            ctr = ctr +1;
        end
    end
    
    figure;
    subplot(2,2,1)
    imagesc(masks(:,:,1));
    axis off;
    title(['Column 1 of A, subpiece #' num2str(n)]);
    
    subplot(2,2,2)
    imagesc(masks(:,:,2));
    axis off;
    title('Column 2 of A');
    
    subplot(2,2,3)
    plot(Y');
    legend('First component', 'Second component')
    xlabel('Frequency index')
    ylabel('Magnitude of Frequency domain')
    
    figure;
    imagesc(masks(:,:,1)<masks(:,:,2))
    title(['Arbitrary choice of mask, subpiece #' num2str(n)]);
end

save(output_filename,'BigA','BigY','-v7.3');

%% Somehow take user input for a point that we know is part of the target ROI
keypoint_x = 174;
keypoint_y = 74;


signal_masks = [];
for n=1:num_pieces-1
    
    %Convert the A vectors into masks
    masks = zeros(size(data,1),size(data,1),2);
    ctr = 1;
    for y=1:size(data,1)
        for x=1:size(data,2)
            masks(y,x,1) = BigA(n,ctr,1);
            masks(y,x,2) = BigA(n,ctr,2);
            ctr = ctr +1;
        end
    end
    
    %Then at the manually defined keypoint, get the index of selecting mask
    
    [t,signal_component_idx] = max( [masks(keypoint_y,keypoint_x,1),masks(keypoint_y,keypoint_x,2)]);
    
    
    
    non_signal_component_idx = 1;
    if non_signal_component_idx == signal_component_idx
        non_signal_component_idx =2;
    end

    
    
    %Create the binary mask
    signal_masks(n,:,:) = masks(:,:,signal_component_idx)> masks(:,:,non_signal_component_idx );
end

figure;
for n=1:num_pieces-1
    subplot(num_pieces-1,1,n)
    imagesc(squeeze(signal_masks(n,:,:)))
    axis off;
    title(sprintf('Mask applied to subset of data %i',n));
end

figure;
summed = squeeze(sum(signal_masks,1));
imagesc(summed)
title(sprintf('Summation of the binary mask from the %i time subsections',num_pieces-1));

%The mask we use is all non-zero pixel values, which means that at any
%subsection of time, the pixel was classified as signal.
final_mask = summed>0;

%% Make a mask of this data

output = data;
noise = output;
for z=1:size(output,3)
    frame = data(:,:,z);
    frame(~final_mask)=0; %keep only the signal
    noise(final_mask)=0;  %keep only the noise
    output(:,:,z) = frame;
end

%If you want, you can save the mask file or the masked data here.
% save3DTif_uint16(output,'Slice4-013-mask');

%% Use the mask data to create a new tif stack

%Compress the 3D structure into a 2D structure of linear indices for (x,y)
%across Z: this is if we want to further segment with
%nmf/kmeans/correlation etc.

[mask_y,mask_x] = find(final_mask);
X_mask = zeros(length(mask_y),size(data,3));
ctr = 1;

for mask_idx = 1:length(mask_y)
    X_mask(mask_idx,:) = data(mask_y(mask_idx),mask_x(mask_idx),:);
end


time_signal = sum(X_mask,1);
%This is a much smaller matrix, print just how much smaller:
percentage_masked = sum(sum(final_mask))/(size(data,1)*size(data,2))*100;
disp(['Created masked matrix which is ' num2str(percentage_masked) '% the size of the original'])

figure;
subplot(4,2,1)
imagesc(mean(data,3));
title('Avg in t-domain of original image')

subplot(4,2,2)
imagesc(final_mask);
title(['Mask of data, ' num2str(percentage_masked) '% original data']);
subplot(4,2,3);
t = ([1:length(time_signal)]-1)*1/Fs;
plot(t,time_signal);
title('T-profile of masked components');
xlabel('time (seconds)');
ylabel('sum of pixel values per frame');

subplot(4,2,4);
t = time_signal - mean(time_signal);
window = hanning(length(t));
t = t.*window';
T = abs(fft(t));
num_bins = length(T);
freq_signal = T(1:num_bins/2);
plot(log(freq_signal))
title('Frequency domain')
xlabel('frequency index: nyquist at 2.2Hz');
ylabel('Log magnitude')

subplot(4,2,5);
size_data = size(data);
time_signal = sum(reshape(data, prod(size_data([1 2])), []));
t = ([1:length(time_signal)]-1)*1/Fs;
plot(t,time_signal);
title('Sum of all pixels in the frame');
xlabel('time (seconds)');
ylabel('sum of pixel values per frame');

subplot(4,2,6);

t = time_signal - mean(time_signal);
window = hanning(length(t));
t = t.*window';
T = abs(fft(t));
freq_signal = T(1:num_bins/2);
plot(log(freq_signal));
title('Frequency domain')
xlabel('frequency index: nyquist at 2.2Hz');
ylabel('Log magnitude')

subplot(4,2,7);
size_data = size(noise);
time_signal = sum(reshape(noise, prod(size_data([1 2])), []));
t = ([1:length(time_signal)]-1)*1/Fs;
plot(t,time_signal);
title('Sum of all excluded pixels in mask');
xlabel('time (seconds)');
ylabel('sum of noise pixel values per frame');

subplot(4,2,8);

t = time_signal - mean(time_signal);
window = hanning(length(t));
t = t.*window';
T = abs(fft(t));
freq_signal = T(1:num_bins/2);
plot(log(freq_signal));
title('Frequency domain')
xlabel('frequency index: nyquist at 2.2Hz');
ylabel('Log magnitude')



%% Use the connected components code in MATLAB to aggregate pixels into discrete units

CC = bwconncomp(final_mask);

component_thistories = []; %zeros(length(CC.PixelIdxList), size(data,3));
component_centroids = [];
component_sizes = [];
ctr_thist = 1;
disp(['Found ' num2str(length(CC.PixelIdxList)) ' components']);
sizes = [];

CC_IGNORE_MIN = 6;
CC_IGNORE_MAX = 2000;

centroids_mask = zeros(size(data,1),size(data,2));
for i = 1:length(CC.PixelIdxList)
    oneD_indices = CC.PixelIdxList{i};
    [y_indices,x_indices] = ind2sub([size(data,1),size(data,2)],oneD_indices);
    
    num_pixels_in_component = length(y_indices);
    
    
    
    sizes(i) = num_pixels_in_component;
    
    disp(['CC index #' num2str(i) ' has ' num2str(num_pixels_in_component) ' pixels']);
    
    if num_pixels_in_component <= CC_IGNORE_MIN || ...
            num_pixels_in_component >= CC_IGNORE_MAX
        disp('Skipping due to size');
        continue
    end
    
    data_subset = zeros(size(y_indices,1),size(data,3));
    for ctr=1:length(y_indices)
        data_subset(ctr,:) = data(y_indices(ctr),x_indices(ctr),:);
        %create the final mask via centroids
        centroids_mask(y_indices(ctr),x_indices(ctr),:)=1;
    end
    
    
    [comp_space,comp_time] = nnmf(data_subset,1);
    
    
    
    %check the sum to make sure it's not an empty vector
    %candidate_vector = squeeze(sum(data_subset,1));
    candidate_vector = comp_time; %Do rank1-approximation to remove some noise
    if sum(diff(candidate_vector))~=0
        component_thistories(ctr_thist,:)= candidate_vector;
        component_centroids(ctr_thist,:) = [mean(y_indices),mean(x_indices)];
        component_sizes(ctr_thist) = num_pixels_in_component;
        ctr_thist = ctr_thist +1;
    else
        disp('Skipping due to lack of activity in the area');
    end
    %     disp(num2str(i));
end

disp(['Found ' num2str(size(component_centroids,1)) ' legitimate components']);

%If you want to explore statistics of the size of the ROIs
% figure;
% cc_sizes = [ [1,5,10], 20:10:400];
% histogram(sizes,cc_sizes);
% title('Distribution of sizes of connected components, in pixels');

figure;
subplot(3,1,1)
imagesc(final_mask);
hold on;
for i = 1:size(component_centroids,1)
    plot(component_centroids(i,2),component_centroids(i,1),'r.');
end
hold off;
title(['Mask with ' num2str(size(component_centroids,1)) ' legitimate components']);

subplot(3,1,2)

imagesc(mean(data,3));colormap gray
hold on;
for i = 1:size(component_centroids,1)
    plot(component_centroids(i,2),component_centroids(i,1),'r.');
end
hold off;

subplot(3,1,3)
imagesc(centroids_mask);
title('Final refined mask after picking regions');


%%
clear data %remove the biggest item
save('ARCONSLICE_ALL_VARS.mat')
