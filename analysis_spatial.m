%Check that you have the ROIs in memory
if ~exist('component_thistories','var')
    %If not, load the .mat file if it's been processed by
    %extractROIsandTraces
    try 
        load('ARCONSLICE_ALL_VARS.mat')
    catch
        %else, run the the full script 
        extractROIsAndTraces;
    end
end

%% Making the matrices
corrMatrix = zeros(size(component_thistories,1),size(component_thistories,1));
distMatrix = zeros(size(component_thistories,1),size(component_thistories,1));
disp('Now calculating time-based correlations between connected components');

max_anticorrelation = 100;
max_correlation = -100;

%setting the lowpass at .1*Fs=33Hz
[b a] = butter(4,.001,'high'); 

for i = size(corrMatrix,1):-1:1
    for j = size(corrMatrix,1):-1:1 %1:i
        %Get the distance between the components
        distMatrix(j,i) = norm(component_centroids(i,:)-component_centroids(j,:),2);
        
        sig_i = filtfilt(b,a,component_thistories(i,:));
        sig_j = filtfilt(b,a,component_thistories(j,:));
        corrValue = corr2(sig_i,sig_j);
        corrMatrix(j,i) = corrValue;
        
        if i~=j && corrValue>max_correlation
            max_correlation = corrValue;
        end
        
        if corrValue < max_anticorrelation 
            max_anticorrelation = corrValue;
        end
    end
end




%% Now let's visualize this scatter plot

%we'll remake the matrices within some parameters
threshold_distance = 40;
linethickness_max = 10.;
figure;
imagesc(mean(data,3)); %or final_mask
hold on;
for i = 1:size(component_centroids,1)
   plot(component_centroids(i,2),component_centroids(i,1),'r.');
end

title('Nearby connections of positive and negative correlation');



for i = size(corrMatrix,1):-1:1
    for j = 1:i-1
        distVal = norm(component_centroids(i,:)-component_centroids(j,:),2);
        
        if distVal>threshold_distance
            continue
        end
        
        corrVal = corr2(component_thistories(i,:),component_thistories(j,:) );
        if corrVal==1
        barf()
        end
        color = 'g';
        thickness = linethickness_max*corrVal/max_correlation;        
        if corrVal <0
            color = 'r';
            thickness = linethickness_max*corrVal/max_anticorrelation;
        end
        centr_i = component_centroids(i,:);
        centr_j = component_centroids(j,:);
        %plot the line between the two: green if positive corr, red if neg
        plot([centr_i(2),centr_j(2)],[centr_i(1),centr_j(1)],color,'LineWidth',thickness);

        
    end
end

hold off;
    

%% Querying specific pairs of puncta

% Coords for anticorr demo
% coordyx_1 = [222,219];
% coordyx_2 = [231,216];
%posscor demo
coordyx_1 = [50.3,237.5];
coordyx_2 = [34.5,248.2];
thumbnail_region = 100;

%Because the centroids are floats
%this requires searchign through for the closest centroids
min_distance = 10000;
min_idx = 0;
for i = 1:size(component_centroids,1)
    distVal = norm(coordyx_1-component_centroids(i,:),2);
    if distVal < min_distance
        min_distance = distVal;
        min_idx = i;
    end
end
centroid_idx_1 = min_idx;

min_distance = 10000;
min_idx = 0;
for i = 1:size(component_centroids,1)
    distVal = norm(coordyx_2-component_centroids(i,:),2);
    if distVal < min_distance
        min_distance = distVal;
        min_idx = i;
    end
end
centroid_idx_2 = min_idx;

figure;
subplot(2,1,1)

history_1 = component_thistories(centroid_idx_1,:);
history_1 = history_1/max(history_1);

history_2 = component_thistories(centroid_idx_2,:);
history_2 = history_2/max(history_2);

t = ([1:length(history_1)]-1)*1/Fs;
plot(t,history_1,'color','g')
hold on;
plot(t,history_2,'color','r')
hold off;
legend('Time History Puncta 1','Time History Puncta 2')
title('Comparative time history of two puncta');
xlabel('Time (seconds)')
ylabel('Sum of pixel values in puncta mask region');

subplot(2,1,2)
min_y = min(coordyx_1(1),coordyx_2(1)) - thumbnail_region;
max_y = max(coordyx_1(1),coordyx_2(1)) + thumbnail_region;
min_x = min(coordyx_1(2),coordyx_2(2)) - thumbnail_region;
max_x = max(coordyx_1(2),coordyx_2(2)) + thumbnail_region;

avg = mean(data,3);
imagesc(avg);
hold on;
plot(component_centroids(centroid_idx_1,2),component_centroids(centroid_idx_1,1),'g.');
plot(component_centroids(centroid_idx_2,2),component_centroids(centroid_idx_2,1),'r.');
hold off;
legend('Puncta 1', 'Puncta 2');
% (min_y:max_y,min_x:max_x)

%across each pixel 


%% Create some correlation-based images to show 
rois_to_compare = [centroid_idx_1,centroid_idx_2,24,32];
for roidx = rois_to_compare
target_signal = component_thistories(roidx,:)';
heatmap = zeros(size(data,1),size(data,2));
for i = 1:size(heatmap,1)
    for j = 1:size(heatmap,2)
        heatmap(i,j) = corr2(target_signal,squeeze(data(i,j,:)));
    end
end
figure; imagesc(heatmap); title(sprintf('ROI idx %i',roidx));
end

%% Outputting everything

%Save every single figure
h = get(0,'children');
for i=1:length(h)
  saveas(h(i), ['Arcon2001_output_' num2str(i)], 'fig');
end

%Make a big pile of results to save to file

% component_thistories
% avg_intensity_projection
% std_intensity_projection
% final_mask
avg_intensity_projection = mean(data,3);
std_intensity_projection = std(data,[],3);



%% Exploration in factoring correlation matrix

Z = linkage(corrMatrix,'complete','correlation');
num_clusters = 2;
c = cluster(Z,'maxclust',num_clusters);
[~,sorted_indices] = sort(c);

values = unique(c);
instances = histc(c,values)

figure;
subplot(3,1,1);
imagesc(corrMatrix); colormap jet; colorbar;
title('Unclustered correlation matrix');
subplot(3,1,2);
imagesc(corrMatrix(sorted_indices,sorted_indices));
colormap jet; colorbar;
title(sprintf('Clustered correlation matrix, C=%i',num_clusters));


color_cues = {'k.','kx','k*','kv','k*','ks','k^','kd'};
subplot(3,1,3)
imagesc(mean(data,3)); 
hold on;
for cluster_idx = 1:num_clusters
    [indices_clustermember] = find(c==cluster_idx);
    for i = 1:size(indices_clustermember)
        plot(component_centroids(indices_clustermember(i),2), ...
             component_centroids(indices_clustermember(i),1),color_cues{cluster_idx});
    end
end
hold off;
title(['Mask with ' num2str(size(component_centroids,1)) ' legitimate components']);

%% Plot the time history of the ones marked 'x'

figure;
cluster_idx = 1; %the one that is anti-correlated with nehibods
[indices_clustermember] = find(c==cluster_idx);
% [b a] = butter(4,.001,'high');
for i = 1:size(indices_clustermember)
    tsignal = component_thistories(indices_clustermember(i),:) ...
        /component_sizes(indices_clustermember(i));

%     tsignal = tsignal - min(tsignal);
%     tsignal = tsignal./max(tsignal);
%     tsignal = filtfilt(b,a,tsignal);

    
    plot(tsignal);
    hold on;
    
end
title(sprintf('Total plot for cluster_idx=%i',cluster_idx));


%% Capture the global spiking events 

global_thistory = sum(component_thistories,1);
[b a] = butter(4,.05,'low'); %setting the lowpass at .1*Fs=33Hz
global_thistory_trend = filtfilt(b,a,global_thistory);
global_pks = getSpikes(global_thistory);

t = ([1:length(global_thistory)]-1)*1/Fs;

figure;
plot(t,global_thistory);
hold on;
plot(t,global_thistory_trend,'r')
plot(t(global_pks),global_thistory(global_pks),'g*'); 
hold off;
title('Sum of all histories compared to a low-pass version of itself');


%% Go through all time histories, and create a spike chart
spike_field = zeros(size(component_thistories)+[1,0]);
spike_field(1,global_pks) = 1;
figure;

for c_idx=1:size(component_thistories,1)
    pks = getSpikes(component_thistories(c_idx,:));
    
    spike_field(c_idx+1,pks) = 1;

    
    plot(t,global_thistory);
    hold on;
    plot(t,global_thistory_trend,'r')
    plot(t(global_pks),global_thistory(global_pks),'g*'); 
    hold off;
    title('Sum of all histories compared to a low-pass version of itself');

    pause
end

%blur each of the row

windowWidth = int16(5);
halfWidth = windowWidth / 2;
gaussFilter = gausswin(5);
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

% Do the blur.
for c_idx=1:size(spike_field,1)
    spike_field(c_idx,:) = conv(spike_field(c_idx,:), gaussFilter,'same');
end

spike_similarity_matrix = zeros(size(spike_field,1),size(spike_field,1));
for i = 1:size(spike_field,1)
    for j = 1:size(spike_field,1)
        if i==j
            continue
        end
        spike_similarity_matrix(i,j) = spike_field(i,:)*spike_field(j,:)';
    end
end
%%
Z = linkage(spike_similarity_matrix,'complete');
num_clusters = 4;
c = cluster(Z,'maxclust',num_clusters);
[~,sorted_indices] = sort(c);

values = unique(c);
instances = histc(c,values)

figure;
subplot(2,1,1);
imagesc(spike_similarity_matrix); colormap jet; colorbar;
title('Unclustered spike similarity');
subplot(2,1,2);
imagesc(spike_similarity_matrix(sorted_indices,sorted_indices));
colormap jet; colorbar;
title(sprintf('Clustered spike similarity matrix, C=%i',num_clusters));



%% Pick a region of time that we know has global spikes:
time_subset_start = 25410;
time_subset_end = 36960;

figure;
plot(t(time_subset_start:time_subset_end), global_thistory(time_subset_start:time_subset_end));
title('Signal sum of all rois');
for c_idx=1:size(component_thistories,1)
   if mod(c_idx,7)==1
    figure;
   end
   subplot(7,1,mod(c_idx,7)+1)
   plot(t(time_subset_start:time_subset_end), component_thistories(c_idx,time_subset_start:time_subset_end))
   title(sprintf('%i',c_idx));
   axis off
end
%% Sorting the ROIs by most similar to the global signal
similarity_to_global = spike_similarity_matrix(1,:);
[~,I] = sort(similarity_to_global,'descend'); 

figure; plot(similarity_to_global(I))
title('Similarity metric (spiking) of vectors compared to the total sum of ROI');
xlabel('Index of ROI');


for sorted_idx=1:10
   c_idx = I(sorted_idx);
   if mod(sorted_idx,5)==1
    figure;
    subplot(6,1,1)
    plot(t(time_subset_start:time_subset_end), global_thistory(time_subset_start:time_subset_end),'r');
    axis off;
   end
   subplot(6,1,mod(sorted_idx,5)+2)
   plot(t(time_subset_start:time_subset_end), component_thistories(c_idx,time_subset_start:time_subset_end))
   pos= component_centroids(c_idx,:);
   title(sprintf('%i. size = %i, location: (%i, %i,)',c_idx, component_sizes(c_idx),pos(1),pos(2)));
   axis off
end
