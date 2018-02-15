%Check that you have the ROIs in memory
if ~exist('final_mask','var')
    %If not, load the .mat file if it's been processed by
    %extractROIsandTraces
    try
        load('ARCONSLICE_ALL_VARS.mat')
    catch
        %else, run the the full script
        extractROIsAndTraces;
    end
end

%We also need the raw data ready:
%Variabels set in the extractROIsAndTraces.m file
data = load3DTif_uint16(input_file);
filtered_data = bsxfun(@times,data,~final_mask);
%% Use the inverse of the mask to get the background information

background_trace = squeeze(mean(mean(filtered_data,2),1));
background_space = reshape(data,[],size(data,3))*background_trace/norm(background_trace)^2;
background_space = max(reshape(background_space,size(data,1),[]),0);

data_matrix = reshape(data,[],size(data,3)) - background_space(:)*background_trace';

%Then use the mask to get all filtered data (ie, signal pixels) in one
%matrix
data_total_masked = compressDataWithMask(final_mask,data);

%%

signal = sum(data_total_masked,1);
mf = medfilt1(signal,100,[],2);

instadff = signal./mf-1;
instadff(1) = instadff(2);
instadff(end) = instadff(end-1);

figure
subplot(3,1,1)
plot(mf);
title('Median filtered time signal, N=100');
subplot(3,1,2)
plot(signal);
title('Raw sum of masked pixel values')
subplot(3,1,3)
plot(instadff);
title('DF/F calculation via median filtering')
%% filtered signals for all centroids

filtered_centroids = zeros(size(component_thistories));

background_summed = background_trace;

for idx = 1:size(component_thistories,1)
    
    signal = component_thistories(idx,:)/component_sizes(idx) - background_summed';
    
    mf = medfilt1(signal,100,[],2);
    
    instadff = signal./mf-1;
    instadff(1) = instadff(2);
    instadff(end) = instadff(end-1);
    
    filtered_centroids(idx,:) = instadff;
end

signal = sum(component_thistories,1)./sum(component_sizes) - background_summed';
mf = medfilt1(signal,100,[],2);

instadff = signal./mf-1;
instadff(1) = instadff(2);
instadff(end) = instadff(end-1);
instadff_global=instadff; %Create the global singla

%Chop off the extremes in case of artifacts from the median filtering
% T = size(filtered_centroids,2);
% filtered_centroids = filtered_centroids(:,.05*T:.95*T);
% instadff_global = instadff_global(.05*T:.95*T);
% T = size(filtered_centroids,2);

%Edge handling:
instadff_global(1) = instadff_global(2); instadff_global(end) = instadff_global(end-1);
filtered_centroids(1,:) = filtered_centroids(2,:); filtered_centroids(end,:) = filtered_centroids(end-1,:);
figure;
imagesc(filtered_centroids); axis off;
% subplot(2,1,2)
% plot(filtered_centroids(idx+1,:));
% xlim([1,length(filtered_centroids(idx+1,:))]);

%%
%For all spikes, make an average of the time history per roi
global_pks = getSpikes(instadff_global);

%take the instadff_global from the previous cell
%And the filtered_centroids from eftychios_script.m
SAMPLES_BEFORE = 10;
SAMPLES_AFTER = 40;
SPIKE_THRESHOLD = .01;

[~,filtered_pk_indices] = find(instadff_global(global_pks)>SPIKE_THRESHOLD);

roi_spike_history = zeros(size(component_thistories,1),length(filtered_pk_indices),...
    SAMPLES_AFTER+SAMPLES_BEFORE);
global_spike_history= zeros(length(filtered_pk_indices),SAMPLES_AFTER+SAMPLES_BEFORE);
ctr = 1;
spike_t_indices = [];
for peak_idx = global_pks
    
    if instadff_global(peak_idx)>SPIKE_THRESHOLD
        start_idx = peak_idx-SAMPLES_BEFORE+1;
        end_idx = peak_idx+SAMPLES_AFTER;
        spike_t_indices(ctr,:) = [start_idx,end_idx];
        global_spike_history(ctr,:) = instadff_global(start_idx:end_idx);
        for roi_idx = 1:size(filtered_centroids,1)
            roi_spike_history(roi_idx,ctr,:) = filtered_centroids(roi_idx,start_idx:end_idx);
        end
        ctr = ctr+1;
    end
end


%%
figure
filename = 'summed_spike_history.gif';
t = ([1:SAMPLES_BEFORE+SAMPLES_AFTER]-1)*1/Fs;
n = size(global_spike_history,1);
for idx = 1:size(filtered_centroids,1)
    plot(t,mean(global_spike_history,1))
%     errorbar(t,mean(global_spike_history,1),std(global_spike_history,1)/sqrt(n));
    hold on;
    plot(t,median(global_spike_history,1),'b--');
    
    plot(t,mean(squeeze(roi_spike_history(idx,:,:)),1),'r');
%     errorbar(t,mean(squeeze(roi_spike_history(idx,:,:)),1),std(squeeze(roi_spike_history(idx,:,:)),1)/sqrt(n),'r')    
    plot(t,median(squeeze(roi_spike_history(idx,:,:)),1),'r--');
    
    
    

    %     plot(t,mean(squeeze(roi_spike_history(idx,:,:)),1))
    hold off;
    title(sprintf('ROI idx = %i',idx));
    xlabel('Time, seconds')
    ylabel('Summed time history');
    legend('Mean dendrite activity','Median dendrite activity', 'mean ROI activity','median ROI activity');
    NumTicks = 10; L =get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))

    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if idx == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
end
%% Make an image of this average

%First interpolate the values a bit:
figure;
sample_signal = sum(squeeze(roi_spike_history(1,:,:)),1);
x = 1:length(sample_signal);
xx = 1:.2:length(sample_signal);
yy = spline(x,sample_signal,xx);
plot(x,sample_signal,'-o',xx,yy);

spatial_map = zeros(size(data,1),size(data,2),length(xx));

%Now a bit of awkwardness to make sure we're getting all the right pixels
%from the connected component analysis earlier
index_to_mycomponents = 1;
image_mask = zeros(size(data,1),size(data,2));

for i = 1:length(CC.PixelIdxList)
    oneD_indices = CC.PixelIdxList{i};
    [y_indices,x_indices] = ind2sub([size(data,1),size(data,2)],oneD_indices);
    
    num_pixels_in_component = length(y_indices);
    
    if num_pixels_in_component <= CC_IGNORE_MIN || ...
       num_pixels_in_component >= CC_IGNORE_MAX 
        disp('Skipping due to size');
         continue 
    end
    
    roi_thistory = sum(squeeze(roi_spike_history(index_to_mycomponents,:,:)),1);
    %sanity check to make sure the size is what we expect
    if (length(y_indices) ~= component_sizes(index_to_mycomponents))
        barf();
    end
    
    roi_thistory_interp = spline(x,roi_thistory,xx);
    for pix_idx = 1:length(y_indices)
        spatial_map(y_indices(pix_idx),x_indices(pix_idx),:) = roi_thistory_interp;
        image_mask(y_indices(pix_idx),x_indices(pix_idx)) = 1;
    end
    
    index_to_mycomponents = index_to_mycomponents +1;
end


filename = 'video_of_behavior.gif';
pix_max = max(max(max(spatial_map))); pix_min = min(min(min(spatial_map)));

figure(2)


grayscale_content = mean(data,3) - min(min(mean(data,3)));
grayscale_content = grayscale_content/max(max(grayscale_content));

spatial_map = spatial_map - pix_min;
spatial_map = spatial_map./pix_max;
for idx = 1:size(spatial_map,3)
    
    
    C = imfuse(grayscale_content,spatial_map(:,:,idx),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    
    %imagesc(spatial_map(:,:,idx),[pix_min, pix_max]);
    imshow(C,'InitialMagnification',300);

    title(sprintf('Frame %i/%i, spike at 50',idx,size(spatial_map,3)));
    colormap gray
    axis off;
    drawnow
    frame = getframe(2);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if idx == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
end

%% Make an image to show the ROI indices in space
figure;
imagesc(mean(data,3))
hold on;
for c_idx=1:size(component_thistories,1)
    text(component_centroids(c_idx,2), component_centroids(c_idx,1), ...
        cellstr(num2str(c_idx)), 'FontSize', 15, 'Color', 'k');
    axis off
end
hold off;
%% Step through all time histories for a specific 
roi_of_interest = 4;

figure(3)
for ctr=1:size(roi_spike_history,2)
    
    plot(squeeze(roi_spike_history(roi_of_interest,ctr,:)));
    pause
end
%% Make a correlation matrix of the averaged time history
corrMatrix = zeros(size(component_thistories,1),size(component_thistories,1));
distMatrix = zeros(size(component_thistories,1),size(component_thistories,1));

max_anticorrelation = 100;
max_correlation = -100;


for i = size(corrMatrix,1):-1:1
    for j = size(corrMatrix,1):-1:1 %1:i
        if i==j
            continue
        end
        %Get the distance between the components
        distMatrix(j,i) = norm(component_centroids(i,:)-component_centroids(j,:),2);
        
        sig_i = sum(squeeze(roi_spike_history(i,:,:)),1);
        sig_j = sum(squeeze(roi_spike_history(j,:,:)),1);
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

figure;
subplot(2,1,1);
% imagesc(corrMatrix); colormap jet; colorbar;
% title('Unclustered correlation matrix');
rgb = zeros(size(data,1),size(data,2),3);
rgb(:,:,1) = mean(data,3)/max(max(max(mean(data,3)))); 
rgb(:,:,2) = rgb(:,:,1); rgb(:,:,3) = rgb(:,:,2);
for i = 1:length(CC.PixelIdxList)
    oneD_indices = CC.PixelIdxList{i};
    [y_indices,x_indices] = ind2sub([size(data,1),size(data,2)],oneD_indices);
    num_pixels_in_component = length(y_indices);
    if num_pixels_in_component <= CC_IGNORE_MIN || ...
       num_pixels_in_component >= CC_IGNORE_MAX 
        continue 
    end

    for ctr=1:length(y_indices)
        rgb(y_indices(ctr),x_indices(ctr),1)=1.;
    end
end
imshow(rgb);
title('Automatically determined ROIs');


% imagesc(mean(data,3)); colormap gray;
% hold on;
% for i = 1:size(component_centroids,1)
%    plot(component_centroids(i,2),component_centroids(i,1),'b.','MarkerSize',10);
% end
% colormap gray
% hold off;
% axis off


Z = linkage(corrMatrix,'complete','correlation');
num_clusters = 4;
c = cluster(Z,'maxclust',num_clusters);
[~,sorted_indices] = sort(c);

values = unique(c);
instances = histc(c,values)

%we'll remake the matrices within some parameters
threshold_distance = 61;
linethickness_max = 2.;
subplot(2,1,2)
rgb = zeros(size(data,1),size(data,2),3);
rgb(:,:,1) = mean(data,3)/max(max(max(mean(data,3)))); 
rgb(:,:,2) = rgb(:,:,1); rgb(:,:,3) = rgb(:,:,2);
imshow(rgb);
colormap gray;
hold on;

title('Correlation of ROI time histories at dendritic spike events');

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
        
        centr_i = component_centroids(i,:);
        centr_j = component_centroids(j,:);
        
%         color = 'g';
        thickness = linethickness_max*corrVal/max_correlation;        
        
        if corrVal <0
            color = [1 .5 0];
            thickness = linethickness_max*corrVal/max_anticorrelation;
            %plot the line between the two: green if positive corr, red if neg
            plot([centr_i(2),centr_j(2)],[centr_i(1),centr_j(1)],'color',[1 .5 0],'LineWidth',thickness);
        else 
            %plot the line between the two: green if positive corr, red if neg
            plot([centr_i(2),centr_j(2)],[centr_i(1),centr_j(1)],'g','LineWidth',thickness);
        end
        
        

        
    end
end



% for i = 1:size(component_centroids,1)
%    plot(component_centroids(i,2),component_centroids(i,1),'r.');
% end


% color_cues = {'k.','kx','k*','kv','k*','ks','k^','kd'};
% for cluster_idx = 1:num_clusters
%     [indices_clustermember] = find(c==cluster_idx);
%     for i = 1:size(indices_clustermember)
%         plot(component_centroids(indices_clustermember(i),2), ...
%              component_centroids(indices_clustermember(i),1),color_cues{cluster_idx});
%     end
% end
hold off;

figure;


imagesc(corrMatrix(sorted_indices,sorted_indices));
colormap jet; colorbar;
title(sprintf('Clustered correlation matrix, C=%i',num_clusters));

imagesc(corrMatrix);
colormap jet; colorbar;
title('Unclustered correlation matrix');



%% Create a sorting based on Victor-Purpura spike distance (cost param unclear)

VPCOST = 1;
global_signal = filtered_centroids(idx+1,:);
global_pks = getSpikes(global_signal);
t = ([1:T]-1)*1/Fs;
% figure;
% plot(t,signal)
% hold on;
% plot(t(global_pks),signal(global_pks),'g*');
% hold off;

similarity_to_global = [];
for c_idx=1:size(component_thistories,1)
    pks = getSpikes(filtered_centroids(c_idx,:));
    
    %the last parameter is the distance
    %For q=0 the distance is equal to the difference in spike counts,
    %while for large q the distance approaches the number of
    %non-coincident spikes, as it becomes more favorable to delete and
    %reinsert all non-coincident spikes rather than shifting them. Thus, by increasing the cost, the distance is transformed from a rate distance to a temporal distance.
    similarity_to_global(c_idx) = spkd(global_pks/Fs,pks/Fs,VPCOST);
end

[~,I] = sort(similarity_to_global,'ascend');

figure;
plot(similarity_to_global(I))
title(['Similarity metric (spiking) of vectors compared to the total sum of ROI cost=' num2str(VPCOST)]);
xlabel('Index of ROI');

figure;
scatter(component_sizes,similarity_to_global)
xlabel('Size of ROI')
ylabel('Correlation to global signal');

%Creating a set of total time plots
% figure;
% for peak_idx = 1:length(global_pks)
%
%    t_idx = global_pks(peak_idx);
%    indices = [max(t_idx-Fs,1):min(t_idx+Fs,size(filtered_centroids,2))];
%
%    plot(t(indices),global_signal(indices),'LineWidth',4);
%    hold on;
%    plot(t(indices),filtered_centroids(I(1:3),indices)');
%    plot(t(indices),global_signal(indices),'LineWidth',4);
%    title('Zoomed time area around a dendritic spike')
%    ylabel('dF/F');
%    xlabel('Time (seconds)')
%    hold off;
%    pause
% end

%Creating a zoom-in to explore certain peaks
for sorted_idx=1:10
    c_idx = I(sorted_idx);
    if mod(sorted_idx,5)==1
        figure;
        subplot(6,1,1)
        plot(t, global_signal,'r'); hold on
        plot(t(global_pks),global_signal(global_pks),'g*'); hold off;
        axis off;
    end
    subplot(6,1,mod(sorted_idx,5)+2)
    
    signal = filtered_centroids(c_idx,:);
    plot(t,signal); hold on;
    pks = getSpikes(signal);
    plot(t(pks),signal(pks),'g*'); hold off;
    
    pos= component_centroids(c_idx,:);
    title(sprintf('distance rank=%i, centroid=%i. size = %i, location: (%i, %i,)',sorted_idx,c_idx, component_sizes(c_idx),pos(1),pos(2)));
    axis off
end

%% Create a matrix of similarities

spk_matrix = zeros(size(component_thistories,1));

for i = 1:size(spk_matrix,1)
    row_pks = getSpikes(filtered_centroids(i,:));
    for j = i+1:size(spk_matrix,1)
        col_pks = getSpikes(filtered_centroids(j,:));
        spk_matrix(i,j) = spkd(row_pks/Fs,col_pks/Fs,VPCOST);
        spk_matrix(j,i) = spk_matrix(i,j);
    end
end

figure; imagesc(spk_matrix)
%% Sort the spike matrix into clusters and view
spk_similarity_matrix = 1./spk_matrix;
for d= 1:size(spk_matrix)
    spk_similarity_matrix(d,d) = min(min(spk_matrix));
end
Z = linkage(spk_similarity_matrix,'complete');
num_clusters = 4;
c = cluster(Z,'maxclust',num_clusters);
[~,sorted_indices] = sort(c);

values = unique(c);
instances = histc(c,values);

figure;
subplot(2,1,1);
imagesc(spk_similarity_matrix); colormap jet; colorbar;
title('Unclustered spike similarity');
subplot(2,1,2);
imagesc(spk_similarity_matrix(sorted_indices,sorted_indices));
colormap jet; colorbar;
title(sprintf('Clustered spike similarity matrix, C=%i',num_clusters));


color_cues = {'k.','kx','k*','kv','k*','ks','k^','kd'};
figure;
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

%% Take time averages of each ROI around the spike



%% Create a sorting based on correlation

global_signal = filtered_centroids(idx+1,:);
% figure;
% plot(t,signal)
% hold on;
% plot(t(global_pks),signal(global_pks),'g*');
% hold off;

similarity_to_global = [];
for c_idx=1:size(component_thistories,1)
    similarity_to_global(c_idx) = corr2(global_signal,filtered_centroids(c_idx,:));
end
[~,I] = sort(similarity_to_global,'descend');

figure; plot(similarity_to_global(I))
title('Correlation of ROIs to global signal');
xlabel('Index of ROI');

figure;
scatter(component_sizes,similarity_to_global)
xlabel('Size of ROI')
ylabel('Correlation to global signal');

for sorted_idx=1:10
    c_idx = I(sorted_idx);
    if mod(sorted_idx,5)==1
        figure;
        subplot(6,1,1)
        plot(t, global_signal,'r'); hold on
        plot(t(global_pks),global_signal(global_pks),'g*'); hold off;
        axis off;
    end
    subplot(6,1,mod(sorted_idx,5)+2)
    
    signal = filtered_centroids(c_idx,:);
    plot(t,signal); hold on;
    pks = getSpikes(signal);
    plot(t(pks),signal(pks),'g*'); hold off;
    
    pos= component_centroids(c_idx,:);
    title(sprintf('distance rank=%i, centroid=%i. size = %i, location: (%i, %i,)',sorted_idx,c_idx, component_sizes(c_idx),pos(1),pos(2)));
    axis off
end

%% plotting a spatial link using the

figure;
imagesc(mean(data,3))
hold on;
for sorted_idx=1:10
    c_idx = I(sorted_idx);
    
    text(component_centroids(c_idx,2), component_centroids(c_idx,1), ...
        cellstr(num2str(sorted_idx)), 'FontSize', 15, 'Color', 'k');
    title(sprintf('distance rank=%i, centroid=%i. size = %i, location: (%i, %i,)',sorted_idx,c_idx, component_sizes(c_idx),pos(1),pos(2)));
    axis off
end

for sorted_idx=size(component_thistories,1)-10:size(component_thistories,1)
    c_idx = I(sorted_idx);
    
    text(component_centroids(c_idx,2), component_centroids(c_idx,1), ...
        cellstr(num2str(sorted_idx)), 'FontSize', 15, 'Color', 'k');
    title(sprintf('distance rank=%i, centroid=%i. size = %i, location: (%i, %i,)',sorted_idx,c_idx, component_sizes(c_idx),pos(1),pos(2)));
    axis off
end


hold off;
