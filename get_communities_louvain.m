function [communities_reshaped,temporal_activation_sequence,dwell_time_total,dwell_time_average,instantaneous_correlation_tensor,community_modularity,transition_cost] = ...
    get_communities_louvain(connectivity_subsampled,connectivity_complete,community_iterations)
% Extracts the main meta-states from the connectivity tensors
% Nunez et al., 2020 Abnormal meta-state activation of dynamic brain
%                    networks across the Alzheimer spectrum
%
%       Input:
%               - connectivity_subsampled: N x N x M subsampled instantaneous
%               connectivity tensor where N is the number of channels (ROIs,
%               electrodes...) and M is the number of data-driven windows.
%               This tensor is built by averaging "connectivity_complete"
%               in the data-driven windows
%               - connectivity_complete: N x N x L weighted connectivity tensor
%               where N is the number of channels (ROIs, electrodes...) and
%               L is the number of temporal samples
%               - community_iterations: number of times to run the Louvain
%               algorithm for  community detection
%
%       Output:
%               - communities_reshaped: K x N x N matrix, where K is the
%               number of communities detected and N is the number of channels
%               - temporal_activation_sequence: array of size M indicating
%               the active state in each temporal sample, where M is the number
%               of temporal samples. Each state corresponds to a
%               community in the "communities_reshaped" matrix. E.g. state
%               number "2" corresponds to communities_reshaped(2,:,:)
%               - dwell_time_total: total dwell time (in samples) of each
%               community (meta-state) in the whole sequence
%               - dwell_time_average: average dwell time (in samples) of each
%               community (meta-state)
%               - instantaneous_correlation_tensor: K x L array indicating
%               the Spearman correlation of each community (meta-state)
%               with the instantaneous connectivity tensor for each temporal
%               sample, where K is the number of communities detected and L
%               is the number of temporal samples
%               - community_modularity: modularity value of the solution
%               with the highest modularity of the Louvain GJA algorithm
%               - transition_cost: representation of the metabolic cost of
%               switching between one connectivity pattern to the next
%               (distance between one instantaneous connectivity pattern
%               and the next). Needed for the leap size
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the number of channels
nChannels=size(connectivity_subsampled,1);

% Get the upper diagonal for each subsampled window. Remove the diagonal
% since it is always 0 or 1
for nWindow=1:size(connectivity_subsampled,3)
    temp=triu(connectivity_subsampled(:,:,nWindow),1);
    At=temp;
    m =(1:size(At,1)).'>=(1:size(At,2));
    connectivity_subsampled_reshaped_triu(:,nWindow)=At(m==0);
end

% Get the upper diagonal for each temporal sample
for nWindow=1:size(connectivity_complete,3)
    temp=triu(connectivity_complete(:,:,nWindow),1);
    
    At=temp;
    m =(1:size(At,1)).'>=(1:size(At,2));
    
    connectivity_reshaped_triu(:,nWindow)=At(m==0);
end

% Compute the recurrence plots (graphs) for the windowed connectivity
% matrix
connectivity_subsampled_reshaped=reshape(squeeze(connectivity_subsampled),[size(squeeze(connectivity_subsampled),1)*size(squeeze(connectivity_subsampled),2) size(squeeze(connectivity_subsampled),3)]);

RP_subsampled=corr(connectivity_subsampled_reshaped_triu,'Type','Spearman');
RP=corr(connectivity_reshaped_triu,'Type','Spearman');

transition_cost=1-diag(RP,1);

% Since Louvain is non-deterministic, run the algorithm a set number of times times and keep
% the solution with the highest modularity
community_modularity_temp=0;

for iterationN=1:community_iterations
    % This function is available in the Brain Connectivity Toolbox: https://sites.google.com/site/bctnet/
    [community,community_modularity] = community_louvain(RP_subsampled,[],[],'negative_asym');
    
    if community_modularity>community_modularity_temp
        community_modularity_temp=community_modularity;
        community_final=community;
    end
    
end

community_modularity=community_modularity_temp;

community=community_final;

% Get the centroids from each community by computing the mean of all points
% in that community and rebuild them
for clusterN=1:length(unique(community))
    
    clusters(clusterN,:)=mean(connectivity_subsampled_reshaped(:,community==clusterN),2);
    clusters_triu(clusterN,:)=mean(connectivity_subsampled_reshaped_triu(:,community==clusterN),2);
    
end

% Assign each sample to a cluster.
[~,temporal_activation_sequence] = pdist2(clusters_triu,connectivity_reshaped_triu.','spearman','Smallest',1);

% Compute and save the correlations with all centroids
D=NaN(size(clusters_triu,1),size(connectivity_reshaped_triu,2));

for ii=1:size(clusters_triu,1)
    D(ii,:)=pdist2(clusters_triu(ii,:),connectivity_reshaped_triu(:,:).','spearman','Smallest',1);
end

instantaneous_correlation_tensor=1-D;

% Back to normal
temporal_activation_sequence=temporal_activation_sequence.';

% Get the total dwell time of each state
for clusterN=1:length(unique(community))
    
    dwell_time(clusterN)=sum(temporal_activation_sequence==clusterN);
    
end

% Get the average dwell time of each state
communities=1:length(unique(community));
for clusterN=1:length(unique(community))
    
    activated_time = (temporal_activation_sequence.'==communities(clusterN)); % 1 if the state is active, 0 else
    nStates = countNumberOfSequences(activated_time); % Number of times the state appears in the sequence
    dwell_time_average(clusterN)=dwell_time(clusterN)/nStates; % Total dwell time / times the state appears
    
end

% Order the clusters in order of importance (more dwell time to less dwell time)
[dwell_time_total,I]=sort(dwell_time,'descend');
clusters_sorted=clusters(I,:);
dwell_time_average=dwell_time_average(I);

% Change temporal activation to reflect the new state order
temporal_activation_sequence_temp=zeros(1,length(temporal_activation_sequence));

for nState=1:length(I)
    
    indices=temporal_activation_sequence==I(nState);
    temporal_activation_sequence_temp(indices)=nState;
    
end

temporal_activation_sequence=temporal_activation_sequence_temp; 
% Now the temporal_activation_sequence has the states in the same order as "clusters_sorted" (in order of dwell time)

% Rebuild the clusters into connectivity matrixes
for nCluster=1:size(clusters_sorted, 1)
    
    communities_reshaped(nCluster,:,:)=reshape(clusters_sorted(nCluster,:),[nChannels nChannels]);
    
end

end