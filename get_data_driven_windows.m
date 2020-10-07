function [peak_locations,number_of_windows,average_window_time,window_time]=get_data_driven_windows(connectivity_complete,mean_freq,fs)
% Computes measures derived from the temporal activation sequence and the
% instantaneous correlation tensor
% Nunez et al., 2020 Abnormal meta-state activation of dynamic brain
%                    networks across the Alzheimer spectrum
%
%       Input:
%               - connectivity_complete: N x N x L weighted connectivity tensor
%               where N is the number of channels (ROIs, electrodes...) and
%               L is the number of temporal samples
%               - mean_freq: mean frequency of the frequency band of
%               interest
%               - fs: sampling frequency
%
%       Output:
%               - peak_locations: locations of the maxima in the gradient matrix of
%               the recurrence plot, corresponding to transitions between states.
%               - number_of_windows: number of data-driven windows (state
%               transitions)
%               - average_window_time: average length of a window (in
%               seconds)
%               - window_time: length of each individual window (in
%               seconds)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the upper diagonal for each temporal sample
for nWindow=1:size(connectivity_complete,3)
    temp=triu(connectivity_complete(:,:,nWindow),1);
    At=temp;
    m =(1:size(At,1)).'>=(1:size(At,2));
    meas(:,nWindow)=At(m==0);
end

% Compute recurrence plot
RP=corr(meas);

% Get data-driven windows. Based on:
% Tewarie et al., 2019 Tracking dynamic brain networks using high temporal
% resolution MEG measures of functional connectivity
% The difference is that here we base the windows on the connectivity
% itself

vert_lns = sum(abs(gradient(RP)));
vert_lns_sm = smooth(vert_lns,round(fs/mean_freq*0.2));
[~, peak_locations] = findpeaks(vert_lns_sm,'MinPeakDistance',round((fs/mean_freq)));

number_of_windows=length(peak_locations)-1;
average_window_time=mean(diff(peak_locations).*(1/fs));
window_time=diff(peak_locations).*(1/fs);

end