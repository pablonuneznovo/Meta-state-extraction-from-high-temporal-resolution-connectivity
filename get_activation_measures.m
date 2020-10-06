function [TAS_complexity,leap_size,instantaneous_correlation_speed,instantaneous_correlation_acceleration] = get_activation_measures(temporal_activation_sequence,instantaneous_correlation_tensor,transition_cost)
% Computes measures derived from the temporal activation sequence and the
% instantaneous correlation tensor
% Nunez et al., 2020 Abnormal meta-state activation of dynamic brain
%                    networks across the Alzheimer spectrum
%
%       Input:
%               - temporal_activation_sequence: array of size M indicating
%               the active state in each temporal sample, where M is the number
%               of temporal samples.
%               - instantaneous_correlation_tensor: K x L array indicating
%               the Spearman correlation of each community (meta-state)
%               with the instantaneous connectivity tensor for each temporal
%               sample, where K is the number of communities detected and L
%               is the number of temporal samples
%               - transition_cost: L-1 array representating of the metabolic cost of
%               switching between one connectivity pattern to the next
%               (distance between one instantaneous connectivity pattern
%               and the next), where L is the number of temporal samples in
%               the instantaneous connectivity tensor
%
%       Output:
%               - TAS_complexity: K x N x N matrix, where K is the
%               number of communities detected and N is the number of channels
%               - leap_size: array with the transition cost for each state
%               transition
%               - instantaneous_correlation_speed: array of length L-1 with
%               the instaneous correlation speed for each temporal sample,
%               whre L is the number of temporal samples in the
%               instantaneous connectivity tensor
%               - instantaneous_correlation_acceleration: array of length L-2 with
%               the instaneous correlation acceleration for each temporal sample,
%               whre L is the number of temporal samples in the
%               instantaneous connectivity tensor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check when a state transition occurs
diffTemp=diff(temporal_activation_sequence);

% The leap size is the transition cost when a state transition occurs
leap_size=transition_cost(diffTemp~=0);

% Compute Lempel-Ziv complexity of the temporal activation sequence
TAS_complexity = compute_LZC(temporal_activation_sequence);

% Speed and acceleration
% Compute average acceleration and speed
[~,instantaneous_correlation_acceleration,~,instantaneous_correlation_speed]=compute_velocity(instantaneous_correlation_tensor);

end
