function [acceleration,instantaneous_correlation_acceleration,velocity,instantaneous_correlation_speed]=compute_velocity(instantaneous_correlation_tensor)
% Computes the instantanaeous acceleration and velocity
% The time unit is the sample (each point in the second dimension is a
% sample)
% The position (r) unit is the correlation in each time point. There are N dimensions to the
% correlation (as much as states)
% Nunez et al., 2020 Abnormal meta-state activation of dynamic brain
%                    networks across the Alzheimer spectrum
%
%       Input:
%               - instantaneous_correlation_tensor: K x L array indicating
%               the Spearman correlation of each community (meta-state)
%               with the instantaneous connectivity tensor for each temporal
%               sample, where K is the number of communities detected and L
%               is the number of temporal samples
%
%       Output:
%               - acceleration: K x L-2 array representing the speed of the
%               ICT for each temporal sample, where K is the number of communities 
%               detected and L is the number of temporal samples in the 
%               instantaneous connectivity tensor
%               - instantaneous_correlation_acceleration: array of length L-2 with
%               the instaneous correlation acceleration for each temporal sample,
%               whre L is the number of temporal samples in the
%               instantaneous connectivity tensor
%               - velocity: K x L-1 array representing the velocity of the
%               ICT for each temporal sample, where K is the number of communities 
%               detected and L is the number of temporal samples in the 
%               instantaneous connectivity tensor
%               - instantaneous_correlation_speed: array of length L-1 with
%               the instaneous correlation speed for each temporal sample,
%               whre L is the number of temporal samples in the
%               instantaneous connectivity tensor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

velocity=diff(instantaneous_correlation_tensor,1,2); % velocity=dr/dt
instantaneous_correlation_speed=sqrt(sum(velocity.^2,1));

acceleration=diff(velocity,1,2); % acceleration=dv/dt
instantaneous_correlation_acceleration=sqrt(sum(acceleration.^2,1));

end