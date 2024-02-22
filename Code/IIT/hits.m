%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHITOOLBOX FROM OIZUMI
% This Toolbox allows to compute phi G ; phi* ; phi H and mutual information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

tic

toolboxpath = "/d/aveyrie/Documents/MATLAB/phi_toolbox_oizumi";
npymatlabpath = "/d/aveyrie/Documents/MATLAB/npy-matlab/npy-matlab"
datapath = "/d/aveyrie/Documents/epochs_hits_around_detection/125Hz/"
addpath(genpath(toolboxpath));
addpath(genpath(npymatlabpath));
addpath(genpath(datapath));

epoch = readNPY(datapath + "epochs_hits_around_detection_values_sujet_xx.npy");

sf = 125;
epoch_time_length = 3; % second
time_sample = sf * epoch_time_length;
time_around_detection = time_sample * 2;
time_lag = time_around_detection - 100; % tau (millisecond)

epoch_shape = size(epoch);
trial_number  = epoch_shape(1);
channel_number  = epoch_shape(2);
time_around_detection = epoch_shape(3);

% target_channels = [1, 35, 20, 49]; % Fz, FCz, Cz, CPz % Sagital Cluster
% target_channels = [6, 36, 25, 37, 21, 54, 10, 50, 20]; % 'FC1','FCz','FC2','C1','Cz','C2','CP1','CPz','CP2'
% target_channels = [33, 55, 7, 22, 38, 51]; % FT7, FT8, T7, T8, TP7, TP8 % Temporal Cluster
% target_channels = [33, 7, 38]; % FT7, T7, TP7 % Left hemisphere
% target_channels = [55, 22, 51]; % FT8, T8, TP8 % Right hemisphere
target_channels = [1, 7, 20, 22, 33, 35, 38, 49, 51, 55]; % Fz, FCz, Cz, CPz, FT7, FT8, T7, T8, TP7, TP8 % Temporo-Sagital Cluster
target_cluster = epoch(:,target_channels,:);

phi_G = zeros(trial_number, time_lag);
phi_star = zeros(trial_number, time_lag);
phi_H = zeros(trial_number, time_lag);
MI = zeros(trial_number, time_lag);

phi_G_mean = 1:time_lag;
phi_star_mean = 1:time_lag;
phi_H_mean = 1:time_lag;
MI_mean = 1:time_lag;

for trial = 1:1:trial_number;

	tmp = target_cluster(trial,:,:);
	tmp = fliplr(tmp);
	N = length(tmp(:,1,:));
	tmp = reshape(tmp, [length(target_channels), N]);
	options.type_of_dist = 'Gauss';
	options.type_of_MIPsearch = 'Queyranne';

	for j=1:1:time_lag;

	tau = j;
	params.tau = tau;

	options.type_of_phi = 'Geo';
	[Z_MIP_G, phi_G(trial,j)] = MIP_search(tmp, params, options)
	phi_G(trial,j) = phi_comp(tmp, Z_MIP_G, params, options)

	options.type_of_phi = 'star';
	[Z_MIP_star, phi_star(trial,j)] = MIP_search(tmp, params, options)
	phi_star(trial,j) = phi_comp(tmp, Z_MIP_star, params, options)

	options.type_of_phi = 'SI';
	[Z_MIP_H, phi_H(trial,j)] = MIP_search(tmp, params, options)
	phi_H(trial,j) = phi_comp(tmp, Z_MIP_H, params, options)

	options.type_of_phi = 'MI';
	[Z_MIP_MI, phi_MI(trial,j)] = MIP_search(tmp, params, options)
	MI(trial,j) = phi_comp(tmp, Z_MIP_MI, params, options)

	end
	
end

for lag = 1:1:time_lag;
	phi_G_mean(lag) = mean(phi_G(:,lag));
	phi_star_mean(lag) = mean(phi_star(:,lag));
	phi_H_mean(lag) = mean(phi_H(:,lag));
	MI_mean(lag) = mean(MI(:,lag));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = cat(1, phi_G_mean, phi_star_mean, phi_H_mean, MI_mean);
chdir('/d/aveyrie/Documents/MATLAB/')
csvwrite('phi_measures_mean_over_trials_hits_sujet_xx.csv', phi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure;
%hold on
%plot(phi_G_mean, 'm');
%plot(phi_star_mean, 'b');
%plot(phi_H_mean, 'g')
%plot(MI_mean, 'r')
%title('Evolution of Phi with time lag Tau')
%xlabel('Tau');
%ylabel('Phi');
%legend('Phi_G', 'Phi *', 'Phi_H', 'MI')
%save('/d/aveyrie/Documents/MATLAB/phi_measures_mean_over_trials_hits_sujet_xx.png')

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
