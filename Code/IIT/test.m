%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHITOOLBOX FROM OIZUMI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

addpath(genpath('/d/aveyrie/Documents/TII/phi_toolbox_last_version'));
addpath(genpath('/d/aveyrie/Documents/125Hz/after/sujet_1'));

epoch = csvread('epoch_0.csv');
epoch = epoch(1:16,:);

N = length(epoch(:,1));
O = 100; % time lag / observations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi*

tic
options.type_of_dist = 'Gauss';
options.type_of_MIPsearch = 'Queyranne'
options.type_of_phi = 'star';
tau = 1;
params.tau = tau;
phi_star = 1:O;
for i=1:1:O;
    tau = i;
    params.tau = tau;
    [Z_MIP, phi_star(i)] = MIP_search(epoch, params, options)
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mutual information

tic
options.type_of_dist = 'Gauss';
options.type_of_MIPsearch = 'Queyranne'
options.type_of_phi = 'MI';
tau = 1;
params.tau = tau;
MI = 1:O;
for i=1:1:O;
    tau = i;
    params.tau = tau;
    [Z_MIP, MI(i)] = MIP_search(epoch, params, options)
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometric phi

tic
options.type_of_dist = 'Gauss';
options.type_of_MIPsearch = 'Queyranne'
options.type_of_phi = 'Geo';
tau = 1;
params.tau = tau;
phi_G = 1:O;
for i=1:1:O;
    tau = i;
    params.tau = tau;
    [Z_MIP, phi_G(i)] = MIP_search(epoch, params, options)
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stochastic interaction

tic
options.type_of_dist = 'Gauss';
options.type_of_MIPsearch = 'Queyranne'
options.type_of_phi = 'SI';
tau = 1;
params.tau = tau;
phi_H = 1:O;
for i=1:1:O;
    tau = i;
    params.tau = tau;
    [Z_MIP, phi_H(i)] = MIP_search(epoch, params, options)
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = cat(1, phi_star, phi_G, phi_H, MI)
chdir('/d/aveyrie/Documents')
csvwrite('phi.csv', phi)

% figure;
% plot(phi_star);
% title('Evolution of Phi* with time lag Tau')
% xlabel('Tau');
% ylabel('Phi*');

% probs = Cov_comp(epoch, 1, true);
% Cov_X = probs.Cov_X;
% Cov_Y = probs.Cov_Y;
% Cov_XY = probs.Cov_XY;
% Cov_X_Y = Cov_cond(Cov_X, Cov_XY, Cov_Y);
% probs = data_to_probs(epoch, params, options);
% phi = phi_comp_probs(options.type_of_dist, options.type_of_phi, Z, probs);
% phi = phi_Gauss(type_of_phi, Z, probs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
