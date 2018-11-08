%==============================analyse_data.m=============================%
%                                                                         %
% Matlab code to analyse vortex configurations output from Markov Chain   %
% Monte Carlo code, 'MCMC_pointvortices.m'.                               %
%                                                                         %
% Classifies vortices using prescription outlined in Section 6 of Valani  %
% et al., 'Einstein-Bose condensation of Onsager vortices', New Journal   %
% of Physics 20, 053038 (2018). Then plots the mean populations of        %
% clustered vortices, dipole vortices and free vortices as a function of  %
% inverse temperature.                                                    %
%                                                                         %
% Copyright (C) 2018  Andrew J. Groszek                                   %
%                                                                         %
% This program is free software: you can redistribute it and/or modify    %
% it under the terms of the GNU General Public License as published by    %
% the Free Software Foundation, either version 3 of the License, or       %
% (at your option) any later version.                                     %
%                                                                         %
% This program is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           %
% GNU General Public License for more details.                            %
%                                                                         %
% You should have received a copy of the GNU General Public License       %
% along with this program.  If not, see <https://www.gnu.org/licenses/>.  %
%                                                                         %
% Contact: Andrew.Groszek@newcastle.ac.uk                                 %
%                                                                         %
%=========================================================================%

close all, clear, clc

% Plot classified vortices?
plot_clusters = 0;

% Number of vortices
N_v = 12;

% Load output file from 'MCMC_pointvortices.m'
dir_name = 'data_beta';
load([dir_name, '/N_', num2str(N_v), '_sweep_data.mat'])

% Initiate temperature counting index
mm = 1;

%% Loop over all sampled temperatures
for beta = beta_vec
    
    % Display current temperature
    disp(['beta = ', num2str(beta*E_o)])
    
    % Initiate observable output vectors for this temperature
    N_clus_v = zeros(1,size(Z,1)); % Number of clustered vortices
    N_dip_v  = N_clus_v;           % Number of dipole vortices
    N_free_v = N_clus_v;           % Number of free vortices
    n_clus_v = N_clus_v;           % Total number of clusters
    N_vc_v   = N_clus_v;           % Mean number of vortices per cluster
    
    % Extract this temperature's list of vortex configurations
    Z_mm = squeeze(Z(:,:,mm));
    S_mm = squeeze(S(:,:,mm));
    
    %% Loop over all saved configurations at this temperature
    % Choose which configurations to analyse
    ii_vec   = 1:size(Z,1);
        
    % Initial configuration counting index
    kk = 1;
    
    % Begin loop
    for ii = ii_vec
        
        % Pull out this configuration
        z_i = Z_mm(ii,:);
        s_i = S_mm(ii,:);

        %% Classify vortices
        c_v_i = classify_vortices(z_i, s_i, plot_clusters, 1);
        
        % Plot classified vortex configuration?
        if plot_clusters == 1
            drawnow
            pause(0.1)
        end

        % Locate clustered vortices, dipoles and free vortices
        ind_clus = find(abs(real(c_v_i))>0);
        ind_dip  = find(imag(c_v_i)>0);
        ind_free = find(isnan(c_v_i));
        
        %% Save vortex observables
        % Number of clustered vortices, dipole vortices, and free vortices
        N_clus_v(kk) = length(ind_clus);
        N_dip_v(kk)  = length(ind_dip);
        N_free_v(kk) = length(ind_free);
        
        kk = kk + 1;
        
    end
    
    %% Save all data for this temperature
    N_clus_mat(mm,:) = N_clus_v;
    N_dip_mat(mm,:)  = N_dip_v;
    N_free_mat(mm,:) = N_free_v;
    
    mm = mm + 1;

end

%% Plot thermometry curves for this vortex number
% Scale temperature with the 'EBC' and 'BKT' critical temperatures in
% negative/positive regions, respectively
beta_vec_scaled                    = beta_vec * E_o;
beta_vec_scaled(beta_vec_scaled<0) = beta_vec_scaled(beta_vec_scaled<0)/4 *N_v_i;
beta_vec_scaled(beta_vec_scaled>0) = beta_vec_scaled(beta_vec_scaled>0)/2;

% Plot mean populations of clustered vortices, dipole vortices and free
% vortices as a function of the scaled inverse temperature
figure
plot(beta_vec_scaled, mean(N_clus_mat,2) / N_v_i)
hold on
plot(beta_vec_scaled, mean(N_dip_mat,2)  / N_v_i)
plot(beta_vec_scaled, mean(N_free_mat,2) / N_v_i)
set(gca, 'XDir', 'Reverse')
xlabel('$\beta$', 'interpreter', 'latex'), ylabel('$p_j$', 'interpreter', 'latex')