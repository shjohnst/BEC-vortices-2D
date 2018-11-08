%==========================MCMC_pointvortices.m===========================%
%                                                                         %
% Matlab code to perform Markov Chain Monte Carlo simulations for point   %
% vortices in a uniform disk configuration. Runs through a chosen number  %
% of temperature sweeps, starting at extreme negative temperatures        %
% (cluster dominated), and moving towards extreme positive temperatures   %
% (dipole dominated).                                                     %
%                                                                         %
% Requirements: 'calc_pointvortex_hamiltonian_disk.m'                     %
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


clear, clc, close all

%% System parameters
N_v_i   = 20; % Number of vortices
R       = 1;  % System radius

% Minimum distance between vortices
core_size = 0.006 * R;

% Constant in front of point-vortex Hamiltonian
E_o = 1;

%% MC sampling settings
n_save  = 1e4; % How often to save configuration (per temperature, per sweep)
n_sweep = 2e2; % Number of temperature sweeps
n_steps = 1e5; % Number of steps (per temperature, per sweep)
n_temp  = 100; % Number of temperature points, either side of beta = 0

%% Set inverse temperature range to sample
beta_max     = 2.5;
beta_min     = -1.5*4/N_v_i;
beta_vec_neg = linspace(beta_min,0,n_temp+1);
beta_vec_pos = linspace(0, beta_max, n_temp+1);

beta_vec = [beta_vec_neg, beta_vec_pos(2:end)]/E_o;

%% Initiate output vectors - fill with NaN
p_v    = NaN * zeros(n_steps/n_save*n_sweep, length(beta_vec));
PVE_v  = p_v;
PVL_v  = p_v;

Z      = (NaN + NaN*1j) * zeros(n_steps/n_save*n_sweep, N_v_i, length(beta_vec));
S      = Z;

% Initiate sweep counting index
nn = 1;

%% Loop over temperature sweeps
for sweep = 1:n_sweep
    
    % Print sweep number
    disp('=====================================================')
    disp(['SWEEP NO. ', num2str(sweep), ' OF ', num2str(n_sweep)])
    disp('=====================================================')
    
    % Each sweep starts at extreme negative temperatures, so input with a highly clustered state
    z_v = [0.3 * rand(1,ceil(N_v_i/2)) .* exp(1i * 2*pi * rand(1, ceil(N_v_i/2))) - 0.5, 0.3 * rand(1,floor(N_v_i/2)) .* exp(1i * 2*pi * rand(1, floor(N_v_i/2))) + 0.5];
    s_v = [ones(1,ceil(N_v_i/2)), -ones(1,floor(N_v_i/2))];

    % Temperature counting index
    mm = 1;
    
    tic

    %% Loop over temperatures in this sweep
    for beta = beta_vec

        % Initial step size
        dr = 0.08 * 2;

        % Clear all variables from previous iteration
        clearvars -except N_v_i R E_o n_save n_sweep n_steps nn mm sweep dr n_trav beta image_show col_pos col_neg z_v s_v core_size theta x_circ y_circ config_fig beta_vec save_fig Z S p_v dipole PVE_v PVL_v

        % Print current temperature
        disp('-----------------------------------------------------')
        disp(['beta = ', num2str(beta*E_o)])
        disp('-----------------------------------------------------')

        %% Initiate loop

        % Current energy, and weighted probability
        PVE = calc_pointvortex_hamiltonian_disk(z_v, s_v);
        p   = exp(-PVE*beta);

        % Initiate rejection counter
        n_reject = 0;

        jj = 1;
        kk = 1;

        while jj <= n_steps

            z_v_prev   = z_v;
            PVE_prev   = PVE;
            p_prev     = p;

            %% Shift a random vortex in a random direction
            % Choose which vortex to shift
            ind = randi(N_v_i);
            
            % Unphysically large values which will ensure the 'while' loop runs
            new_z_v      = 10*R;
            min_distance = core_size/10;

            % Shift the chosen vortex, but make sure it doesn't get too
            % close to any other vortex or the boundary.
            while abs(new_z_v) > (R - core_size) || min_distance < core_size
                % Choose displacement vector, shift chosen vortex
                dR      = dr * exp(1i * 2*pi * rand);
                new_z_v = z_v(ind) + dR;

                % Minimum distance between new vortex and all other vortices
                min_distance = min(abs(z_v_prev(1:end ~= ind) - new_z_v));
            end

            % Replace old vortex position with the new one
            z_v(ind) = new_z_v;

            %% Calculate energy and change in energy / weighting
            % Calculate energy of new state
            PVE = calc_pointvortex_hamiltonian_disk(z_v, s_v);

            % Calculate change in energy and weighting compared to previous state
            dE = PVE - PVE_prev;
            dp = exp(-dE*beta);

            % Current weighting
            p  = exp(-PVE*beta);

            %% Choose whether to accept/reject new configuration based on Markov Chain Monte Carlo rule
            if (dp < 1) && (rand > dp)
                %% Case 1: Reject new vortex configuration
                
                % Reset vortex positions, energy, and weighting
                z_v    = z_v_prev;
                PVE    = PVE_prev;
                p      = p_prev;

                % Add 1 to the rejection counter
                n_reject = n_reject + 1;

                % If too many consecutive rejections occur, halve the step size
                if n_reject > 100
                    disp('Too many rejections. Step size halved.')
                    disp(' ')
                    dr = 0.5 * dr;

                    % Reset rejection counter
                    n_reject = 0;
                end

            else
                %% Case 2: Accept new vortex configuration
                n_reject = 0;

                % Calculate and save observables. Only save the second half
                % of the chain to remove burn-in time
                if mod(jj, n_save/2) == 0 && jj > n_steps/2

                    % Save current vortex locations and charges
                    Z(kk + (nn-1)*n_steps/n_save,:,mm) = z_v;
                    S(kk + (nn-1)*n_steps/n_save,:,mm) = s_v;

                    % Save current point-vortex energy and momentum
                    PVE_v(kk + (nn-1)*n_steps/n_save,mm)  = calc_pointvortex_hamiltonian_disk(z_v, s_v);
                    PVL_v(kk + (nn-1)*n_steps/n_save,mm)  = sum(s_v .* (abs(z_v)).^2);

                    % Save current weighting
                    p_v(kk + (nn-1)*n_steps/n_save,:,mm) = p;

                    % Display current progress through this temperature
                    disp(['Step no. ', num2str(jj), ' (', num2str(100*jj/n_steps), '%)'])
                    disp(' ')

                    kk = kk + 1;
                    
                end

                jj = jj + 1;

            end

        end

        mm = mm + 1;

    end
    
    % Print time taken for each temperature sweep
    t_sim = toc;
    disp(['Time taken for temperature sweep: ', num2str(t_sim), ' seconds'])
   
    nn = nn + 1;
    
end

% Create folder if it doesn't already exist
dir_name = 'data_beta';
if not(exist(dir_name, 'dir') == 7)
    mkdir(dir_name);
end

% Save data to output file
save([dir_name, '/N_', num2str(N_v_i), '_sweep_data.mat'], '-v7.3');