%===================calc_pointvortex_hamiltonian_disk.m===================%
%                                                                         %
% Calculate point vortex Hamiltonian for a given vortex configuration in  %
% a uniform, disk-shaped trap with hard walls.                            %
%                                                                         %
% Inputs:                                                                 %
% z (an N_v*1 vector of 2D vortex positions, in complex coordinates,      %
%    z = x + iy)                                                          %
% s (an N_v*1 vector of vortex circulations, (s = +/-1 for                %
%    vortices/antivortices, respectively)                                 %
%                                                                         %
% Outputs:                                                                %
% Energy of the vortex configuration                                      %
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


function H = calc_pointvortex_hamiltonian_disk(z, s)

    % Number of vortices
    N_v = length(s);

    % Constant in front of the Hamiltonian
    E_o = 1;

    % Create meshgrids from vortex location/sign vectors
    [z_x, z_y] = meshgrid(z,z);
    [s_x, s_y] = meshgrid(s,s);
    z_diff     = abs(z_x - z_y);
    
    %% 1. Interaction energy
    H_int = - 2 * E_o .* s_x .* s_y .* log(z_diff); % Construct interaction energy matrix
    H_int = triu(H_int,1);                        % Upper triangular part (j > i)
    H_int(isnan(H_int)) = 0;                      % Remove NaNs on diagonal (if they exist)
    
    %% 2. Image interaction energy
    % Terms which go into the image interaction Hamiltonian
    term_1 = 2 * (real(z_x) .* real(z_y) + imag(z_x) .* imag(z_y));
    term_2 = ( abs(z_x) .* abs(z_y) ).^2;
    
    % Other vortices' images
    H_img  = E_o * s_x .* s_y .* log(1 - term_1 + term_2);
    
    % Divide self-image terms by a factor of 2
    H_img  = triu(H_img, 0) .* (1 - 0.5.*eye(N_v)); % Pull out j >= i, and divide diagonal by 2
    
    %% Add the two contributions
    H      = sum(H_int(:) + H_img(:));

end