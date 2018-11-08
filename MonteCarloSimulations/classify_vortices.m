%===========================classify_vortices.m===========================%
%                                                                         %
% Vortex classification algorithm. Input a vortex configuration with      %
% locations z_v and signs s_v, and the algorithm will detect vortex       %
% clusters, dipoles and free vortices.                                    %
% The output tells you which cluster (if any) each vortex belongs to.     %
%                                                                         %
% Copyright (C) 2016  Andrew J. Groszek                                   %
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



function c_v = classify_vortices(z_v, s_v, plot_bool, fig_handle)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vortex classification algorithm. Input a vortex configuration with locations z_v and
    % signs s_v, and the algorithm will detect vortex clusters, dipoles and free vortices.
    % The output tells you which cluster (if any) each vortex belongs to.
    %
    % For details, see Section 6 of Valani et al., 'Einstein-Bose condensation of Onsager vortices',
    % New Journal of Physics 20, 053038 (2018).
    %
    %
    % Author: Andrew Groszek
    % Date: 9/3/16
    %         Edit: Fixed issue where having only one vortex would result in an error
    %         Edit: Fixed issue where all vortices being the same sign would give an error (19/7/16)
    %         Edit: Changed code to take figure handle as input (20/7/16)
    %
    % Inputs:
    % z_v:        A vector of complex co-ordinates corresponding to 2D vortex locations
    %             NOTE: z_v needs to be normalised such that |z_v| < 1
    % s_v:        A vector of integers corresponding to vortex charges (usually +/-1)
    % plot_bool:  Boolean variable - produce plots (true/1) or not (false/0)?
    % fig_handle: The handle of the figure you want to plot in
    %
    % Outputs:
    % c_v: a vector of labels corresponding to each vortex's cluster status.
    %      - Positive clusters are labelled as +1, +2, +3, ...
    %      - Negative clusters are labelled as -1, -2, -2, ...
    %      - Dipoles are labelled as 1i, 2i, 3i, ...
    %      - Vortices that aren't part of a cluster or dipole are labelled NaN
    %        e.g. if c_v = [+1, +1, 1i, 1i, +1, -1, 2i, -1, 2i, +2, +2], then:
    %        - Vortices 1, 2 and 5 are in positive cluster #1
    %        - Vortices 10 and 11 are in positive cluster #2
    %        - Vortices 6 and 8 are in negative cluster #1
    %        - Vortices 3 and 4 form dipole pair #1
    %        - Vortices 7 and 9 form dipole pair #2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Initiate variables
    
    % Number of vortices
    N_v = length(z_v);
    
    % If z_v isn't normalised
    if max(abs(z_v)>1)
        warning('z_v not normalised!')
    end

    % Sort z_v and s_v so that all the + vortices are on the left, all - on the right
    [s_v, ind_sort] = sort(s_v,'descend');
    z_v = z_v(ind_sort);

    z_pos = z_v(s_v > 0);
    z_neg = z_v(s_v < 0);

    % If there are fewer than two vortices, don't perform cluster detection
    if N_v < 2
        c_v           = nan * ones(1,N_v); % Label the one vortex as a free vortex
        dipole        = [];
        
        if plot_bool == true
            cluster_plot(z_v, c_v, dipole, z_pos, z_neg, fig_handle)
        end
        
        return
    % If all vortices are the same sign, immediately put them all in a single cluster
    elseif isempty(z_pos) || isempty(z_neg)
        c_v   = s_v(1) * ones(1,N_v); % Label all vortices as belonging to cluster +/-1
        
        n_pos = length(z_pos);
        n_neg = length(z_neg);
        
        dipole = [];
        
        if plot_bool == true
            cluster_plot(z_v, c_v, dipole, z_pos, z_neg, fig_handle)
        end
        
        return
    end

    % Define labelling vectors for each vortex
    c_v = NaN * ones(1,N_v); % Which cluster does vortex i belong to?

    [z_X, z_Y] = meshgrid(z_v);

    % Calculate |z_i - z_j|
    z_diff = (abs(z_X - z_Y));

    % % Set diagonal and lower triangle to infinity so they're not counted
    z_diff(z_diff==0) = inf; % Set diagonal to infinity

    %% Calculate the order of vortex distances j from each vortex i

    % Number of positive vortices
    n_pos  = sum(s_v>0);
    n_neg  = sum(s_v<0);
    n_diff = n_pos - n_neg;

    % Find the order of distances from each vortex to each positive (negative) vortex
    [z_diff_pos_sort, ind_sort_pos] = sort(z_diff(:,s_v>0),2);
    [z_diff_neg_sort, ind_sort_neg] = sort(z_diff(:,s_v<0),2);

    % The indices of the negative vortices need to be offset by n_pos (because s_v = [1, 1, ..., 1, -1, -1, ..., -1])
    ind_sort_neg = ind_sort_neg + n_pos;

    % If there are more of one sign than the other, pad the smaller matrix with
    % infinities so they can still be compared
    if n_pos > n_neg
        z_diff_neg_sort = [z_diff_neg_sort, inf .* ones(N_v,abs(n_diff))];
    elseif n_neg > n_pos
        z_diff_pos_sort = [z_diff_pos_sort, inf .* ones(N_v,abs(n_diff))];
    end

    %% Find dipoles

    % Find which vortices have nearest neighbours which are the opposite sign
    list = (1:N_v).';
    v_p = z_diff_pos_sort(:,1) < z_diff_neg_sort(:,1);
    v_n = z_diff_neg_sort(:,1) < z_diff_pos_sort(:,1);

    % Concatenate list of vortices with sorted indices, pulling out only those
    % for which the closest vortex is the opposite sign of the chosen vortex.
    % e.g.:
    % M_pn = [1, 0]  (closest vortex same sign)
    %        [2, 18] (closest vortex is z_v(18), and is opposite sign)
    %        [3, 0]  ... etc
    M_pn = [list, v_p .* ind_sort_pos(:,1)];
    M_np = [list, v_n .* ind_sort_neg(:,1)];

    M_np = M_np(1:n_pos,:);
    M_pn = M_pn(n_pos+1:end,:);

    % Check if vortices are mutual nearest neighbours
    dipole_count = 1;

    for i = 1:n_pos

        j = M_np(i,2);

        if j > 0

            % If f(i) = j and f^-1(j) = i, then i and j are mutual nearest neighbours
            if M_pn(j-n_pos,2) == i

                % Label these vortices as a dipole, increase dipole counter
                c_v([i,j])   = dipole_count * 1i;
                
                % Save vortex indices (for plotting purposes)
                dipole(dipole_count,:) = [i, j];
                
                dipole_count = dipole_count + 1;

            end

        end

    end
    
    % If there were no dipoles detected, create an empty vector
    if dipole_count == 1
        dipole = [];
    end

    %% Find n same-sign vortices closer than an opposite sign vortex

    % Loop over each vortex to check its neighbours
    for i = 1:N_v

        if (i == 1) || (i == (n_pos + 1))
            cluster_count = 1;
        end

        if i <= n_pos
            % How many + vortices are closer to this + vortex than the nearest -?
            n_closer = sum(z_diff_pos_sort(i,:) < min(z_diff_neg_sort(i,:)));
                                                 % ^ may be problematic if all inf
        else
            % How many - vortices are closer to this - vortex than the nearest +?
            n_closer = sum(z_diff_neg_sort(i,:) < min(z_diff_pos_sort(i,:)));
        end

        % Generate lists of nearest like-sign neighbours to each vortex
        if n_closer > 0
            if i <= n_pos
                list_i = ind_sort_pos(i,1:n_closer);
            else
                list_i = ind_sort_neg(i,1:n_closer);
            end
        else
            list_i = [];
        end

        % Label lists as positive/negative accordingly
        if i <= n_pos
            list_cell_pos{i} = list_i;
        else
            list_cell_neg{i} = list_i;
        end

        % Check members of this list to see if there are any mutual-sign neighbours
        for j = list_i

            % If this vortex has already been checked
            if j < i

                % If vortex i is in vortex j's list, then they are mutual neighbours
                if i <= n_pos
                    check = sum(ismember(list_cell_pos{j},i));
                else
                    check = sum(ismember(list_cell_neg{j},i));
                end

                if check

                    % If both i and j are already labelled, and their labels don't match,
                    % set all vortices in the two clusters to be in the same cluster.
                    % Choose the smaller (absolute) label.
                    if ~isnan(c_v(j)) && ~isnan(c_v(i)) && c_v(j) ~= c_v(i)
                        if abs(c_v(i)) < abs(c_v(j))
                            c_v(c_v==c_v(j)) = c_v(i);
                        else
                            c_v(c_v==c_v(i)) = c_v(j);
                        end
                    end

                    % If vortex j has already been labelled and i hasn't, label i as the
                    % same, and vice versa. Else, label them both something new.
                    if ~isnan(c_v(j)) && isnan(c_v(i))
                        c_v(i) = c_v(j);
                    elseif ~isnan(c_v(i)) && isnan(c_v(j))
                        c_v(j) = c_v(i);
                    elseif isnan(c_v(i)) && isnan(c_v(j))
                        c_v([i,j]) = cluster_count * s_v(i);
                        cluster_count = cluster_count + 1;
                    end

                end  
            end

        end

    end
    
    %% Fix labels to avoid skipping numbers
    
    % List of unique cluster indices, e.g. [1, 2, 5, 6, 7, 10, ...]
    cluster_list_pos = unique(sort( real(c_v(c_v>0))));
    cluster_list_neg = unique(sort(-real(c_v(c_v<0))));
    
    % Find the range of cluster labels used
    labelled_num_pos = max(cluster_list_pos) - min(cluster_list_pos) + 1;
    labelled_num_neg = max(cluster_list_neg) - min(cluster_list_neg) + 1;
    
    % Find actual number of clusters of each sign (i.e. number of indices used)
    actual_num_pos   = length(cluster_list_pos);
    actual_num_neg   = length(cluster_list_neg);
    
    if labelled_num_pos ~= actual_num_pos
        % Relabel clusters so that no numbers are skipped
        for kk = 1:actual_num_pos

            ind_label_old             = cluster_list_pos(kk);
            ind_label_new             = kk;
            c_v(c_v == ind_label_old) = ind_label_new;

        end
    end
    
    if labelled_num_neg ~= actual_num_neg
        % Relabel clusters so that no numbers are skipped
        for kk = 1:actual_num_neg

            ind_label_old              = cluster_list_neg(kk);
            ind_label_new              = -kk;
            c_v(c_v == -ind_label_old) = ind_label_new;

        end
    end
    
    %% Plot
    if plot_bool == true
        cluster_plot(z_v, c_v, dipole, z_pos, z_neg, fig_handle)
    end
    
    % If data was initially unsorted, sort back to original format
    [~, ind_sort_sort] = sort(ind_sort);
    c_v = c_v(ind_sort_sort);

end

%% Plot clusters. Using minimum spanning tree algorithm to find joining lines
function cluster_plot(z_v, c_v, dipole, z_pos, z_neg, fig_handle)

    % Plotting colours
    col_pos = [0, 80, 180]/255;  % Blue
    col_neg = [80, 200, 0]/255;  % Green
    col_dip = [230, 40, 20]/255; % Red
    
    lw = 1.5;

    % Initialise figure
    figure(fig_handle)
    clf
    set(gca, 'fontsize',12)
    set(gcf, 'position',[500,200,600,600],'color','w')
    axis off
    hold on

    %% Plot boundary of system
    theta  = linspace(0, 2*pi, 10000);
    x_circ = cos(theta);
    y_circ = sin(theta);
    plot(x_circ, y_circ, 'k','linewidth',1.5)
    axis(1.2*[-1 1 -1 1])
    axis square
    
    %% Plot dipoles
    for i = 1:size(dipole,1)

        x12 = [real(z_v(dipole(i,1))), real(z_v(dipole(i,2)))];
        y12 = [imag(z_v(dipole(i,1))), imag(z_v(dipole(i,2)))];

        plot(x12, y12, 'color', col_dip,'linewidth',lw)

    end
    
    %% Plot clusters
    % Number of positive/negative clusters
    max_cluster_pos =  max(real(c_v));
    max_cluster_neg = -min(real(c_v));

    % Plot lines between vortices in positive clusters
    for ii = 1:max_cluster_pos
        
        % List of vortices in this cluster
        z_v_ii = (z_v(c_v == ii)).';
            
        % If there are >2 vortices and they don't have the same x/y positions*, calculate minimum spanning tree
        % *This is a problem if the vortices exist on a numerical grid - the edge coords can't be defined
        if length(z_v_ii) > 2 && range(real(z_v_ii)) > 0 && range(imag(z_v_ii)) > 0
            
            % Calculate Delaunay triangulation between vortices in this cluster
            DT = delaunayTriangulation(real(z_v_ii), imag(z_v_ii));

            % Co-ordinates of each edge
            edge_coords  = edges(DT);

            % Weights of each edge (i.e. the distance between vortices)
            edge_weights = abs(z_v_ii(edge_coords(:,1)) - z_v_ii(edge_coords(:,2)));

            % Find minimum spanning tree
            [~, ST] = kruskal(edge_coords, edge_weights);
            
            % Plot each line (between each set of two vortices)
            for kk = 1:length(ST)
               plot([real(z_v_ii(ST(kk,1))), real(z_v_ii(ST(kk,2)))], [imag(z_v_ii(ST(kk,1))), imag(z_v_ii(ST(kk,2)))], '-', 'color', col_pos, 'linewidth', lw);
            end
        % If there are only 2 vortices, or if n vortices had the same x/y position, join them together
        else %if z_v == 2 
            for i = 1:length(z_v_ii)-1
                plot([real(z_v_ii(i)), real(z_v_ii(i+1))], [imag(z_v_ii(i)), imag(z_v_ii(i+1))], '-', 'color', col_pos, 'linewidth', lw)
            end
        end
        
    end
    
    % Plot lines between vortices in negative clusters
    for ii = 1:max_cluster_neg
        
        % List of vortices in this cluster
        z_v_ii = (z_v(c_v == -ii)).';
        
        % If there are >2 vortices and they don't have the same x/y positions, calculate minimum spanning tree
        % *This is a problem if the vortices exist on a numerical grid - the edge coords can't be defined
        if length(z_v_ii) > 2 && range(real(z_v_ii)) > 0 && range(imag(z_v_ii)) > 0
            
            % Calculate Delaunay triangulation between vortices in this cluster
            DT = delaunayTriangulation(real(z_v_ii), imag(z_v_ii));

            % Co-ordinates of each edge
            edge_coords  = edges(DT);

            % Weights of each edge (i.e. the distance between vortices)
            edge_weights = abs(z_v_ii(edge_coords(:,1)) - z_v_ii(edge_coords(:,2)));

            % Find minimum spanning tree
            [~, ST] = kruskal(edge_coords, edge_weights);
            
            % Plot each line (between each set of two vortices)
            for kk = 1:length(ST)
               plot([real(z_v_ii(ST(kk,1))), real(z_v_ii(ST(kk,2)))], [imag(z_v_ii(ST(kk,1))), imag(z_v_ii(ST(kk,2)))], '-', 'color', col_neg, 'linewidth', lw)
            end
        % If there are only 2 vortices, or if n vortices had the same x/y position, join them together
        else %if z_v == 2
            for i = 1:length(z_v_ii)-1
                plot([real(z_v_ii(i)), real(z_v_ii(i+1))], [imag(z_v_ii(i)), imag(z_v_ii(i+1))], '-', 'color', col_neg, 'linewidth', lw)
            end
        end
        
    end
    
    %% Plot vortices last
    plot(real(z_pos),imag(z_pos),'o','MarkerEdgeColor','k','MarkerFaceColor',col_pos,'MarkerSize',6,'linewidth',1)
    plot(real(z_neg),imag(z_neg),'o','MarkerEdgeColor','k','MarkerFaceColor',col_neg,'MarkerSize',6,'linewidth',1)

end

%% The remaining functions are from
%% https://au.mathworks.com/matlabcentral/fileexchange/41963-kruskal-s-algorithm
%% Reproduced under the following license conditions:

%% Copyright (c) 2014, Georgios Papachristoudis
%% All rights reserved.
%%
%% Redistribution and use in source and binary forms, with or without
%% modification, are permitted provided that the following conditions are met:
%%
%% * Redistributions of source code must retain the above copyright notice, this
%%   list of conditions and the following disclaimer.
%%
%% * Redistributions in binary form must reproduce the above copyright notice,
%%   this list of conditions and the following disclaimer in the documentation
%%   and/or other materials provided with the distribution
%% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
%% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Kruskal algorithm - needed for plotting (original author in description)
function [w_st, ST, X_st] = kruskal(X, w)
%
% This function finds the minimum spanning tree of the graph where each
% edge has a specified weight using the Kruskal's algorithm.
% 
% Assumptions
% -----------
%     N:  1x1  scalar      -  Number of nodes (vertices) of the graph
%    Ne:  1x1  scalar      -  Number of edges of the graph
%   Nst:  1x1  scalar      -  Number of edges of the minimum spanning tree
% 
% We further assume that the graph is labeled consecutively. That is, if
% there are N nodes, then nodes will be labeled from 1 to N.
%
% INPUT
% 
%     X:  NxN logical      -  Adjacency matrix
%             matrix          If X(i,j)=1, this means there is directed edge
%                             starting from node i and ending in node j.
%                             Each element takes values 0 or 1.
%                             If X symmetric, graph is undirected.
% 
%  or     Nex2 double      -  Neighbors' matrix
%              matrix         Each row represents an edge.
%                             Column 1 indicates the source node, while
%                             column 2 the target node.
% 
%     w:  NxN double       -  Weight matrix in adjacency form
%             matrix          If X symmetric (undirected graph), w has to
%                             be symmetric.
% 
%  or     Nex1 double      -  Weight matrix in neighbors' form
%              matrix         Each element represents the weight of that
%                             edge.
% 
% 
% OUTPUT
% 
%  w_st:    1x1 scalar     -  Total weight of minimum spanning tree
%    ST:  Nstx2 double     -  Neighbors' matrix of minimum spanning tree
%               matrix
%  X_st:  NstxNst logical  -  Adjacency matrix of minimum spanning tree
%                 matrix      If X_st symmetric, tree is undirected.
% 
% EXAMPLES
%
% Undirected graph
% ----------------
% Assume the undirected graph with adjacency matrix X and weights w:
%
%         1   
%       /   \
%      2     3
%     / \
%    4 - 5
% 
% X = [0 1 1 0 0;
%      1 0 0 1 1;
%      1 0 0 0 0;
%      0 1 0 0 1;
%      0 1 0 1 0];
%  
% w = [0 1 2 0 0;
%      1 0 0 2 1;
%      2 0 0 0 0;
%      0 2 0 0 3;
%      0 1 0 3 0];
% 
% [w_st, ST, X_st] = kruskal(X, w);
% The above function gives us the minimum spanning tree.
% 
% 
% Directed graph
% ----------------
% Assume the directed graph with adjacency matrix X and weights w:
%
%           1
%        / ^ \
%       / /   \
%      v       v
%       2 ---> 3
% 
% X = [0 1 1
%      1 0 1
%      0 0 0];
%  
% w = [0 1 4;
%      2 0 1;
%      0 0 0];
% 
% [w_st, ST, X_st] = kruskal(X, w);
% The above function gives us the minimum directed spanning tree.
% 
% 
% Author: Georgios Papachristoudis
% Copyright 2013 Georgios Papachristoudis
% Date: 2013/05/26 12:25:18

    isUndirGraph = 1;
    
    % Convert logical adjacent matrix to neighbors' matrix    
    if size(X,1)==size(X,2) && sum(X(:)==0)+sum(X(:)==1)==numel(X)        
        if any(any(X-X'))
            isUndirGraph = 0;
        end
        ne = cnvrtX2ne(X,isUndirGraph);
    else
        if size(unique(sort(X,2),'rows'),1)~=size(X,1)
            isUndirGraph = 0;
        end
        ne = X;
    end
    
    % Convert weight matrix from adjacent to neighbors' form
    if numel(w)~=length(w)
        if isUndirGraph && any(any(w-w'))
            error('If it is an undirected graph, weight matrix has to be symmetric.');
        end
        w = cnvrtw2ne(w,ne);
    end
    
    N    = max(ne(:));   % number of vertices
    Ne   = size(ne,1);   % number of edges    
    lidx = zeros(Ne,1);  % logical edge index; 1 for the edges that will be
                         % in the minimum spanning tree                         
    % Sort edges w.r.t. weight
    [w,idx] = sort(w);
    ne      = ne(idx,:);
    
    % Initialize: assign each node to itself
    [repr, rnk] = makeset(N);
    
    % Run Kruskal's algorithm
    for k = 1:Ne
        i = ne(k,1);
        j = ne(k,2);
        if fnd(i,repr) ~= fnd(j,repr)
            lidx(k) = 1;
            [repr, rnk] = union(i, j, repr, rnk);
        end
    end
    
    % Form the minimum spanning tree
    treeidx = find(lidx);
    ST      = ne(treeidx,:);
    
    % Generate adjacency matrix of the minimum spanning tree
    X_st = zeros(N);
    for k = 1:size(ST,1)
        X_st(ST(k,1),ST(k,2)) = 1;
        if isUndirGraph,  X_st(ST(k,2),ST(k,1)) = 1;  end
    end
    
    % Evaluate the total weight of the minimum spanning tree
    w_st = sum(w(treeidx));
end

function ne = cnvrtX2ne(X, isUndirGraph)
    if isUndirGraph
        ne = zeros(sum(sum(X.*triu(ones(size(X))))),2);
    else
        ne = zeros(sum(X(:)),2);
    end
    cnt = 1;
    for i = 1:size(X,1)
        v       = find(X(i,:));        
        if isUndirGraph
            v(v<=i) = [];
        end
        u       = repmat(i, size(v));
        edges   = [u; v]';
        ne(cnt:cnt+size(edges,1)-1,:) = edges;
        cnt = cnt + size(edges,1);
    end
end

function w = cnvrtw2ne(w,ne)
    tmp = zeros(size(ne,1),1);
    cnt = 1;
    for k = 1:size(ne,1)
        tmp(cnt) = w(ne(k,1),ne(k,2));
        cnt = cnt + 1;
    end
    w = tmp;
end

function [repr, rnk] = makeset(N)
    repr = (1:N);
    rnk  = zeros(1,N);
end

function o = fnd(i,repr)
    while i ~= repr(i)
        i = repr(i);
    end
    o = i;
end

function [repr, rnk] = union(i, j, repr, rnk)
    r_i = fnd(i,repr);
    r_j = fnd(j,repr);
    if rnk(r_i) > rnk(r_j)
        repr(r_j) = r_i;
    else
        repr(r_i) = r_j;
        if rnk(r_i) == rnk(r_j)
            rnk(r_j) = rnk(r_j) + 1;
        end
    end
end