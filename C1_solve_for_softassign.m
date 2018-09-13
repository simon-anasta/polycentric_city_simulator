function [assignment, alpha_list, costs, log_prob, log_obj, nn, c ] = C1_solve_for_softassign(agent_info_file, locat_info_file, linear_rewards_file, interaction_reward_file )
%
% Progressively tightens agents preference for the option that gives them
% the best outcome. Outer function for softassign heuristic algorithm.
%

%% Preliminaries & Parameters

% problem size
[nn, c] = C2_load_global_data(agent_info_file, locat_info_file, linear_rewards_file, interaction_reward_file );
% convergence accuracy
req_accur = 0.01;
% display iterations
disp_iters = 0;
plot_iters = 1; % warning! very slow

% alpha parameter values
alpha_start = 0.5;
alpha_scale = 1.05;
alpha_iters = 180;

% gamma parameter values
gamma = 0.5;

% global access
global probs
global indiv_reward

%% Checks

%% Initialize

% lists of alpha and values
alpha_list = alpha_start * alpha_scale.^(0:alpha_iters);
% alpha_list = sort([alpha_list,alpha_list]);

% initial probabilities
probs = ones(nn,nn) / nn;
% initial costs
costs = zeros(1,nn);

% log for probabilities
% log_prob = zeros(nn,nn,alpha_iters+1); % record all probs - small models only
log_prob = zeros(1,alpha_iters+1);     % record non-integrality - any size model
% log for objective value
log_obj = zeros(1,alpha_iters+1);

%% Iterative solver

for ii = 1:(alpha_iters+1)
    % current alpha value
    current_alpha = alpha_list(ii);
    
    % get rewards to individual agents
    C3_get_current_rewards(nn, c);
    % current assignment reward
    obj = sum(indiv_reward(:) .* probs(:));
    
    % resolve markets
    [costs] = C4_solve_for_probs(current_alpha, gamma, costs, nn);
    
    % concentrate agents
    C5_concentrate_identical_agents(nn, c);
    
    % plot location preferences (optional & slow)
    if plot_iters == 1 && mod(ii,10) == 0
        location_preference_plot(nn,c);
    end
    
    % variation from assignment
    nonInt = closeness_to_assign(nn);
    
    % store logs
%     log_prob(:,:,ii) = probs; % record all probs - small models only
    log_prob(ii) = nonInt;    % record non-integrality - any size model
    log_obj(ii) = obj;
    
    % inform user
    if disp_iters == 1
        fprintf('iteration %4d complete, alpha %5.2f, from objective %8.4f, non-integrality %6.5f\n',ii,current_alpha,obj,nonInt);
    end
    
    % exit condition
    if nonInt < req_accur
        break
    end
end

%% Prep output

% clear extra storage
alpha_list(ii+1:end) = [];
log_obj(ii+1:end) = [];
% log_prob(:,:,ii+1:end) = []; % if recorded all probs - small models only
log_prob(ii+1:end) = [];    % if recorded non-integrality - any size model

% summary of assignment
[~,assignment] = max(probs,[],2);

end

%% Sub function: closeness to assignment
function [nonInt] = closeness_to_assign(nn)
%
% computes how close the prob matrix is to an assignment matrix.
%

% get global
global probs
% initialize
nonInt = 0;

% iterate through rows of probs
for ii = 1:nn
    % values to compare
    this_row = probs(ii,:);
    this_max = max(this_row);
    % find strongest assignment
    i_max = find(this_row == this_max,1);
    % build closest assignment vector
    assign_vec = zeros(1,nn);
    assign_vec(i_max) = 1;
    % compute and add on distance
    nonInt = nonInt + sqrt( sum( (this_row - assign_vec).^2 ));
end

end

%% sub function: location preference plot
function [] = location_preference_plot(nn,c)

% location size
locat_size = 1;

% access globals
global probs
global agent_info
global locat_info

% city limits
xmin = min(locat_info(:,c.LI.x_coord)) - locat_size / 2;
xmax = max(locat_info(:,c.LI.x_coord)) + locat_size / 2;
ymin = min(locat_info(:,c.LI.y_coord)) - locat_size / 2;
ymax = max(locat_info(:,c.LI.y_coord)) + locat_size / 2;

% values to plot
values_to_plot = sum(probs(agent_info(:,c.AI.type) == 20,:));
values_to_plot = max(values_to_plot,0);
values_to_plot = min(values_to_plot,1);

% setup figure
AX = figure;
set(AX,'Visible','off')
axis([xmin, xmax, ymin, ymax])
axis equal; axis off
% iterate through agents
for this_locat = 1:nn
    % location coords
    this_locat_x = locat_info(this_locat,c.LI.x_coord);
    this_locat_y = locat_info(this_locat,c.LI.y_coord);
    % cost colour scale
    RGB = [0, values_to_plot(this_locat), values_to_plot(this_locat)];
    % plotting position
    position = [this_locat_x - locat_size/2, this_locat_y - locat_size/2, locat_size, locat_size];
    % plot
    rectangle('Position',position,'Curvature',[1,1],'FaceColor',RGB);
end
% update display
title('Location of agents')
set(AX,'Visible','on')

end