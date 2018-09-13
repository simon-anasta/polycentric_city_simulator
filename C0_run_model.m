function [] = C0_run_model()
%
% Runs the entire model, solving and producing output.
%

%% Parameters

% files
agent_info_file = 'agent_info.csv';
linear_rewards_file = 'linear_rewards.csv';
interaction_reward_file = 'interaction_reward.csv';
locat_info_file = 'locat_info.csv';

% location size
locat_size = 1;

%% Run model

% time
t1 = tic();
% run
[assignment, alpha_list, costs, log_prob, log_obj, nn, c ] = C1_solve_for_softassign(agent_info_file, locat_info_file, linear_rewards_file, interaction_reward_file );
% inform user
fprintf('Model size %d complete in %2.4f minutes\n',nn,toc(t1)/60)

%% Plot performance output

% time
t2 = tic();

figure;

% integrality
subplot(2,2,1);plot(log_prob)
xlabel('iteration'); ylabel('non-assignment')
title('Distance from integral assignment over time')

% integrality vs alpha
subplot(2,2,3);plot(alpha_list,log_prob)
xlabel('alpha value'); ylabel('non-integrality')
title('Distance from integral assignment w.r.t. alpha')

% objective value
subplot(2,2,2);plot(log_obj)
xlabel('iteration'); ylabel('objective value')
title('Total value of city over time')

% objective value vs. gamma
subplot(2,2,4);plot(alpha_list,log_obj)
xlabel('alpha value'); ylabel('objective value')
title('Total value of city w.r.t. alpha')

%% Map setup

% access global data
global agent_info
global locat_info

% city limits
xmin = min(locat_info(:,c.LI.x_coord)) - locat_size / 2;
xmax = max(locat_info(:,c.LI.x_coord)) + locat_size / 2;
ymin = min(locat_info(:,c.LI.y_coord)) - locat_size / 2;
ymax = max(locat_info(:,c.LI.y_coord)) + locat_size / 2;

%% Map of which agents are where

% setup figure
AX = figure;
set(AX,'Visible','off')
axis([xmin, xmax, ymin, ymax])
axis equal; axis off
% iterate through agents
for ii = 1:nn
    % location
    this_locat = assignment(ii);
    % location coords
    this_locat_x = locat_info(this_locat,c.LI.x_coord);
    this_locat_y = locat_info(this_locat,c.LI.y_coord);
    % agent colour
    RGB = [agent_info(ii,c.AI.red),agent_info(ii,c.AI.green),agent_info(ii,c.AI.blue)] / 255;
    % plotting position
    position = [this_locat_x - locat_size/2, this_locat_y - locat_size/2, locat_size, locat_size];
    % plot
    rectangle('Position',position,'Curvature',[1,1],'FaceColor',RGB);
end
% update display
title('Location of agents')
set(AX,'Visible','on')

%% Map of cost of location

% range of costs
min_cost = min(costs);
max_cost = max(costs);

% setup figure
AX = figure;
set(AX,'Visible','off')
axis([xmin, xmax, ymin, ymax])
axis equal; axis off
% iterate through locations
for this_locat = 1:nn
    % location coords
    this_locat_x = locat_info(this_locat,c.LI.x_coord);
    this_locat_y = locat_info(this_locat,c.LI.y_coord);
    % cost colour scale
    RGB = [0, (costs(this_locat) - min_cost) / (max_cost - min_cost), 0];
    % plotting position
    position = [this_locat_x - locat_size/2, this_locat_y - locat_size/2, locat_size, locat_size];
    % plot
    rectangle('Position',position,'Curvature',[1,1],'FaceColor',RGB);
end
% update display
title('Cost of locations')
set(AX,'Visible','on')

% inform user
fprintf('Plotting of model complete in %2.4f minutes\n',toc(t2)/60)

end
