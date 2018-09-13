function [nn, c] = C2_load_global_data(agent_info_file, locat_info_file, linear_rewards_file, interaction_reward_file )
%
% Loads data sources for city problem. Inputs are strings for files to
% load. Output is problem size, column headers and creation of global
% variables.
%
% Global variables used as a way to save memory. Passing matrices as
% agruments can create duplicate copies of the objects filling more memory
% that desired.
%

%% Preliminaries & Parameters

% global access
clear('agent_info', 'locat_info', 'linear_rewards')
global agent_info
global locat_info
global linear_rewards
global interaction_reward

%% Load data

% agent information list
agent_info = importdata(agent_info_file, ',', 1);
c.AI = column_names_to_structure(agent_info.textdata);
agent_info = agent_info.data;

% problem size
nn = size(agent_info,1);
% check problem fits in memory
mm = memory();
assert(sqrt(mm.MaxPossibleArrayBytes/24) >= nn, 'problem size is larger than RAM can handle')
% check no missing values
assert(~any(isnan(agent_info(:))), 'missing values in agent info')

% location information list
locat_info = importdata(locat_info_file, ',', 1);
c.LI = column_names_to_structure(locat_info.textdata);
locat_info = locat_info.data;
% check size match
assert(size(locat_info,1) == nn,'identical numbers of agents and locations required')
% check no missing values
assert(~any(isnan(locat_info(:))), 'missing values in location info')

% linear rewards
linear_rewards = importdata(linear_rewards_file, ',', 1);
% drop first column
linear_rewards = linear_rewards.data(:,2:end);
% add tie breaking
linear_rewards = linear_rewards + rand(nn,nn);
% check size match
assert(size(linear_rewards,1) == nn,'rows of linear rewards do not correspond to agent info')
assert(size(linear_rewards,2) == nn,'columns of linear rewards do not correspond to location info')
% check no missing values
assert(~any(isnan(linear_rewards(:))), 'missing values in linear rewards')

% interaction reward
interaction_reward = importdata(interaction_reward_file, ',', 1);
c.IR = column_names_to_structure(interaction_reward.textdata);
interaction_reward = interaction_reward.data;
% check no missing values
assert(~any(isnan(interaction_reward(:))), 'missing values in interaction reward')

end

%% Sub function
function [c] = column_names_to_structure(column_names)
%
% Loop through column names building a structure with the number of the
% column.
%

% keep only first row
column_names = column_names(1,:);

for ii = 1:length(column_names)
    eval(['c.',column_names{ii},' = ',num2str(ii),';'])
end

end