function C3_get_current_rewards(nn, c)
%
% Given each agents estimates of the probability of the location of all
% other agents, this function computes the expected individual reward to
% each agent of selecting each location.
%

%% Preliminaries & Parameters

% globals
global agent_info
global locat_info
global linear_rewards
global interaction_reward
global probs
global indiv_reward

% list of agent types
agent_types = unique(interaction_reward(:,c.IR.to));
kk = length(agent_types);

% number of distance cost intervals
ss = 3; % NOTE: CONSEQUENCES OF ss ARE HARD CODED DUE TO FUNCTIONAL FORM

%% Checks

% check for probability
assert(all(probs(:) <= 1),'all probs must be less than one')
assert(all(probs(:) >= 0),'all probs must be greater than zero')
assert(~any(isnan(probs(:))),'all probs must be numeric')

%% Initialize

%% Fill output

% location rewards
indiv_reward = linear_rewards;

% sum by agent types
probs_by_type = zeros(kk,nn);
% loop through agent types
for ii = 1:kk
    % select agents
    i_this_type = agent_info(:,c.AI.type) == agent_types(ii);
    % aggregate agents by location type
    probs_by_type(ii,:) = sum(probs(i_this_type,:),1);
end

% loop through location
for uu = 1:nn
    % current location
    this_x_coord = locat_info(uu,c.LI.x_coord);
    this_y_coord = locat_info(uu,c.LI.y_coord);
    % distance
    this_dist = sqrt( (this_x_coord - locat_info(:,c.LI.x_coord)).^2 + (this_y_coord - locat_info(:,c.LI.y_coord)).^2 );
    this_dist = this_dist';
    
    % loop through agent types
    for ij = 1:(kk^2)
        % agent types
        type_agent_from = interaction_reward(ij,c.IR.from);
        type_agent_to   = interaction_reward(ij,c.IR.to);
        
        % gather by distance classes
        reward_by_class = zeros(ss,1);
        % identify correct row
        i_row = agent_types == type_agent_from;
        
        % first distance class
        aa = (this_dist - interaction_reward(ij,c.IR.x1)) / (interaction_reward(ij,c.IR.x2) - interaction_reward(ij,c.IR.x1));
        reward = (1 - aa) * interaction_reward(ij,c.IR.y1) .* probs_by_type(i_row,:) + aa * interaction_reward(ij,c.IR.y2) .* probs_by_type(i_row,:);
        reward_by_class(1) = sum( reward(0 <= aa & aa < 1) );
        % second distance class
        aa = (this_dist - interaction_reward(ij,c.IR.x2)) / (interaction_reward(ij,c.IR.x3) - interaction_reward(ij,c.IR.x2));
        reward = (1 - aa) * interaction_reward(ij,c.IR.y2) .* probs_by_type(i_row,:) + aa * interaction_reward(ij,c.IR.y3) .* probs_by_type(i_row,:);
        reward_by_class(2) = sum( reward(0 <= aa & aa < 1) );
        
        % compute rewards
        this_reward = sum(reward_by_class);
        % assign rewards
        i_this_type = agent_info(:,c.AI.type) == type_agent_to;
        indiv_reward(i_this_type,uu) = indiv_reward(i_this_type,uu) + this_reward;
    end
    
    % add bonus reward
%     indiv_reward(:,uu) = indiv_reward(:,uu) + locat_info(:,c.AI.bonus);
end

% scale down by overall size
indiv_reward = indiv_reward / nn;

end