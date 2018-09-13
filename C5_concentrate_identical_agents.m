function [] = C5_concentrate_identical_agents(nn, c)
%
% To improve convergence of our soft assign algorithm, where there are many
% identical agents this function concentrates identical agents together in
% order to improve integrality.
%
% E.g. two agents each 0.5 on two locations gets transform so each agent
% fully occupies a single location.
%

%% Preliminaries & Parameters

% globals
global probs
global agent_info
global interaction_reward

% list of agent types
agent_types = unique(interaction_reward(:,c.IR.to));
kk = length(agent_types);

%% Concentrate

% loop through agent types
for ii = 1:kk
    % select agents
    agent_nums = find(agent_info(:,c.AI.type) == agent_types(ii));
    % sum probs for these agents
    probs_by_type = sum(probs(agent_nums,:),1);
    % clear current probs
    probs(agent_nums,:) = 0;
    
    % mid check
    assert(abs(sum(probs_by_type) - length(agent_nums)) <= 1E-8,'agents and areas occupied do not coincide')
    
    % current agent and location
    c_agent = 1;
    c_locat = 1;
    % fill in concentrated probabilities
    while c_agent <= length(agent_nums) && c_locat <= nn
        % current prob awaiting assignment for this agent and location
        ca_prob = 1 - sum(probs(agent_nums(c_agent),:));
        cl_prob = probs_by_type(c_locat) - sum(probs(agent_nums,c_locat));
        % assign prob to current agent/location pair
        probs(agent_nums(c_agent),c_locat) = min(cl_prob, ca_prob);
        % increment location once full
        if cl_prob == min(cl_prob, ca_prob)
            c_locat = c_locat + 1;
        end
        % increment agent once full
        if ca_prob == min(cl_prob, ca_prob)
            c_agent = c_agent + 1;
        end
    end
end






end