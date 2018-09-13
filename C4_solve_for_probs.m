function [costs] = C4_solve_for_probs(alpha, gamma, costs, nn)
%
% Applies Economic varient of Sinkhorn's iterative method to produce a
% doubly stochastic matrix for assignment of agents to locations.
%
% Economic varient chosen for superior performance over traditional.
%

%% Preliminaries & Parameters

% globals
global probs
global indiv_reward

% convergence accuracy
req_accur = 1E-3;
% inform user
disp_iters = 1;
% maximum iterations
max_iters = 100;

%% Checks

% checks for dimension
assert(length(alpha) == 1,'alpha must be a scalar')
% check for numeric
assert(isnumeric(alpha),'alpha must be numeric')
assert(~any(isnan(indiv_reward(:))),'rewards must all be numeric')

%% Self-consistency reward

% average reward
% avg_reward = mean(indiv_reward,2);
% reward self-consistency
indiv_reward = indiv_reward + gamma*probs;
% indiv_reward = indiv_reward + gamma*(avg_reward * ones(1,nn)).*probs;

%% Normalize by costs

for ii = 1:max_iters
    % current probs
    % scale by alpha
    probs = alpha * (indiv_reward - ones(nn,1) * costs);
    % resolve w.r.t. the largest reward
    % (otherwise vulnerable to overflow and underflow errors)
    probs = probs - ( max(probs,[],2) * ones(1,nn) ) + 1;
    % exponentiate
    probs = exp(probs);
    
    % normalize by row
    probs = probs ./ ( sum(probs,2) * ones(1,nn) );
    % current location sums
    sLocations = sum(probs,1);
    % error
    this_error = max(abs(sLocations - 1));
    
    % check exit condition
    if all(this_error < req_accur)
        break
    end
    
    % desired decrease in cost
%     dCosts = log(sLocations) / alpha;
    dCosts = sign(sLocations-1) .* log(1+abs(sLocations-1)) / alpha;
    % take step
    costs = costs + dCosts;
end

%% Output

% warn user if exit due to maximum iteration reached
if ii == max_iters
    fprintf('max interations reached by Sinkhorn, alpha %4.3f, worst error %9.8f\n',alpha,max(abs(sLocations - 1)));
% inform user of total number of iterations
elseif disp_iters == 1
    fprintf('Sinkhorn converged after %5d iterations\n',ii);
end

end