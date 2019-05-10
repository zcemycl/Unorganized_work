function [v, pi] = valueIteration(model, maxit)

% initialize the value function
v = zeros(model.stateCount, 1);

for i = 1:maxit
    v_ = zeros(model.stateCount, 1);
    pi_= ones(model.stateCount, 1);

    for s = 1:model.stateCount
        tmpa = zeros(4,1);
        for a = 1:4 
            tmpa(a) = model.P(s,:,a)*(model.R(s,a)+model.gamma*v(:));    
        end

        [v_(s),pi_(s)] = max(tmpa);
    end
    
    if max(v_-v)<0.01
        break;
    else 
        v = v_;
        pi = pi_;
    end
    
end
end
