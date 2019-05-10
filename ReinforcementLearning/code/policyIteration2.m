function [v, pi] = policyIteration2(model, maxit)

% initialize the value function
v = zeros(model.stateCount, 1);
pi = ones(model.stateCount, 1);
for i = 1:maxit
    v_ = zeros(model.stateCount, 1);
    
    for j = 1:maxit
        for s = 1:model.stateCount
            tmpaweight = model.P(s,:,pi(s))*(model.R(s,pi(s))+model.gamma*v(:));
            v_(s) = sum(tmpaweight);       
        end
        
        if max(v_-v)<0.001 %norm or max
            break
        else
            v = v_;
        end
    end
    
    % for policy 
    for s = 1:model.stateCount
        tmpa = zeros(4,1);
        for a = 1:4
            tmpa(a) = model.P(s,:,a)*(model.R(s,a)+model.gamma*v(:));
        end

        [V,choice] = max(tmpa);
        tmpV =  model.P(s,:,pi(s))*(model.R(s,pi(s))+model.gamma*v(:));
        if V >= tmpV
            pi(s) = choice;
        end
    end

end
end
