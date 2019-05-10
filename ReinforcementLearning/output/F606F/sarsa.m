function [v, pi, all] = sarsa(model, maxit, maxeps, atype,alpha,etype,ep)
% initialize the value function
suceps = 100;
Q = zeros(model.stateCount, 4);
all = zeros(maxeps/suceps,1);
accre = 0;
for i = 1:maxeps
    if etype == "iter"
        ep = 1;
    end
    % every time we reset the episode, start at the given startState
    s = model.startState;
    probability = rand;
    if probability < 1-ep
        [ga,a] = max(Q(s,:));
    else
        a = randi(4);
    end
    
    for j = 1:maxit 
        if atype == "iter"
            alpha = 1/j;
        end
        if etype == "iter" 
            ep = 1/j;
        end
        accre = accre + model.gamma*model.R(s,a);
        
        p = 0;
        r = rand;

        for s_ = 1:model.stateCount
            p = p + model.P(s, s_, a);
            if r <= p
                break;
            end
        end

        % s_ should now be the next sampled state.
        % IMPLEMENT THE UPDATE RULE FOR Q HERE.
        probability = rand;
        if probability < 1-ep
            [ga_,a_] = max(Q(s_,:));
        else
            a_ = randi(4);
        end

        delta = model.R(s,a)+model.gamma*Q(s_,a_)-Q(s,a);
        Q(s,a) = Q(s,a)+alpha*delta;
        % SHOULD WE BREAK OUT OF THE LOOP?
        s = s_;
        a = a_;
                        
        
        if s == model.goalState
            break
        end

    end

    all(ceil(i/suceps)) = all(ceil(i/suceps))+accre;
    accre = 0;

end

% REPLACE THESE
[v,pi] = max(Q,[],2);


end


