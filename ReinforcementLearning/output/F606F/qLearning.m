function [v, pi, all] = qLearning(model, maxit, maxeps, atype, alpha, etype, ep)
% initialize the value function
suceps = 100;
Q = zeros(model.stateCount, 4);
all = zeros(maxeps/suceps,1);
accre = 0;
for i = 1:maxeps
    % every time we reset the episode, start at the given startState
    s = model.startState;

    for j = 1:maxit
        if atype == "iter"
            alpha = 1/j;
        end
        if etype == "iter" 
            ep = 1/j;
        end
        probability = rand;
        if probability < 1-ep
            [ga,a] = max(Q(s,:));
        else
            a = randi(4);
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
        delta = model.R(s,a)+model.gamma*max(Q(s_,:))-Q(s,a);
        Q(s,a) = Q(s,a)+alpha*delta;
        % SHOULD WE BREAK OUT OF THE LOOP?
        s = s_;
        
        if s == model.goalState
            break
        end

    end
    
%     if rem(i,suceps) == 0
%         all(i/suceps) = accre/suceps;
%         accre = 0;
%     end
    all(ceil(i/suceps)) = all(ceil(i/suceps))+accre;
    accre = 0;
end

% REPLACE THESE
[v,pi] = max(Q,[],2);

end

