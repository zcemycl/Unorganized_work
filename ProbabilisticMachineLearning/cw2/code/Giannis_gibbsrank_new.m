load tennis_data

randn('seed',27); % set the pseudo-random number generator seed

M = size(W,1);            % number of players
N = size(G,1);            % number of games in 2011 season 

pv = 0.5*ones(M,1);           % prior skill variance 

w = zeros(M,1);               % set skills to prior mean

for i = 1:500

  % First, sample performance differences given the skills and outcomes
  
  t = nan(N,1); % contains a t_g variable for each game
  
  for g = 1:N   % loop over games
    s = w(G(g,1))-w(G(g,2));  % difference in skills
    t(g) = randn()+s;         % performace difference sample
    while t(g) < 0  % rejection sampling: only positive perf diffs accepted
      t(g) = randn()+s; % if rejected, sample again
    end
  end 
 
  
  % Second, jointly sample skills given the performance differences
  
  m = zeros(M,1);  % container for the mean of the conditional 
                 % skill distribution given the t_g samples
  for p = 1:M
      for g=1:N
          if(G(g,1)==p)
            m(p)=m(p)+t(g) ; % (***TO DO***) complete this line
          end
          if(G(g,2)==p)
              m(p)=m(p)-t(g) ; 
          end
      end
  end
  
  iS = zeros(M,M); % container for the sum of precision matrices contributed
                   % by all the games (likelihood terms)
for p=1:M
    for p2=1:M
        if(p==p2) %%ii case
            for g=1:N
                if(G(g,1)==p)
                iS(p,p2)=iS(p,p2)+1;
                elseif(G(g,2)==p)
                iS(p,p2)=iS(p,p2)+1;
                end
            end
        elseif(p~=p2)%%ij case
              for g=1:N
                if(G(g,1)==p && G(g,2)==p2)
                iS(p,p2)=iS(p,p2)-1;
                elseif(G(g,2)==p && G(g,1)==p2)
                iS(p,p2)=iS(p,p2)-1;
                end
              end
        end
    end
end
    

  iSS = diag(1./pv) + iS; % posterior precision matrix
  % prepare to sample from a multivariate Gaussian
  % Note: inv(M)*z = R\(R'\z) where R = chol(M);
  iR = chol(iSS);  % Cholesky decomposition of the posterior precision matrix
  mu = iR\(iR'\m); % equivalent to inv(iSS)*m but more efficient
    
  % sample from N(mu, inv(iSS))
  w = mu + iR\randn(M,1);
  
  Ms = w;
  cw2;
  drawnow;
  
end
