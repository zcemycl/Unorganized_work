clear all; 
load tennis_data

randn('seed',27); % set the pseudo-random number generator seed

M = size(W,1);            % number of players 107
N = size(G,1);            % number of games in 2011 season 1801

pv = 0.5*ones(M,1);           % prior skill variance 

w = zeros(M,1);               % set skills to prior mean
figure(1);
w_all = [];
for i = 1:1100

  % First, sample performance differences given the skills and outcomes
  
  t = nan(N,1); % contains a t_g variable for each game
  for g = 1:N   % loop over games
    s = w(G(g,1))-w(G(g,2));  % difference in skills
    t(g) = randn()+s;         % p10erformace difference sample
    while t(g) < 0  % rejection sampling: only positive perf diffs accepted
      t(g) = randn()+s; % if rejected, sample again
    end
  end 

  % Second, jointly sample skills given the performance differences
  
  m = nan(M,1);  % container for the mean of the conditional 
                 % skill distribution given the t_g samples
  for p = 1:M  
    countw = 1;   
    countl = 1;
    % --------------------------------------------------------- %
    % General for all players
    wlist = []; llist = [];
    indices = find(G==p);
    for k = 1:size(indices)
        index = indices(k);
        if index <= 1801
            wlist(countw) = t(index);
            countw = countw + 1;
        else
            index = index - 1801;
            llist(countl) = -t(index);
            countl = countl + 1;
        end
    end
    m(p) = sum(wlist)+sum(llist);
    % --------------------------------------------------------- %     
  end
  
  iS = zeros(M,M); % container for the sum of precision matrices contributed
                       % by all the games (likelihood terms)
  for g = 1:N
      
          
      % --------------------------------------------------------- %
      % General for all players
      players_involved = G(g,:);
      index_win = players_involved(1);
      index_lose= players_involved(2);

      % diagonals
      iS(index_win,index_win) = iS(index_win,index_win)+1;
      iS(index_lose,index_lose) = iS(index_lose,index_lose)+1;
      % non-diagonals
      iS(index_win,index_lose) = iS(index_win,index_lose)-1;
      iS(index_lose,index_win) = iS(index_lose,index_win)-1;

      % --------------------------------------------------------- %
  end

  iSS = diag(1./pv) + iS; % posterior precision matrix
  % prepare to sample from a multivariate Gaussian
  % Note: inv(M)*z = R\(R'\z) where R = chol(M);
  iR = chol(iSS);  % Cholesky decomposition of the posterior precision matrix
  mu = iR\(iR'\m); % equivalent to inv(iSS)*m but more efficient
    
  % sample from N(mu, inv(iSS))
  w = mu + iR\randn(M,1);
  
  % --------------------------------------- %
  % Optimize better ranking difference
  w_copy = w;
  num = randi([1,M],1,2);
  num1 = num(1); num2 = num(2);
  while num2 == num1
      num2 = randi([1,M],1,1);
  end
  w_copy(num1) = w_copy(num2);
  % w = w_copy;
  w_all(1:M,i) = w_copy;
  w = sum(w_all,2)/i;
  % --------------------------------------- %
  Ms = w;
  cw2;
  hold on;
  title(['Iteration:', num2str(i)])
  drawnow;
  hold off;
end


















