% clear all; 
load tennis_data

randn('seed',100); % set the pseudo-random number generator seed

M = size(W,1);            % number of players 107
N = size(G,1);            % number of games in 2011 season 1801

pv = 0.1*ones(M,1);           % prior skill variance 

w = 0.5*ones(M,1);               % set skills to prior mean
figure(1);
w_all = [];
w_sum_all = [];
var_all = [];
for i = 1:11000

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
  
  % Sampling all the thing and store them in the array.
  w_all(:,i) = w;
  w_sum_all(:,i) = sum(w_all,2)/i;
  % Variance
  var_all(:,i) = sum((w_all-w_sum_all).^2,2)/i;
  
  
%   figure(1)
%   Ms = w;
%   cw2;
%   hold on;
%   title(['Iteration:', num2str(i)])
%   drawnow;
%   hold off;
%   
%   figure(2)
%   Ms = w_sum_all(:,i);
%   cw2;
%   hold on;
%   title(['Iteration:', num2str(i)])
%   drawnow;
%   hold off;
  
end


% plot one player
% means and samples of the player
figure(3)
plot(w_all(16,:));
hold on; plot(w_sum_all(16,:));

% variance of the player skill
figure(4)
plot(var_all(16,:));

% Take every 5th sample
g_sample = [];
index = [];
for i = 1:220
    index(i) = 5*i;
    g_sample(:,i) = w_all(:,5*i);
end

% MEAN 
MEAN = sum(g_sample,2)/220;
% VARIANCE
VAR = [];
VAR = sum((g_sample - MEAN).^2,2)/220;

