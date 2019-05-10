% e1
% Plot
figure(1)
hold on
for i = 1:20
plot(topic_mixpro(i,:))
end
legend(num2str(linspace(1,20,20)'),'Location','bestoutside')
xlabel('iteration')
ylabel('mixing proportion')
xlim([1,100])

% Perplexity for documents in B


% Word entropy of topic 
% as a funciton of iterations
figure(2)
hold on 
for j=1:20
    plot(entropy(j,:));
end
legend(num2str(linspace(1,20,20)'),'Location','bestoutside')
xlabel('iteration')
ylabel('word entropy')
xlim([1,20])






