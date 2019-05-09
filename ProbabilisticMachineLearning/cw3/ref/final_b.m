% final b
% load final_a.mat
final_a;
alpha = 0.1;
B1 = B(:,1);
B2 = B(:,2);
B3 = B(:,3);

probDid = (alpha+countid)./(alpha*max(unique_id) + tot);

figure(2)
subplot(1,2,1)
plot(unique_id, probid)
subplot(1,2,2)
plot(unique_id, probDid)

% ------------------------------- %
% Only for part b
% ------------------------------- %
figure(3)
% Top 3 common words
probDidsort = sort(probDid, 'descend');
probDid3 = probDidsort(1:3);
dimD = size(probDidsort);
top3w_id = [];
for k = 1:3
    top3w_id(k) = unique_id(find(probDid == probDid3(k)));
end

% Top 3 rare words
probDidlast3 = probDidsort(dimD(2));
last3w_id = [];
inter = find(probDid == probDidlast3);
last3w_id = inter(1:3);
probDidlast3 = probDidlast3*ones(1,3);

% Non-existing words x3
C = setdiff(B2,A2);
probnon = (alpha/(alpha*max(unique_id)+tot))*ones(1,3);

proball = [probDid3, probDidlast3, probnon];
index = [top3w_id, last3w_id, C(1:3)'];
barh(proball)
set(gca,'YTickLabel',V(index),'YTick',1:9,'FontSize',5)
xlabel('Probability')



