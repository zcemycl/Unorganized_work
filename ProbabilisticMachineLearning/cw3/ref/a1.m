load kos_doc_data.mat

W = max([A(:,2); B(:,2)]);  % number of unique words

total = sum(A(:,3));
freq = A(:,3)/total;

for i=1:W
   freq_unique(i) = sum(A(A(:,2)==i,3))/total; 
end
[sort_freq, ind] = sort(freq_unique,'descend');
sort_wordID = A(ind,2);
twentyID = sort_wordID(1:20);

for j = 1:W
   l(j) = log(freq_unique(j)).* sum(A(A(:,2)==j,3));
end
l_total = nansum(l);

np=20;
barh(sort_freq(np:-1:1))
set(gca,'YTickLabel',V(ind(np:-1:1)),'YTick',1:np,'FontSize',6)
%set(gca,'YTickLabel',V(s),'Ytick',1:np)
%axis([0 1 0.5 np+0.5])
title('Frequency of 20 most probable words','FontSize',10);
xlabel('Frequency','FontSize',8)