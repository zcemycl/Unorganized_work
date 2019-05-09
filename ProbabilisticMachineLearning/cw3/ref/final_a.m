% Final a
close all; clear all; clc;
load kos_doc_data.mat

A1 = A(:,1);
A2 = A(:,2);
A3 = A(:,3);
tot = sum(A3);

% extract word id
unique_id = unique(A(:,2));
dim_uid = size(unique_id);

% count of words with the same 
% id, independent of document
countid = [];
for i = 1:dim_uid(1)
    countid(i) = sum(A3(find(A2 == unique_id(i))));
end

% probability of word id
probid = countid/tot;

% ------------------------------ %
% Only for part a histogram
% ------------------------------ %
probidsort = sort(probid, 'descend');
probid20 = probidsort(1:20);
top20w_id = [];
for j = 1:20
    top20w_id(j) = unique_id(find(probid == probid20(j)));
end

figure(1)
barh(probid20)
set(gca,'YTickLabel',V(top20w_id),'YTick',1:20,'FontSize',6)
xlabel('Probability')

% Perplexity 
p = 1/6906;


