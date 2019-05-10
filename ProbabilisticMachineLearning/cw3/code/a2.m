% a2
% find the highest probability of words in B
clear all;clc;
% find the index 
load kos_doc_data.mat

total = sum(B(:,3));
unique_id = unique(B(:,2),'first');
dim_uid = size(unique_id);

countid = []; % independent of the document id
A3 = B(:,3);
for i = 1:dim_uid(1)
    rows = find(B(:,2)==unique_id(i));
    countid(i) = sum(A3(rows));
end

probid = countid/total;
figure(3)
plot(unique_id, probid);

probidsort = sort(probid,'descend');
probidtop20 = probidsort(1:20);
top20_wid = [];
for j = 1:20
    index = find(probid == probidtop20(j));
    top20_wid(j) = unique_id(index);
end

figure(4)
barh(probidtop20)
set(gca,'YTickLabel',V(top20_wid),'YTick',1:20,'FontSize',6)

% Set difference between A and B
C = setdiff(A(:,2),B(:,2));
V(C)
D = setdiff(B(:,2),A(:,2));
V(D)
