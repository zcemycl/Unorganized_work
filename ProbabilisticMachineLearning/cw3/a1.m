% find the highest probability of words in A
% clear all;clc;
% find the index 
load kos_doc_data.mat

total = sum(A(:,3));
unique_id = unique(A(:,2),'first');
dim_uid = size(unique_id);

countid = []; % independent of the document id
A3 = A(:,3);
for i = 1:dim_uid(1)
    rows = find(A(:,2)==unique_id(i));
    countid(i) = sum(A3(rows));
end

probid = countid/total;
figure(1)
plot(unique_id, probid);

probidsort = sort(probid,'descend');
probidtop20 = probidsort(1:20);
top20_wid = [];
for j = 1:20
    index = find(probid == probidtop20(j));
    top20_wid(j) = unique_id(index);
end

figure(2)
barh(probidtop20)
set(gca,'YTickLabel',V(top20_wid),'YTick',1:20,'FontSize',6)


