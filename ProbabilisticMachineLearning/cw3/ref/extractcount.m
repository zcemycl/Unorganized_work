function [index,count] = extractcount(data,index)
% extract document 2001
target = data(find(data(:,1)==index),:);

index = target(:,2);
count = target(:,3);


end