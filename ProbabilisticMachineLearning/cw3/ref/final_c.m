% final_c
final_b;
tot_B = sum(B3);


% extract document 2001
B_2001 = B(find(B(:,1)==2001),:);
tot_B2001 = sum(B_2001(:,3));
B2_2001 = B_2001(:,2);
B3_2001 = B_2001(:,3);

% check whether document 2001 has the words that 
% training set A does not have
D = intersect(B2_2001, C);

% extract probabilities 
probw_2001 = [];
dim_d2001 = size(B2_2001);
for i = 1:dim_d2001(1)
    probw_2001(i) = probDid(find(unique_id==B2_2001(i)));
end

% --------------------------------------------------- %
% Log probabilities
Lc_without = sum(B3_2001.*log(probw_2001'));
Lc_with = logfact(tot_B2001)-sum(logfact(B3_2001))+Lc_without;
% --------------------------------------------------- %
% Per-word perplexity for the document ID 2001
PwP2001wo = exp(-Lc_without/tot_B2001);
PwP2001wi = exp(-Lc_with/tot_B2001);

% =================================================== %
uni_doc_id = unique(B(:,1));
% =================================================== %
% General
% group the words with the same index
% independent of the document
unique_idB = unique(B2);
probB = [];
countB = [];
dimB = size(unique_idB);
for z = 1:dimB(1)
    rows = find(B2==unique_idB(z));
    countB(z) = sum(B3(rows));
    row = find(unique_id==unique_idB(z));
    if isempty(row)
        probB(z) = max(probnon);
    else
        probB(z) = probDid(row);
    end
end

% Log probabilities
Lc_withoutB = sum(countB.*log(probB));
Lc_withB = Lc_withoutB+logfact(tot_B)-sum(logfact(countB));
% Per-word perplexity
PwPBwo = exp(-Lc_withoutB/tot_B);
PwPBwi = exp(-Lc_withB/tot_B);

% Find the no. of words that does not exist
% in A each document has
no_non = [];
countloop = 1;
for c = min(B1):max(B1)
    settmp = B(find(B(:,1)==c),:);
    settmp2 = settmp(:,2);
    dimension = size(intersect(settmp2, C));
    no_non(countloop) = dimension(1);
    countloop = countloop+1;
end



figure(4)
plot(linspace(2001,3430,1430),no_non);
ylabel('Number of non-existing words in A')
xlabel('Document index')
xlim([2001 3430])

figure(5)
subplot(1,2,1)
[index2001,count2001] = extractcount(B,2001);
stem(index2001,count2001)
xlim([1 6906])
title('Document 2001')
xlabel('Word index')
ylabel('Number of occurence')

subplot(1,2,2)
[index2100,count2100] = extractcount(B,2100);
stem(index2100,count2100)
xlim([1 6906])
title('Document 2100')
xlabel('Word index')
ylabel('Number of occurence')



