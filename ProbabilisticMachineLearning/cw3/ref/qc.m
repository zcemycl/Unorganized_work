% c 
load kos_doc_data.mat
 
 W = max([A(:,2); B(:,2)]);  % number of unique words
 for i=1:W
   count_unique(i) = sum(A(A(:,2)==i,3)); 
end
 
 wordID = B(B(:,1)==2001, 2);
 wordCount = B(B(:,1)==2001, 3);
 alpha = 0.1;
 
 pi = ((alpha+count_unique(wordID))./(alpha*W+sum(A(:,3))))';
 ll = sum(wordCount.*log(pi));
 perplexity = exp(-(ll/sum(wordCount)));
 
 %%%%%%%%%%%%%% All documents in B %%%%%%%%%%%%%%%%%%%%
 allwordID = unique(B(:, 2));
 allwordCount = zeros(1,max(allwordID));
 for j = 1:max(allwordID)
     allwordCount(j) = sum(B(B(:,2)==j,3));
 end
 allwordCount = allwordCount';
 all_pi = ((alpha+count_unique(allwordID))./(alpha*W+sum(A(:,3))))';
 all_ll = (allwordCount(allwordID))'*log(all_pi);
 all_perplexity = exp(-(all_ll/sum(allwordCount)));
 
 
 
 