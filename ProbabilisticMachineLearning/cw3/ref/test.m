% e
load kos_doc_data.mat
union=[A;B]; % for union of A and B
doc_union=length(unique(union(:,1)));
words_union=sum(union(:,3)); 
unique_words_union=length(unique(union(:,2)));

for m=1:unique_words_union
    c(m)=sum(A(A(:,2)==m,3)); %total count of vocabulary m in A
end

words_ID_B=unique(B(:,2));
for i=1:max(words_ID_B)
    words_count_B(i)=sum(B(B(:,2)==i,3)); %calculate total counts for a word in B
end
predictive_prob_B=(alpha+c(words_ID_B))./(alpha*unique_words_union+sum(c));
l_B=log(predictive_prob_B).*words_count_B(words_ID_B); %calculate log-likelihood for test B
perplexity_B=exp(-l_B/sum(words_count_B)); %calculate per word perplexity for test B