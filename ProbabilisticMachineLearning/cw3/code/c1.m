% c1
% run a1 first
load kos_doc_data.mat

B_2001 = B(find(B(:,1)==2001),:);
total_B2001 = sum(B_2001(:,3));

diffA_B2001 = setdiff(B_2001(:,2),A(:,2));

unique_idc = unique(B_2001(:,2)); 

dim_uidc = size(unique_idc);
pro_uidc = [];
for i = 1:dim_uidc(1)
    index = find(unique_id==unique_idc(i));
    pro_uidc(i) = probid(index);
end

% Log probability
log_pc = sum(log(pro_uidc'.^B_2001(:,3)));
% with factor
log_pcf= log(factorial(total_B2001))-log(prod(factorial(B_2001(:,3))))+log_pc;
mnpdf(B_2001(:,3)',pro_uidc)