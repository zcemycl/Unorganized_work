% x_cov 
c_all = [];
for i = 1:107
[c, lg] = xcov(w_all(i,:), 100, 'coeff');
c_all(i,:) = c;
end

% for i = 1:107
%     hold on 
%     plot(c_all(i,:));
%     hold off
% end

c_sum = [];
c_sum = sum(c_all, 2);

num_ind_sample = 1100./c_sum;
plot(num_ind_sample)

index_16 = [];
ind_16 = [];
for i = 1:275
    index_16(i) = 4*i;
    ind_16(i) = w_all(16,4*i);
end



