% difference between algorithms
diff = [];count_no = 1;
for i = 1:107
    if ii_gibbs(i) ~= ii_mess(i)
        diff(count_no) = i; count_no = count_no+1;
    end
end

diff2 = []; count_no_2 = 1;
for i = 1:107
    if ii_gibbs_mean(i) ~= ii_mess(i)
        diff2(count_no_2) = i; count_no_2 = count_no_2+1;
    end
end