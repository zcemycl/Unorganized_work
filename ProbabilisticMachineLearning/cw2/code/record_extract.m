function [common,no_1,no_2] = record_extract(G, player_1,player_2)
% player 1 vs 16 record
match_1 = [];
match_16 = [];
match_1 = find(G==player_1);
match_16 = find(G==player_2);

for i = 1:size(match_1)
    if match_1(i) > 1801
        match_1(i) = match_1(i)-1801;
    end
end
for i = 1:size(match_16)
    if match_16(i) > 1801
        match_16(i) = match_16(i)-1801;
    end
end

common = intersect(match_1,match_16);
no_1 = sum(G(common)==player_1);
no_2 = sum(G(common)==player_2);


end