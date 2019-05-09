% my code for xcov
% process cross covariance
[c, lg] = xcov(w_all,53, 'coeff');
% [c, lg] = xcov(Ms_all,53, 'coeff');
c_update = [];

for i = 1:1100
    c_update(:,i) = sum(c(:,...
        (1100*(i-1)+1):1100*i),2)/i;
end

figure(1)
for i = 1:107
    hold on
    plot(c_update(i,:));
    xlim([1,1100]);
    hold off
    
end

figure(2);
c_update_2 = [];
c_update_2 = sum(c);
plot(c_update_2);
% for i = 1:1100
%     c_update(:,i) = sum(c(:,...
%         (1100*(i-1)+1):1100*i),2)/i;
% end
% 
% for i = 1:3

%     c_update(:,i) = sum(c(:,...
%         (3*(i-1)+1):3*i),2)/i;
% end
% 
% for i = 1:3
%     hold on
%     plot(c_update(i,:));
%     xlim([1,3]);
%     hold off
%     
% end


% % High ranking players
% figure(2)
% plot(c_update(1,:));
% hold on
% plot(c_update(5,:));
% plot(c_update(11,:));
% plot(c_update(16,:));
% title('High ranking players')
% legend('1','5','11','16');
% xlim([0,1100]);
% % axis([0, 1100, -400, 400])
% hold off;
% 
% % Low ranking players
% figure(3)
% plot(c_update(105,:));
% hold on
% plot(c_update(75,:));
% plot(c_update(97,:));
% plot(c_update(89,:));
% title('Low ranking players')
% legend('105','75','97','89');
% xlim([0,1100]);
% % axis([0, 1100, -400, 400])
% hold off;
% 
% % Middle ranking players
% figure(4)
% plot(c_update(53,:));
% hold on
% plot(c_update(14,:));
% plot(c_update(64,:));
% plot(c_update(82,:));
% title('Mid ranking players')
% legend('53','14','64','82');
% xlim([0,1100]);
% % axis([0, 1100, -400, 400])
% hold off;




