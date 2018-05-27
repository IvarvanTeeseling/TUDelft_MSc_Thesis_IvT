clear all; close all; clc;

rr = load('matlab_data_1.mat');
cc = load('matlab_data_2.mat');

N = rr.data1(:,1);

figure(1)
hold on
plot(rr.data1(:,1), rr.data1(:,2),'r')
plot(cc.data1(:,1), cc.data1(:,2),'b')
title('dbdN-N')
legend('R-R','C-C')
hold off

figure(2)
hold on
h(1) = plot(rr.data1(:,1),rr.data1(:,1),'DisplayName','Ivar');
h(2) = plot(cc.data1(:,2),cc.data1(:,1),'DisplayName','Roland');
h(3) = plot(cc.data1(:,2),cc.data1(:,1)*2,'DisplayName','Larissa');
cnt = 0;
for i = 1:3
    if i == 1 || i == 3
        cnt = cnt + 1;
        H(cnt) = h(i);
    end
end
legend([h])
hold off

% figure(2)
% hold on
% plot(cc.data1(:,1), cc.data1(:,3),'r')
% plot(rr.data1(:,1), rr.data1(:,3),'b')
% title('b-N')
% legend('R-R','C-C')
% hold off
% 
% figure(3)
% hold on
% plot(cc.data1(:,3), cc.data1(:,2),'r')
% plot(rr.data1(:,3), rr.data1(:,2),'b')
% title('dbdN-b')
% legend('R-R','C-C')
% hold off