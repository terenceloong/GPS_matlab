t = filter_ta - filter_ta(1);

figure('Position', [400, 80, 650, 800]);

subplot(3,1,1)
plot(t,filter_gps(:,1))
hold on
grid on
plot(t,filter_nav(:,1))
set(gca,'Xlim',[t(1),t(end)])
title('纬度')

subplot(3,1,2)
plot(t,filter_gps(:,2))
hold on
grid on
plot(t,filter_nav(:,2))
set(gca,'Xlim',[t(1),t(end)])
title('经度')

subplot(3,1,3)
plot(t,filter_gps(:,3))
hold on
grid on
plot(t,filter_nav(:,3))
set(gca,'Xlim',[t(1),t(end)])
title('高度')