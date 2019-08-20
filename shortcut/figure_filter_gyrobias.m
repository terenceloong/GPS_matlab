t = filter_ta - filter_ta(1);

figure('Position', [400, 80, 650, 800]);

subplot(3,1,1)
plot(t,filter_imu(:,1))
hold on
grid on
plot(t,filter_bias(:,1), 'LineWidth',2)
set(gca,'Xlim',[t(1),t(end)])
set(gca,'Ylim',[-1.5,1.5])
title('xÍÓÂİÁãÆ«')

subplot(3,1,2)
plot(t,filter_imu(:,2))
hold on
grid on
plot(t,filter_bias(:,2), 'LineWidth',2)
set(gca,'Xlim',[t(1),t(end)])
set(gca,'Ylim',[-1.5,1.5])
title('yÍÓÂİÁãÆ«')

subplot(3,1,3)
plot(t,filter_imu(:,3))
hold on
grid on
plot(t,filter_bias(:,3), 'LineWidth',2)
set(gca,'Xlim',[t(1),t(end)])
set(gca,'Ylim',[-1.5,1.5])
title('zÍÓÂİÁãÆ«')