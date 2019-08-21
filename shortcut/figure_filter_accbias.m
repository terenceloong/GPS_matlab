t = filter_ta - filter_ta(1);

figure('Position', [400, 80, 650, 800]);

subplot(3,1,1)
plot(t,filter_bias(:,4), 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('x加速度计零偏')

subplot(3,1,2)
plot(t,filter_bias(:,5), 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('y加速度计零偏')

subplot(3,1,3)
plot(t,filter_bias(:,6), 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('z加速度计零偏')