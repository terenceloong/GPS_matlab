t = filter_ta - filter_ta(1);

figure('Position', [400, 80, 650, 800]);

subplot(2,1,1)
plot(t,filter_gps(:,7))
hold on
grid on
plot(t,filter_dtr, 'LineWidth',2)
set(gca,'Xlim',[t(1),t(end)])
title('÷”≤Ó')

subplot(2,1,2)
plot(t,filter_gps(:,8))
hold on
grid on
plot(t,filter_dtv, 'LineWidth',2)
set(gca,'Xlim',[t(1),t(end)])
title('÷”∆µ≤Ó')