t = filter_ta - filter_ta(1);

figure('Position', [400, 80, 650, 800]);

subplot(3,1,1)
plot(t,filter_gps(:,9))
hold on
grid on
plot(t,filter_nav(:,7))
set(gca,'Xlim',[t(1),t(end)])
title('º½Ïò½Ç')

subplot(3,1,2)
plot(t,filter_gps(:,10))
hold on
grid on
plot(t,filter_nav(:,8))
set(gca,'Xlim',[t(1),t(end)])
title('¸©Ñö½Ç')

subplot(3,1,3)
plot(t,filter_nav(:,9))
grid on
set(gca,'Xlim',[t(1),t(end)])
title('¹ö×ª½Ç')