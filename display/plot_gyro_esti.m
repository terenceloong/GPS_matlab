function plot_gyro_esti(t, gyro_bias, gyro_esti, P)

P = P/pi*180;
P = P*3;

%% »­ÍÓÂÝÒÇÁãÆ«¹À¼Æ
figure
subplot(3,1,1)
plot(t, gyro_esti(:,1), 'LineWidth',2)
hold on
grid on
axis manual
plot(t, gyro_bias(:,1)+P(:,1), 'Color','r', 'LineStyle','--')
plot(t, gyro_bias(:,1)-P(:,1), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\it\epsilon_x\rm(\circ/s)')
title('ÍÓÂÝÒÇÁãÆ«¹À¼Æ')

subplot(3,1,2)
plot(t, gyro_esti(:,2), 'LineWidth',2)
hold on
grid on
axis manual
plot(t, gyro_bias(:,2)+P(:,2), 'Color','r', 'LineStyle','--')
plot(t, gyro_bias(:,2)-P(:,2), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\it\epsilon_y\rm(\circ/s)')

subplot(3,1,3)
plot(t, gyro_esti(:,3), 'LineWidth',2)
hold on
grid on
axis manual
plot(t, gyro_bias(:,3)+P(:,3), 'Color','r', 'LineStyle','--')
plot(t, gyro_bias(:,3)-P(:,3), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\it\epsilon_z\rm(\circ/s)')

end