function plot_gyro_esti(t, gyro_bias, gyro_esti, P)
% 仿真时，画陀螺仪零偏估计结果
% t为时间序列
% gyro_bias为仿真设置的零偏真值序列，deg/s
% gyro_esti为估计的零偏序列，deg/s
% P为理论估计标准差，rad/s

P = P/pi*180;
P = P*3;

%% 画陀螺仪零偏估计
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
title('陀螺仪零偏估计')

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