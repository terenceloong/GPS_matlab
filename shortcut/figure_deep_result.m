t = output_ta(:,1) - output_ta(1,1);

%% 位置
figure('Position', [200, 80, 650, 600]);

subplot(3,1,1)
plot(t,output_filter(:,1), 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('纬度')

subplot(3,1,2)
plot(t,output_filter(:,2), 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('经度')

subplot(3,1,3)
plot(t,output_filter(:,3), 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('高度')

%% 速度
figure('Position', [200, 80, 650, 600]);

subplot(3,1,1)
plot(t,output_pos(:,4))
hold on
grid on
plot(t,output_filter(:,4))
set(gca,'Xlim',[t(1),t(end)])
title('北向速度')

subplot(3,1,2)
plot(t,output_pos(:,5))
hold on
grid on
plot(t,output_filter(:,5))
set(gca,'Xlim',[t(1),t(end)])
title('东向速度')

subplot(3,1,3)
plot(t,output_pos(:,6))
hold on
grid on
plot(t,-output_filter(:,6))
set(gca,'Xlim',[t(1),t(end)])
title('天向速度')

%% 姿态
figure('Position', [200, 80, 650, 600]);

subplot(3,1,1)
plot(t,output_Rx(:,1))
hold on
grid on
plot(t,output_filter(:,7))
set(gca,'Xlim',[t(1),t(end)])
title('航向角')

subplot(3,1,2)
plot(t,output_Rx(:,2))
hold on
grid on
plot(t,output_filter(:,8))
set(gca,'Xlim',[t(1),t(end)])
title('俯仰角')

subplot(3,1,3)
plot(t,output_filter(:,9), 'Color',[0.85,0.325,0.098])
grid on
set(gca,'Xlim',[t(1),t(end)])
title('滚转角')

%% 时钟
figure('Position', [200, 80, 650, 600]);

subplot(2,1,1)
plot(t,output_pos(:,7))
grid on
set(gca,'Xlim',[t(1),t(end)])
title('钟差')

subplot(2,1,2)
plot(t,output_df, 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('钟频差')

%% 陀螺仪零偏
figure('Position', [200, 80, 650, 600]);

subplot(3,1,1)
plot(t,output_imu(:,1))
hold on
grid on
plot(t,output_bias(:,1), 'LineWidth',2)
set(gca,'Xlim',[t(1),t(end)])
set(gca,'Ylim',[-1.5,1.5])
title('x陀螺零偏')

subplot(3,1,2)
plot(t,output_imu(:,2))
hold on
grid on
plot(t,output_bias(:,2), 'LineWidth',2)
set(gca,'Xlim',[t(1),t(end)])
set(gca,'Ylim',[-1.5,1.5])
title('y陀螺零偏')

subplot(3,1,3)
plot(t,output_imu(:,3))
hold on
grid on
plot(t,output_bias(:,3), 'LineWidth',2)
set(gca,'Xlim',[t(1),t(end)])
set(gca,'Ylim',[-1.5,1.5])
title('z陀螺零偏')

%% 加速度计零偏
figure('Position', [200, 80, 650, 600]);

subplot(3,1,1)
plot(t,output_bias(:,4), 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('x加速度计零偏')

subplot(3,1,2)
plot(t,output_bias(:,5), 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('y加速度计零偏')

subplot(3,1,3)
plot(t,output_bias(:,6), 'LineWidth',2)
grid on
set(gca,'Xlim',[t(1),t(end)])
title('z加速度计零偏')