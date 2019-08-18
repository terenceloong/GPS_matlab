function plot_nav_error(t, traj, nav, P)
% 仿真时，画导航误差
% t为时间序列
% traj为真实轨迹
% nav为导航结果
% [lat,lon,h, ve,vn,vd, psi,theta,gamma]
% [deg,deg,m,    m/s,        deg]
% P为理论估计标准差，1~3为姿态角，4~6为速度，7~9为位置
% 失准角到姿态角的标准差变换在外面做

nav_error = nav - traj; %导航误差
nav_error(:,7:9) = mod(nav_error(:,7:9)+180,360)-180; %姿态误差转到±180

P(:,1:3) = P(:,1:3)/pi*180;
P(:,7:8) = P(:,7:8)/pi*180;
P = P*3;

%% 画位置误差
figure
subplot(3,1,1)
plot(t, nav_error(:,1), 'LineWidth',2)
hold on
grid on
axis manual
plot(t,  P(:,7), 'Color','r', 'LineStyle','--')
plot(t, -P(:,7), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\itL\rm(\circ)')
title('位置误差')

subplot(3,1,2)
plot(t, nav_error(:,2), 'LineWidth',2)
hold on
grid on
axis manual
plot(t,  P(:,8), 'Color','r', 'LineStyle','--')
plot(t, -P(:,8), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\lambda(\circ)')

subplot(3,1,3)
plot(t, nav_error(:,3), 'LineWidth',2)
hold on
grid on
axis manual
plot(t,  P(:,9), 'Color','r', 'LineStyle','--')
plot(t, -P(:,9), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\ith\rm(m)')

%% 画速度误差
figure
subplot(3,1,1)
plot(t, nav_error(:,4), 'LineWidth',2)
hold on
grid on
axis manual
plot(t,  P(:,4), 'Color','r', 'LineStyle','--')
plot(t, -P(:,4), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\itv_n\rm(m/s)')
title('速度误差')

subplot(3,1,2)
plot(t, nav_error(:,5), 'LineWidth',2)
hold on
grid on
axis manual
plot(t,  P(:,5), 'Color','r', 'LineStyle','--')
plot(t, -P(:,5), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\itv_e\rm(m/s)')

subplot(3,1,3)
plot(t, nav_error(:,6), 'LineWidth',2)
hold on
grid on
axis manual
plot(t,  P(:,6), 'Color','r', 'LineStyle','--')
plot(t, -P(:,6), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\itv_d\rm(m/s)')

%% 画姿态误差
figure
subplot(3,1,1)
plot(t, nav_error(:,7), 'LineWidth',2)
hold on
grid on
axis manual
plot(t,  P(:,1), 'Color','r', 'LineStyle','--')
plot(t, -P(:,1), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\psi(\circ)')
title('姿态误差')

subplot(3,1,2)
plot(t, nav_error(:,8), 'LineWidth',2)
hold on
grid on
axis manual
plot(t,  P(:,2), 'Color','r', 'LineStyle','--')
plot(t, -P(:,2), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\theta(\circ)')

subplot(3,1,3)
plot(t, nav_error(:,9), 'LineWidth',2)
hold on
grid on
axis manual
plot(t,  P(:,3), 'Color','r', 'LineStyle','--')
plot(t, -P(:,3), 'Color','r', 'LineStyle','--')
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\gamma(\circ)')

end