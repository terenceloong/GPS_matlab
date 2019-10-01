function vel_yaw()
% 深组合程序运行后，比较速度方向和航向
% 只画图，无输出

%% 导入数据
vn  = evalin('base', 'output_filter(:,4)'); %北向速度
ve  = evalin('base', 'output_filter(:,5)'); %东向速度
yaw = evalin('base', 'output_filter(:,7)'); %航向

n = length(vn);
t = (1:n)'*0.01;

%% 计算
v = sqrt(vn.^2 + ve.^2);
vel_yaw = NaN(n,1);

index = find(v>0.4);
vel_yaw(index) = atan2d(ve(index),vn(index)); %速度方向
yaw_error = vel_yaw - yaw; %速度方向与航向之差
yaw_error = mod(yaw_error+180,360) - 180;

%% 画图
figure
plot(t,v, 'LineWidth',1)
grid on
title('水平速度')

figure
plot(t,vel_yaw)
hold on
grid on
plot(t,yaw)
title('航向')

figure
plot(t,yaw_error)
grid on
title('航向误差')

for k=2:n
    if yaw(k)-yaw(k-1)<-180
        yaw(k:end) = yaw(k:end) + 360;
    elseif yaw(k)-yaw(k-1)>180
        yaw(k:end) = yaw(k:end) - 360;
    end
end

figure
stackedplot(t, [v, yaw_error, yaw]);
grid on

end