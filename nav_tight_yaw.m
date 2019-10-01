function nav_tight_yaw()
% 紧组合后处理程序，使用航向角作为量测，包含加速度计零偏，17维模型
% 标$的换数据时需要修改
% 看伪距、伪距率噪声运行receiver_noise.m
% 看相位差测量噪声，看相位差曲线

% ** 先运行att_measure.m
% ** 运行前注意：数据范围、滤波器参数

% 结果：加计零偏估计不稳定，静止时水平姿态可能不准
%       加计零偏估计的问题在于，当动起来时加计零偏可估，P值收敛，但是估计不准，当停下来再启动时，加计零偏是不准的，但增益却已经没了
% 0726数据验证通过

%% 导入数据 ($)
imu_data = evalin('base', 'imu_data'); %IMU数据
ta = evalin('base', 'output_ta(:,1)'); %时间序列
output_sv = evalin('base', 'output_sv_A'); %卫星信息，[x,y,z, rho, vx,vy,vz, rhodot]
BLs = evalin('base', 'BLs'); %姿态测量结果

%% 数据范围 ($)
range = 1:length(ta); %所有点
% range = 1:40000;

%% 输入数据截取
ta = ta(range);
output_sv = output_sv(:,:,range);
BLs = BLs(range,:);
index = find(round(imu_data(:,1)*1e4)==round(ta(1)*1e4),1); %IMU数据索引
imu_data = imu_data(index+(1:length(range))-1,2:7); %删除第一列时间

%% 申请存储空间
% 因为要用统一的画图程序，不需要的存储空间还是要有
n = length(ta); %数据点数

filter_nav    = zeros(n,9);  %滤波器导航输出
filter_bias   = zeros(n,6);  %零偏估计
filter_dtr    = zeros(n,1);  %钟差估计
filter_dtv    = zeros(n,1);  %钟频差估计
filter_tau    = zeros(n,1);  %路径差估计
filter_P      = zeros(n,17); %理论估计标准差
filter_gps    = zeros(n,11); %直接卫星解算，最后三列是航向角、俯仰角和路径差

%% 初始位置
sv = output_sv(:,:,1);
sv(isnan(sv(:,1)),:) = [];
p0 = pos_solve(sv);
lat = p0(1);
lon = p0(2);
h = p0(3);

%% 初始化滤波器 ($)
dt = 0.01;
a = 6371000; %地球半径
para.P = diag([[1,1,1]*1 /180*pi, ...     %初始姿态误差，rad
               [1,1,1]*1, ...             %初始速度误差，m/s
               [1/a,secd(lat)/a,1]*5, ... %初始位置误差，[rad,rad,m]
               2e-8 *3e8, ...             %初始钟差距离，m
               3e-9 *3e8, ...             %初始钟频差速度，m/s
               [1,1,1]*0.2 /180*pi, ...   %初始陀螺仪零偏，rad/s
               [1,1,1]*2e-3 *9.8])^2;     %初始加速度计零偏，m/s^2
para.Q = diag([[1,1,1]*0.15 /180*pi, ...
               ... %姿态一步预测不确定度，rad/s（取陀螺仪噪声标准差）
               [1,1,1]*1.5e-3 *9.8, ...
               ... %速度一步预测不确定度，m/s/s（因为存在零偏，取加速度计噪声标准差的数倍）
               [1/a,secd(lat)/a,1]*1.5e-3 *9.8 *(dt/1), ...
               ... %位置一步预测不确定度，m/s（取速度不确定度的积分或半积分）
               0.03e-9 *3e8 *(dt/1), ...
               ... %钟差距离一步预测不确定度，m/s（取钟频差速度漂移的积分或半积分）
               0.03e-9 *3e8, ...
               ... %钟频差速度漂移，m/s/s（需根据所用晶振和估计曲线精心调节）
               [1,1,1]*0.01 /180*pi, ...
               ... %陀螺仪零偏漂移，rad/s/s（需根据估计曲线精心调节）
               [1,1,1]*0.1e-3 *9.8])^2 * dt^2;
                   %加速度计零偏漂移，m/s^2/s（需根据估计曲线精心调节）
para.R_rho    = 4^2;
para.R_rhodot = 0.04^2;
para.R_psi    = (0.1 /180*pi)^2;
NF = navFilter_tight_open_yaw([lat, lon, h], ...
                    [0, 0, 0] + [1, -1, 0.5]*0, ...
                    [BLs(1,1), BLs(1,2), 0] + [1, -1, 1]*0, ...
                    dt, para);

%% 计算
for k=1:n
    % 更新导航滤波器
    sv = output_sv(:,:,k);
    sv(isnan(sv(:,1)),:) = []; %删除无数据的卫星
    NF = NF.update(imu_data(k,:), sv, BLs(k,1));
    
    % 存储导航结果
    filter_nav(k,1:3) = NF.pos;
    filter_nav(k,4:6) = NF.vel;
    filter_nav(k,7:9) = NF.att;
    filter_bias(k,:) = NF.bias;
    filter_dtr(k) = NF.dtr;
    filter_dtv(k) = NF.dtv;
    filter_gps(k,:) = [pos_solve(sv), BLs(k,[1,2,4])];
    
    % P阵
    filter_P(k,:) = sqrt(diag(NF.Px)');
    Cnb = angle2dcm(NF.att(1)/180*pi, NF.att(2)/180*pi, NF.att(3)/180*pi);
    C = zeros(3);
    C(1,1) = -Cnb(1,3)*Cnb(1,1) / (Cnb(1,1)^2+Cnb(1,2)^2);
    C(1,2) = -Cnb(1,3)*Cnb(1,2) / (Cnb(1,1)^2+Cnb(1,2)^2);
    C(1,3) = 1;
    C(2,1) = -Cnb(1,2) / sqrt(1-Cnb(1,3)^2);
    C(2,2) =  Cnb(1,1) / sqrt(1-Cnb(1,3)^2);
    C(2,3) = 0;
    C(3,1) = (Cnb(2,2)*Cnb(3,3)-Cnb(3,2)*Cnb(2,3)) / (Cnb(3,3)^2+Cnb(2,3)^2);
    C(3,2) = (Cnb(3,1)*Cnb(2,3)-Cnb(2,1)*Cnb(3,3)) / (Cnb(3,3)^2+Cnb(2,3)^2);
    C(3,3) = 0;
    P = C*NF.Px(1:3,1:3)*C';
    filter_P(k,1:3) = sqrt(diag(P)');
end

%% 输出数据
assignin('base', 'filter_nav',  filter_nav)
assignin('base', 'filter_bias', filter_bias)
assignin('base', 'filter_dtr',  filter_dtr)
assignin('base', 'filter_dtv',  filter_dtv)
assignin('base', 'filter_tau',  filter_tau)
assignin('base', 'filter_P',    filter_P)
assignin('base', 'filter_gps',  filter_gps)
filter_ta = ta;
assignin('base', 'filter_ta',   filter_ta)
filter_imu = imu_data;
assignin('base', 'filter_imu',  filter_imu)

%% 画图
figure_tight_result;

end