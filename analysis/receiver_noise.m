function receiver_noise
% 分析接收机伪距、伪距率测量噪声
% 接收机静止，给定接收机参考位置，计算卫星到参考位置的距离和相对速度，与测量值相比较
% 标$的换数据时需要修改

%% 导入数据 ($)
svList = evalin('base', 'svList'); %卫星编号列表
sv_info = evalin('base', 'output_sv_A(:,1:8,:)'); %卫星信息，[x,y,z, rho, vx,vy,vz, rhodot]

%% 参考坐标 ($)
p0 = [45.73104, 126.62481, 209]; %经纬度保留小数点后5位，高度保留到整数
rp = lla2ecef(p0); %ecef

%% 数据范围 ($)
% range = 1:size(sv_info,3); %所有点
range = 1:40000; %从第几个点到第几个点

%% 输入数据截取
sv_info = sv_info(:,:,range);

%% 申请存储空间
n = size(sv_info,3); %数据点数
svN = length(svList); %卫星个数

error_rho = zeros(n,svN); %伪距误差，每一列为一颗卫星
error_rhodot = zeros(n,svN); %伪距率误差

%% 计算
for k=1:n
    rs = sv_info(:,1:3,k);
    rsp = ones(svN,1)*rp - rs;
    rho = sum(rsp.*rsp, 2).^0.5;
    rspu = rsp ./ (rho*[1,1,1]);
    vs = sv_info(:,5:7,k);
    vsp = 0 - vs;
    rhodot = sum(vsp.*rspu, 2);
    error_rho(k,:) = (sv_info(:,4,k) - rho)'; %测量值减计算值
    error_rhodot(k,:) = (sv_info(:,8,k) - rhodot)';
end

%% 变量输出
% assignin('base', 'error_rho', error_rho);
% assignin('base', 'error_rhodot', error_rhodot);

%% 画图
std_rho = zeros(svN,1); %伪距噪声标准差
std_rhodot = zeros(svN,1); %伪距率噪声标准差
for k=1:svN
    x = find(~isnan(error_rho(:,k)));
    y = error_rho(x,k); %提取伪距差不为NaN的点，用来做直线拟合
    if ~isempty(x) %没有值的不画
        figure
        %----画伪距误差
        subplot(2,2,1)
        plot(error_rho(:,k))
        grid on
        set(gca, 'xlim', [0,n])
        hold on
        p = polyfit(x, y, 1); %直线拟合
        error_rho_fit = polyval(p, (1:n)'); %拟合后的点
        plot(error_rho_fit) %画拟合直线
        title(['SV ',num2str(svList(k)),',  column ',num2str(k)]) %卫星编号和列号
        %----画伪距噪声
        subplot(2,2,2)
        plot(error_rho(:,k)-error_rho_fit)
        grid on
        set(gca, 'xlim', [0,n])
        std_rho(k) = std(error_rho(:,k)-error_rho_fit, 'omitnan');
        title(['\sigma=',num2str(std_rho(k),'%.2f')])
        %----画伪距率误差
        subplot(2,2,3)
        plot(error_rhodot(:,k))
        grid on
        set(gca, 'xlim', [0,n])
        std_rhodot(k) = std(error_rhodot(:,k), 'omitnan');
        title(['\mu=',num2str(mean(error_rhodot(:,k),'omitnan'),'%.2f'),',  \sigma=',num2str(std_rhodot(k),'%.3f')])
    end
end
% 把所有卫星噪声标准差画在一起比较
figure
subplot(2,1,1)
bar(1:svN, std_rho)
grid on
title('伪距噪声标准差')
subplot(2,1,2)
bar(1:svN, std_rhodot)
grid on
title('伪距率噪声标准差')

end