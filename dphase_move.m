function dphase_move()
% 根据深组合导航输出的姿态计算相位差
% 处理所有数据
% 运行前注意基线信息的修改

%% 导入数据 ($)
svList = evalin('base', 'svList'); %卫星编号列表
output_filter = evalin('base', 'output_filter'); %组合导航滤波器输出
output_sv = evalin('base', 'output_sv_A'); %卫星信息，[x,y,z, rho, vx,vy,vz, rhodot]
PDm = evalin('base', 'output_dphase'); %相位差

%% 基线信息 ($)
bl = 1.3; %基线长度，m
lamda = 299792458 / 1575.42e6; %波长

%% 申请存储空间
n = size(output_filter,1); %数据点数
svN = length(svList); %卫星个数

PDa = NaN(n,svN); %使用姿态算的相位差

%% 计算
for k=1:n
    att = output_filter(k,7:9);
    rb = [cosd(att(2))*cosd(att(1)); ...
          cosd(att(2))*sind(att(1)); ...
         -sind(att(2))] * bl;
    pos = output_filter(k,1:3);
    Cen = dcmecef2ned(pos(1), pos(2));
    rp = lla2ecef(pos);
    rs = output_sv(:,1:3,k);
    rsp = ones(svN,1)*rp - rs;
    rho = sum(rsp.*rsp, 2).^0.5;
    rspu = rsp ./ (rho*[1,1,1]);
    A = rspu * Cen';
    PDa(k,:) = (A*rb/lamda)';
end

%% 画图
colorTable = [    0, 0.447, 0.741;
              0.850, 0.325, 0.098;
              0.929, 0.694, 0.125;
              0.494, 0.184, 0.556;
              0.466, 0.674, 0.188;
              0.301, 0.745, 0.933;
              0.635, 0.078, 0.184;
                  0,     0,     1;
                  1,     0,     0;
                  0,     1,     0;
                  0,     0,     0];
figure
hold on
grid on
legend_str = [];
%----虚线画测量的相位差
for k=1:svN
    if sum(~isnan(PDm(:,k)))~=0
        plot(PDm(:,k), 'Color',colorTable(k,:), 'LineWidth',1, 'LineStyle','--')
        eval('legend_str = [legend_str; string(num2str(svList(k)))];')
    end
end
%----实线画基线算的相位差
for k=1:svN
    if sum(~isnan(PDm(:,k)))~=0
        plot(PDa(:,k), 'Color',colorTable(k,:), 'LineWidth',1)
    end
end
legend(legend_str)
set(gca,'Xlim',[0,n])
title('实测相位差与计算相位差')

% 相位差测量误差
dPD = PDm - PDa;
figure
hold on
grid on
for k=1:svN
    if sum(~isnan(PDm(:,k)))~=0
        plot(dPD(:,k), 'Color',colorTable(k,:))
    end
end
legend(legend_str)
set(gca,'Xlim',[0,n])
title('相位差测量误差')

end