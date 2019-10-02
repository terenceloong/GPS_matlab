function dphase_trend()
% 已知基线，看长时间相位差测量趋势对不对
% 处理所有数据
% 只画图不输出
% 先运行姿态测量程序
% 标$的换数据时需要修改

%% 导入数据 ($)
svList = evalin('base', 'svList'); %卫星编号列表
output_pos = evalin('base', 'output_pos(:,1:3)'); %接收机测量的位置，[lat,lon,h]，deg
output_sv = evalin('base', 'output_sv_A(:,1:8,:)'); %卫星信息，[x,y,z, rho, vx,vy,vz, rhodot]
PDm = evalin('base', 'output_dphase_modified'); %整周修正后的相位差

%% 参考基线 ($)
BLr = [7.8, -0.5, 1.3];
Rr = [cosd(BLr(2))*cosd(BLr(1)); cosd(BLr(2))*sind(BLr(1)); -sind(BLr(2))] * BLr(3);
lamda = 299792458 / 1575.42e6; %波长

%% 申请存储空间
n = size(output_pos,1); %数据点数
svN = length(svList); %卫星个数

PDr = NaN(n,svN); %使用参考基线算的相位差

%% 计算
for k=1:n
    pos = output_pos(k,:);
    Cen = dcmecef2ned(pos(1), pos(2));
    rp = lla2ecef(pos);
    rs = output_sv(:,1:3,k);
    rsp = ones(svN,1)*rp - rs;
    rho = sum(rsp.*rsp, 2).^0.5;
    rspu = rsp ./ (rho*[1,1,1]);
    A = rspu * Cen';
    PDr(k,:) = (A*Rr / lamda)';
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
%----实线画测量的相位差
for k=1:svN
    if sum(~isnan(PDm(:,k)))~=0
        plot(PDm(:,k), 'Color',colorTable(k,:), 'LineWidth',1)
        eval('legend_str = [legend_str; string(num2str(svList(k)))];')
    end
end
%----虚线画参考基线算的相位差
for k=1:svN
    if sum(~isnan(PDm(:,k)))~=0
        plot(PDr(:,k), 'Color',colorTable(k,:), 'LineWidth',1, 'LineStyle','--')
    end
end
legend(legend_str)
set(gca,'Xlim',[0,n])
title('实测相位差与参考相位差')

end