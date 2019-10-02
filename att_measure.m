function att_measure()
% 双天线姿态测量程序
% 处理所有数据
% 输出基线测量结果：BLs
% 输出修正整周模糊度的相位差：output_dphase_modified（仍包含路径差）
% 可以修改权值确定方法
% 标$的换数据时需要修改

%% 导入数据 ($)
svList = evalin('base', 'svList'); %卫星编号列表
output_pos = evalin('base', 'output_pos(:,1:3)'); %接收机测量的位置，[lat,lon,h]，deg
output_sv = evalin('base', 'output_sv_A(:,1:8,:)'); %卫星信息，[x,y,z, rho, vx,vy,vz, rhodot]
output_dphase = evalin('base', 'output_dphase'); %原始测量的相位差

%% 基线信息 ($)
bl = 1.3; %大致基线长度
br = 0.02; %基线长度范围
tr = [-5,5]; %俯仰角范围，deg
lamda = 299792458 / 1575.42e6; %波长

% 相位差数值范围
circ_limit = 1000;
circ_half = circ_limit/2;

%% 申请存储空间
n = size(output_pos,1); %数据点数
svN = length(svList); %卫星个数

BLs = NaN(n,4);   %基线测量结果，[航向角、俯仰角、基线长度、路径差]，[deg,deg,m,circ]
PDb = NaN(n,svN); %根据基线算的相位差（理论相位差）
PDm = NaN(n,svN); %整周模糊度修正后的相位差（实测相位差）

%% 计算
N = NaN(svN,1); %所有通道的相位整周误差，相位差修正时减去该值
for k=1:n
    pos = output_pos(k,:);
    Cen = dcmecef2ned(pos(1), pos(2));
    rp = lla2ecef(pos);
    rs = output_sv(:,1:3,k);
    rsp = ones(svN,1)*rp - rs;
    rho = sum(rsp.*rsp, 2).^0.5;
    rspu = rsp ./ (rho*[1,1,1]);
    A = rspu * Cen';
    An = [A/lamda, ones(svN,1)]; %单位矢量转换成载波周数，最后添加全为1的一列
    p = output_dphase(k,:)'; %原始相位差
    pm = p - N; %整周修正后的相位差，如果某卫星没确定整周模糊度，计算后相位差为NaN
    pm = mod(pm + circ_half, circ_limit) - circ_half; %归到0附近
    %----计算基线矢量
    if sum(~isnan(pm))<4 %修正后的相位差数量小于4，不能定姿
        if sum(~isnan(p))>=5 %模糊度搜索
            % 从A、p中取
            index = find(~isnan(p)); %有相位差的索引号
            Ac = A(index,:);
            pc = p(index);
            pc = mod(pc,1); %取小数部分
            %----可以排除掉某些卫星
%             Ac(1,:) = [];
%             pc(1) = [];
            Rx = IAR(Ac, pc, lamda, bl+[-br,br], tr);
        else
            Rx = NaN(4,1);
        end
    else %修正后的相位差数量大于等于4，可以直接定姿
        % 从An、pm中取
        index = find(~isnan(pm)); %有相位差的索引号
        Ac = An(index,:);
        pc = pm(index);
        %----可以排除掉某些卫星
%         Ac(1,:) = [];
%         pc(1) = [];
        %----最小二乘
        % Rx = (Ac'*Ac) \ (Ac'*pc);
        %----加权最小二乘
        W = diag(Ac(:,3).^3); %高度角越高权值越大
        Rx = (Ac'*W*Ac) \ (Ac'*W*pc);
    end
    N = round(p-An*Rx); %重新计算所有通道的整周模糊度，新捕获通道的值会被直接计算，中断通道的值会被销毁
    %----存储结果
    L = norm(Rx(1:3));         %基线长度
    psi = atan2d(Rx(2),Rx(1)); %基线航向角
    theta = -asind(Rx(3)/L);   %基线俯仰角
    BLs(k,:) = [psi,theta,L,Rx(4)];
    PDb(k,:) = (A*Rx(1:3) / lamda)';
    PDm(k,:) = (mod((p-N) + circ_half, circ_limit) - circ_half)';
end

%% 输出数据
assignin('base', 'BLs', BLs)
assignin('base', 'output_dphase_modified', PDm)

%% 画基线测量结果
figure
subplot(3,1,1)
plot(BLs(:,1))
grid on
set(gca,'Xlim',[0,n])
title('航向角')
subplot(3,1,2)
plot(BLs(:,2))
grid on
set(gca,'Xlim',[0,n])
title('俯仰角')
subplot(3,1,3)
plot(BLs(:,3))
grid on
set(gca,'Xlim',[0,n])
set(gca,'Ylim',[bl-br,bl+br])
title('基线长度')

%% 画相位差
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
%----虚线画基线算的相位差
for k=1:svN
    if sum(~isnan(PDm(:,k)))~=0
        plot(PDb(:,k), 'Color',colorTable(k,:), 'LineWidth',1, 'LineStyle','--')
    end
end
legend(legend_str)
set(gca,'Xlim',[0,n])
title('实测相位差与计算相位差')
% 正确情况是实测相位差与计算相位差只差常值路径差

%% 画实测相位差与计算相位差之差 (路径差)
% 最好所有通道重合，并且在一条直线上
% 重合表示路径差一致，呈直线表示路径差不变
figure
hold on
grid on
for k=1:svN
    if sum(~isnan(PDm(:,k)))~=0
        plot(PDm(:,k)-PDb(:,k), 'Color',colorTable(k,:))
    end
end
legend(legend_str)
set(gca,'Xlim',[0,n])
set(gca,'Ylim',[-0.5,0.5])
title('实测相位差 - 计算相位差')

end