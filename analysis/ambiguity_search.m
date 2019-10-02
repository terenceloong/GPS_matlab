function ambiguity_search()
% 所有点做模糊度搜索，看模糊度搜索成功率，同时标定基线
% 不输出结果，只画图
% 可以设置数据范围
% 标$的换数据时需要修改
% 运行前注意：基线长度、数据范围

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

%% 数据范围 ($)
% range = 1:size(output_pos,1); %所有点
range = 1:2000; %从第几个点到第几个点

%% 输入数据截取
output_pos = output_pos(range,:);
output_sv = output_sv(:,:,range);
output_dphase = output_dphase(range,:);

%% 申请存储空间
n = size(output_pos,1); %数据点数
svN = length(svList); %卫星个数

BLs = NaN(n,4); %基线测量结果，[航向角、俯仰角、基线长度、路径差]，[deg,deg,m,circ]

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
    p = output_dphase(k,:)'; %原始相位差
    %----计算基线矢量
    if sum(~isnan(p))>=5
        index = find(~isnan(p)); %有相位差的索引号
        Ac = A(index,:);
        pc = p(index);
        pc = mod(pc,1); %取小数部分
        %----可以排除掉某些卫星
%         Ac(1,:) = [];
%         pc(1) = [];
        Rx = IAR(Ac, pc, lamda, bl+[-br,br], tr);
    else
        Rx = NaN(4,1);
    end
    %----存储结果
    L = norm(Rx(1:3));         %基线长度
    psi = atan2d(Rx(2),Rx(1)); %基线航向角
    theta = -asind(Rx(3)/L);   %基线俯仰角
    BLs(k,:) = [psi,theta,L,Rx(4)];
end

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

end