classdef navFilter
    % 卫星的测量进来之前校正（钟差、钟频差、路径不等长）
    % 因为它们可以在单独运行接收机时估计出来
    % 这些参数在接收机单独运行时通过积分器修正，因为误差是直接量测，噪声大
    % 在加导航滤波器时直接修正，因为这时的误差是通过滤波器估计出来的
    % 惯导的测量进来之后校正
    % 18维模型
    
    properties (Access = public)
        % 导航状态、卫星校正数据、惯导零偏在执行更新后更新
        % 导航状态（行向量）
        pos     %位置，[lat,lon,h]，[deg,deg,m]
        vel     %速度，[ve,vn,vd]，m/s
        att     %姿态，[psi,theta,gamma]，deg
        % 输出的卫星校正数据
        dtr     %钟差，s
        dtv     %钟频差，s/s
        tau     %路径不等长，circ
        % 惯导零偏
        bias    %[gyro,acc]，[deg/s,g]，行向量
        % 滤波器参数
        T       %更新周期，s
        lamda   %波长，m
        bl      %基线长度，m
        Px
        Qx
        R_rho
        R_drho
        R_phase
    end
    
    properties (Access = private)
        % 惯导解算用的变量（列向量）
        latx    %纬度，rad
        lonx    %经度，rad
        hx      %高度，m
        vx      %速度，[ve;vn;vd]，m/s（列向量）
        qx      %四元数，[q0,q1,q2,q3]（行向量）
    end
    
    methods
        %--------初始化
        function obj = navFilter(p0, v0, a0, T, lamda, bl)
            % p0：初始位置，deg
            % v0：初始速度，m/s
            % a0：初始姿态，deg
            %----设置输出
            obj.pos = p0;
            obj.vel = v0;
            obj.att = a0;
            obj.dtr = 0;
            obj.dtv = 0;
            obj.tau = 0;
            obj.bias = [0,0,0,0,0,0];
            %----设置参数
            obj.latx = p0(1)/180*pi;
            obj.lonx = p0(2)/180*pi;
            obj.hx = p0(3);
            obj.vx = v0'; %转成列向量
            obj.qx = angle2quat(a0(1)/180*pi, a0(2)/180*pi, a0(3)/180*pi);
            obj.T = T;
            obj.lamda = lamda;
            obj.bl = bl;
            %----P阵初值
            a = 6371000; %地球半径
            lat = obj.latx; %纬度，rad
            obj.Px = diag([[1,1,1]*1 /180*pi, ...     %初始姿态误差，rad
                           [1,1,1]*1, ...             %初始速度误差，m/s
                           [1/a,sec(lat)/a,1]*5, ...  %初始位置误差，[rad,rad,m]
                           5, ...                     %初始钟差距离，m
                           0.1, ...                   %初始钟频差速度，m/s
                           0.1, ...                   %初始路径不等长载波周数，circ
                           [1,1,1]*0.2 /180*pi, ...   %初始陀螺仪零偏，rad/s
                           [1,1,1]*2 *0.001*9.8])^2;  %初始加速度计零偏，m/s^2
            %----过程噪声方差阵
            obj.Qx = diag([[1,1,1]*0.15 /180*pi, ...                       %陀螺仪噪声，rad/s，看静止时陀螺仪输出得到
                           [1,1,1]*1.5 *0.001*9.8, ...                    %加速度计噪声，m/s^2，看静止时加速度计输出得到
                           [1/a,sec(lat)/a,1]*(T/2)*1.5 *0.001*9.8, ...   %位置漂移，m/s
                           0.01*(T/2), ...                                 %钟差距离漂移，m/s
                           0.01, ...                                       %钟频差速度漂移，m/s/s
                           0.01, ...                                       %路径不等长载波周数漂移，circ/s
                           [1,1,1]*0.02 /180*pi, ...                       %陀螺仪零偏漂移，rad/s/s
                           [1,1,1]*0.2 *0.001*9.8])^2 * T^2;               %加速度计零偏漂移，m/s^2/s
            %----量测噪声方差
            obj.R_rho   = 8^2;                       %伪距噪声方差，m，看直接定位噪声得到
            obj.R_drho  = 0.1^2;                     %伪距率噪声方差，m/s，看直接测速噪声得到
            obj.R_phase = 0.004^2;                   %相位差噪声方差，circ，看相位差曲线得到
        end
        
        %--------更新
        function obj = update(obj, imu, sv)
            % imu惯导输出，[deg/s, g]，行向量
            % sv = [x,y,z, rho, vx,vy,vz, drho, phaseDiff]
            % 正常情况下进来的卫星测量数据（伪距、伪距率、相位差）都应该是无误差的
            % 没有相位差的phaseDiff为NaN
            
            %% 惯导解算
            % 1. 零偏补偿
            imu = imu - obj.bias;
            % 2. 地球参数
            lat = obj.latx;
            lon = obj.lonx;
            h = obj.hx;
            [g, Rm, Rn] = earthPara(lat, h); %重力和曲率半径
            % 3. 姿态解算
            wb = imu(1:3) /180*pi; %deg/s => rad/s
            q = obj.qx; %行向量
            q = q + (0.5*[  0,   -wb(1), -wb(2), -wb(3);
                          wb(1),    0,    wb(3), -wb(2);
                          wb(2), -wb(3),    0,    wb(1);
                          wb(3),  wb(2), -wb(1),    0]*q'*obj.T)';
            % 4. 速度解算
            Cnb = quat2dcm(q);
            Cbn = Cnb';
            fb = imu(4:6)' *g; %列向量
            fn = Cbn*fb;
            v0 = obj.vx;
            v = v0 + (fn + [0;0;g])*obj.T; %列向量
            % 5. 位置解算
            lat = lat + (v0(1)+v(1))/2/(Rm+h)*obj.T;
            lon = lon + (v0(2)+v(2))/2/(Rn+h)*sec(lat)*obj.T;
            h = h - (v0(3)+v(3))/2*obj.T;
            
            %% 状态方程
            A = zeros(18);
            A(1:3,13:15) = -Cbn;
            A(4:6,1:3) = [0,-fn(3),fn(2); fn(3),0,-fn(1); -fn(2),fn(1),0];
            A(4:6,16:18) = Cbn;
            A(7:9,4:6) = diag([1/(Rm+h), sec(lat)/(Rn+h), -1]);
            A(10,11) = 1;
            Phi = eye(18) + A*obj.T;
            
            %% 量测方程
            % 1. 量测维数
            index = find(isnan(sv(:,9))==0)';   %有相位差的行号
            n1 = size(sv,1);                    %伪距、伪距率量测个数
            n2 = length(index);                 %相位差量测个数
            % 2. 计算伪距
            rp = lla2ecef([lat/pi*180, lon/pi*180, h]); %接收机位置矢量（ecef），行向量
            rs = sv(:,1:3);                             %卫星位置矢量  （ecef），行向量
            rsp = ones(n1,1)*rp - rs;                   %卫星指向接收机的位置矢量（此处使用矩阵乘法比repmat更快）
            rho = sum(rsp.*rsp, 2).^0.5;                %计算的伪距
            rspu = rsp ./ (rho*[1,1,1]);                %卫星指向接收机的单位矢量（ecef）
            % 3. 计算伪距率
            Cen = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
                            -sin(lon),           cos(lon),         0;
                   -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];
            vp = v'*Cen;                                %接收机速度矢量（ecef），行向量
            vs = sv(:,5:7);                             %卫星速度矢量  （ecef），行向量
            vsp = ones(n1,1)*vp - vs;                   %接收机相对卫星的速度矢量
            drho = sum(vsp.*rspu, 2);                   %计算的伪距率
            % 4. 量测矩阵
            f = 1/298.257223563;
            F = [-(Rn+h)*sin(lat)*cos(lon), -(Rn+h)*cos(lat)*sin(lon), cos(lat)*cos(lon);
                 -(Rn+h)*sin(lat)*sin(lon),  (Rn+h)*cos(lat)*cos(lon), cos(lat)*sin(lon);
                   (Rn*(1-f)^2+h)*cos(lat),             0,                 sin(lat)    ];
            Ha = rspu*F;
            Hb = rspu*Cen'; %各行为地理系下卫星指向接收机的单位矢量
            Hc = zeros(n2,3);
            U = Hb(index,:);
            for k=1:n2
                Hc(k,:) = cross(U(k,:),Cnb(1,:));
            end
            Hc = Hc * obj.bl/obj.lamda;
            H = zeros(2*n1+n2, 18);
            H(1:n1,7:9) = Ha;
            H(1:n1,10) = -ones(n1,1);
            H((n1+1):(2*n1),4:6) = Hb;
            H((n1+1):(2*n1),11) = -ones(n1,1);
            H((2*n1+1):end,1:3) = Hc;
            H((2*n1+1):end,12) = -ones(n2,1);
            % 5. 量测量
            Z = [ rho - sv(:,4); ...                         %伪距差（计算减量测）
                 drho - sv(:,8); ...                         %伪距率差
                 U*Cbn(:,1)*obj.bl/obj.lamda - sv(index,9)]; %相位差的差
            % 6. 量测噪声方差阵
            R = diag([ones(1,n1)*obj.R_rho, ...
                      ones(1,n1)*obj.R_drho, ...
                      ones(1,n2)*obj.R_phase]);
            
            %% 滤波更新
            P = obj.Px;
            Q = obj.Qx;
            P = Phi*P*Phi' + Q;
            K = P*H' / (H*P*H'+R);
            X = K*Z;
            P = (eye(length(X))-K*H)*P;
            obj.Px = (P+P')/2;
            
            %% 导航修正
            % 修姿态
            if norm(X(1:3))>1e-6
                phi = norm(X(1:3));
                qc = [cos(phi/2), X(1:3)'/phi*sin(phi/2)];
                q = quatmultiply(qc, q);
                obj.qx = quatnormalize(q);
            end
            % 修速度
            obj.vx = v - X(4:6);
            % 修位置
            obj.latx = lat - X(7);
            obj.lonx = lon - X(8);
            obj.hx = h - X(9);
            % 修惯导零偏
            obj.bias(1:3) = obj.bias(1:3) + X(13:15)'/pi*180;
            obj.bias(4:6) = obj.bias(4:6) + X(16:18)'/g;
            
            %% 输出
            obj.pos = [obj.latx/pi*180, obj.lonx/pi*180, obj.hx]; %deg
            obj.vel = obj.vx; %m/s
            [r1,r2,r3] = quat2angle(obj.qx);
            obj.att = [r1,r2,r3]/pi*180; %deg
            obj.dtr = X(10)/299792458; %s
            obj.dtv = X(11)/299792458; %s/s
            obj.tau = X(12); %circ
            
        end %end function
    end %end methods
    
end %end classdef

function [g, Rm, Rn] = earthPara(lat, h)
% 重力模型参见WGS84手册4-1
% 跟gravitywgs84计算结果完全一样
% 纬度单位：rad

% a = 6378137;
% f = 1/298.257223563;

% Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5;
% Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5;

sin_lat_2 = sin(lat)^2;
Rm = 6335439.32729282 / (1-0.00669437999014*sin_lat_2)^1.5;
Rn = 6378137.00000000 / (1-0.00669437999014*sin_lat_2)^0.5;

% w = 7.292115e-5;
% GM = 3.986004418e14;
% re = 9.7803253359;
% rp = 9.8321849378;

% b = (1-f)*a;
% k = b*rp/(a*re)-1;
% m = w*w*a*a*b/GM;
% e2 = f*(2-f);

% b = 6356752.3142;
% k = 0.00193185265241;
% m = 0.00344978650684;
% e2 = 6.69437999014e-3;

% r = re * (1+k*sin(lat)^2) / (1-e2*sin(lat)^2)^0.5;
% g = r * (1 - 2/a*(1+f+m-2*f*sin(lat)^2)*h + 3/a^2*h^2);

r = 9.7803253359 * (1+0.00193185265241*sin_lat_2) / (1-0.00669437999014*sin_lat_2)^0.5;
g = r * (1 - 3.135711885774796e-07*(1.006802597171588-0.006705621329495*sin_lat_2)*h + 7.374516772941995e-14*h^2);
    
end