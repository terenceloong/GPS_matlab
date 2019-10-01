classdef navFilter_deep_dphase
    % 深组合导航滤波器，相位差作为量测
    % 1、17维状态变量
    % 2、接收机在滤波器外部校正时钟，惯导零偏在滤波器内部校正
    % 3、相位差需要做双差和水平投影，至少需要3个相位差
    
    % 导航状态、卫星校正数据、惯导零偏在执行update后更新
    properties (Access = public)
        % 导航状态（行向量）
        pos     %位置，[lat,lon,h]，[deg,deg,m]
        vel     %速度，[ve,vn,vd]，m/s
        att     %姿态，[psi,theta,gamma]，deg
        % 输出的卫星校正数据
        dtr     %钟差，s
        dtv     %钟频差，s/s
        % 惯导零偏
        bias    %[gyro,acc]，[deg/s,g]，行向量
        % 滤波器参数
        T       %更新周期，s
        lamda   %波长，m
        bl      %基线长度，m
        Px
        Qx
        cnt     %滤波器运行计数
    end
    
    properties (Access = private)
        % 惯导解算用的变量
        latx    %纬度，rad
        lonx    %经度，rad
        hx      %高度，m
        vx      %速度，[ve;vn;vd]，m/s（列向量）
        qx      %四元数，[q0,q1,q2,q3]（行向量）
    end
    
    methods
        %--------初始化
        function obj = navFilter_deep_dphase(p0, v0, a0, T, lamda, bl, para)
            % p0：初始位置，deg
            % v0：初始速度，m/s
            % a0：初始姿态，deg
            %----设置输出
            obj.pos = p0;
            obj.vel = v0;
            obj.att = a0;
            obj.dtr = 0;
            obj.dtv = 0;
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
            obj.Px = para.P;
            obj.Qx = para.Q;
            obj.cnt = 0;
        end
        
        %--------更新
        function [obj, rho, rhodot] = update(obj, imu, sv, sv_m, sigma)
            % imu：惯导输出，[deg/s, g]，行向量
            % sv = [x,y,z, rho, vx,vy,vz, rhodot]，第4列和第8列在滤波器中用不到
            % sv_m = [rho_m, rhodot_m, dphase_m]，伪距、伪距率、相位差量测量，它们的有效维数可以都不同
            % sigma = [sigma_rho, sigma_rhodot, sigma_dphase]，量测量标准差
            % rho, rhodot：使用滤波后的位置、速度计算的伪距、伪距率
            
            obj.cnt = obj.cnt + 1;
            
            %==========惯导解算============================================%
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
            
            %==========状态方程============================================%
            A = zeros(17);
            A(1:3,12:14) = -Cbn;
            A(4:6,1:3) = [0,-fn(3),fn(2); fn(3),0,-fn(1); -fn(2),fn(1),0];
            A(4:6,15:17) = Cbn;
            A(7:9,4:6) = diag([1/(Rm+h), sec(lat)/(Rn+h), -1]);
            A(10,11) = 1;
            Phi = eye(17) + A*obj.T;
            
            %==========量测维数============================================%
            index1 = find(isnan(sv_m(:,1))==0)';  %伪距量测索引
            index2 = find(isnan(sv_m(:,2))==0)';  %伪距率量测索引
            index3 = find(isnan(sv_m(:,3))==0)';  %相位差量测索引
            n1 = length(index1);                  %伪距量测个数
            n2 = length(index2);                  %伪距率量测个数
            n3 = length(index3);                  %相位差量测个数
            n  = size(sv,1);                      %通道个数
            
            if n1<1 %量测数量不够，只做预测
                P = Phi*obj.Px*Phi' + obj.Qx;
                obj.Px = (P+P')/2;
                X = zeros(17,1);
            else %做量测更新
                %==========量测方程========================================%
                % 1. 计算伪距（所有通道都算，不管有没有数）
                rp = lla2ecef([lat/pi*180, lon/pi*180, h]); %接收机位置矢量（ecef，行向量）
                rs = sv(:,1:3);                             %卫星位置矢量（ecef，行向量）
                rsp = ones(n,1)*rp - rs;                    %卫星指向接收机的位置矢量（此处使用矩阵乘法比repmat更快）
                rho = sum(rsp.*rsp, 2).^0.5;                %计算的伪距
                rspu = rsp ./ (rho*[1,1,1]);                %卫星指向接收机的单位矢量（ecef）
                % 2. 计算伪距率（所有通道都算，不管有没有数）
                Cen = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
                                -sin(lon),           cos(lon),         0;
                       -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];
                vp = v'*Cen;                                %接收机速度矢量（ecef，行向量）
                vs = sv(:,5:7);                             %卫星速度矢量（ecef，行向量）
                vsp = ones(n,1)*vp - vs;                    %接收机相对卫星的速度矢量
                rhodot = sum(vsp.*rspu, 2);                 %计算的伪距率
                % 3. 量测矩阵、量测量、量测噪声方差阵
                f = 1/298.257223563;
                F = [-(Rn+h)*sin(lat)*cos(lon), -(Rn+h)*cos(lat)*sin(lon), cos(lat)*cos(lon);
                     -(Rn+h)*sin(lat)*sin(lon),  (Rn+h)*cos(lat)*cos(lon), cos(lat)*sin(lon);
                       (Rn*(1-f)^2+h)*cos(lat),             0,                 sin(lat)    ];
                HA = rspu*F;
                HB = rspu*Cen'; %各行为地理系下卫星指向接收机的单位矢量
                Ha = HA(index1,:); %取有效的行
                Hb = HB(index2,:);
                U  = HB(index3,:);
                if n3<3 %无法构造相位差量测
                    H = zeros(n1+n2, 17);
                    H(1:n1,7:9) = Ha;
                    H(1:n1,10) = -ones(n1,1);
                    H((n1+1):end,4:6) = Hb;
                    H((n1+1):end,11) = -ones(n2,1);
                    Z = [   rho(index1) - sv_m(index1,1); ... %伪距差（计算减量测）
                         rhodot(index2) - sv_m(index2,2)];    %伪距率差
                    R = diag([sigma(index1,1)', ...
                              sigma(index2,2)'])^2;
                else %可以构造相位差量测
                    E = [U, -ones(n3,1)];
                    B = [-ones(n3-1,1), eye(n3-1)]; %使E的最后一列为0，得到D
                    D = B*E;
                    C = [-D(2:end,3), eye(n3-2)*D(1,3)]; %使D的第三列为0
                    J = C*B; %相位差变换阵，J按行求和为0
                    Hc = cross(U,ones(n3,1)*Cnb(1,:), 2); %沿行求叉乘
                    Hc = Hc * obj.bl/obj.lamda;
                    H = zeros(n1+n2+n3-2, 17);
                    H(1:n1,7:9) = Ha;
                    H(1:n1,10) = -ones(n1,1);
                    H((n1+1):(n1+n2),4:6) = Hb;
                    H((n1+1):(n1+n2),11) = -ones(n2,1);
                    H((n1+n2+1):end,1:3) = J*Hc;
                    Z = [   rho(index1) - sv_m(index1,1); ... %伪距差（计算减量测）
                         rhodot(index2) - sv_m(index2,2); ... %伪距率差
         J*(U*Cbn(:,1)*obj.bl/obj.lamda - sv_m(index3,3))];   %相位差的差
                    R = diag([sigma(index1,1)', ...
                              sigma(index2,2)', ...
                              zeros(1,n3-2)])^2;
                    R((n1+n2+1):end,(n1+n2+1):end) = J*diag(sigma(index3,3)')^2*J';
                end
                
                %==========滤波更新========================================%
                P = obj.Px;
                Q = obj.Qx;
                if n1<4 %卫星少时将钟频差对应的Q设为0，防止P阵对应的值过分增大
                    Q(11,11) = 0;
                end
                if norm(v)<0.2 %静止时，水平加速度计对应的Q设为0
                    Q(15,15) = 0;
                    Q(16,16) = 0;
                end
                
                P = Phi*P*Phi' + Q;
                K = P*H' / (H*P*H'+R);
                X = K*Z;
                P = (eye(length(X))-K*H)*P;
                obj.Px = (P+P')/2;
                
                if norm(v)<0.2 %静止时不估水平加速度计零偏
                    Y = zeros(2,17);
                    Y(1,15) = 1;
                    Y(2,16) = 1;
                    X = X - P*Y'/(Y*P*Y')*Y*X;
                end
            end
            
            %==========导航修正============================================%
            % 修姿态
            if norm(X(1:3))>1e-6
                phi = norm(X(1:3));
                qc = [cos(phi/2), X(1:3)'/phi*sin(phi/2)];
                q = quatmultiply(qc, q);
            end
            obj.qx = quatnormalize(q);
            % 修速度
            obj.vx = v - X(4:6);
            % 修位置
            obj.latx = lat - X(7);
            obj.lonx = lon - X(8);
            obj.hx = h - X(9);
            % 修惯导零偏
            obj.bias(1:3) = obj.bias(1:3) + X(12:14)'/pi*180;
            obj.bias(4:6) = obj.bias(4:6) + X(15:17)'/g;
            
            %==========输出================================================%
            obj.pos = [obj.latx/pi*180, obj.lonx/pi*180, obj.hx]; %deg
            obj.vel = obj.vx'; %m/s
            [r1,r2,r3] = quat2angle(obj.qx);
            obj.att = [r1,r2,r3]/pi*180; %deg
            obj.dtr = X(10)/299792458; %s
            obj.dtv = X(11)/299792458; %s/s
            
            %==========根据滤波结果计算伪距、伪距率==========================%
            lat = obj.latx; %rad
            lon = obj.lonx;
            h = obj.hx;
            v = obj.vx;
            % 计算伪距
            rp = lla2ecef([lat/pi*180, lon/pi*180, h]);
            rs = sv(:,1:3);
            rsp = ones(n,1)*rp - rs;
            rho = sum(rsp.*rsp, 2).^0.5;
            rspu = rsp ./ (rho*[1,1,1]);
            % 计算伪距率
            Cen = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
                            -sin(lon),           cos(lon),         0;
                   -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];
            vp = v'*Cen;
            vs = sv(:,5:7);
            vsp = ones(n,1)*vp - vs;
            rhodot = sum(vsp.*rspu, 2);
            
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