% 根据imu输出计算姿态
% 需要静止一段时间，用来校正零偏
% 用IMU数据（变量名imu_data）进行姿态解算，姿态初值都是0，用最前面一段静止数据计算零偏，运行完看nav变量

% 用前面的点计算零偏
m = 1:4000;
imu_data(:,2) = imu_data(:,2) - mean(imu_data(m,2));
imu_data(:,3) = imu_data(:,3) - mean(imu_data(m,3));
imu_data(:,4) = imu_data(:,4) - mean(imu_data(m,4));

n = size(imu_data,1);
nav = zeros(n,3);

q = [1;0;0;0]; %姿态初值，全0

for k=1:n
    wb = imu_data(k,2:4) /180*pi;
    q = q + 0.5*[  0,   -wb(1), -wb(2), -wb(3);
                 wb(1),    0,    wb(3), -wb(2);
                 wb(2), -wb(3),    0,    wb(1);
                 wb(3),  wb(2), -wb(1),    0]*q*0.01;
    q = quatnormalize(q')';
    
    [r1,r2,r3] = quat2angle(q');
    nav(k,:) = [r1,r2,r3] /pi*180;
end