function [rho, rhodot, Cen, rspu, A] = cal_rho_rhodot(rs, vs, pos, vel)
% 根据卫星和接收机位置、速度计算相对距离和相对速度
% pos：接收机位置，纬经高，行向量，deg
% vel：接收机速度，地理系，行向量
% rho：相对距离
% rhodot：相对速度，距离增大为正
% Cen：ecef到地理系坐标变换阵
% A：地理系下卫星指向接收机单位矢量

n = size(rs,1); %卫星数量

lat = pos(1) /180*pi;
lon = pos(2) /180*pi;
Cen = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
                -sin(lon),           cos(lon),         0;
       -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];

rp = lla2ecef(pos);
rsp = ones(n,1)*rp - rs;
rho = sum(rsp.*rsp, 2).^0.5;
rspu = rsp ./ (rho*[1,1,1]);

vp = vel * Cen;
vsp = ones(n,1)*vp - vs;
rhodot = sum(vsp.*rspu, 2);

A = rspu * Cen';

end