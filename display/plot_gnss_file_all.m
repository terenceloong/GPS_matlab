function plot_gnss_file_all(file_path)
% 浏览文件中的所有数据，每1000个点（0.25ms）画一个点

%% 内存映射版本
% 使用内存映射运行时，内存占用急剧升高，可能是画图语句中需要将全部文件都读进来再取数
% m = memmapfile(file_path, 'Format','int16'); %文件内存映射
% 
% di = 2000; %每隔1000个点取1个，di的数值要乘2
% dt = 0.25e-3; %时间间隔
% n = length(m.Data)/di; %总共画多少个点
% 
% t = (1:n)*dt; %时间轴坐标
% figure
% plot(t, m.Data((1:n)*di-1))
% hold on
% plot(t, m.Data((1:n)*di))
% xlabel('\itt\rm(s)')
% grid on

%% 正常读文件版本
listing = dir(file_path); %获取文件信息
di = 1000; %每隔1000个点取1个
dt = 0.25e-3; %时间间隔
n = listing.bytes/4/di; %总共画多少个点
data = zeros(2,n); %存储待画的点，使用int16存储时画图会报莫名其妙的错

fileID = fopen(file_path, 'r');
for k=1:n
    temp = fread(fileID, [2,di], 'int16');
    data(:,k) = double(temp(:,end));
end
fclose(fileID);

t = (1:n)*dt; %时间轴坐标
figure
plot(t, data(1,:)) %实部
hold on
plot(t, data(2,:)) %虚部
xlabel('\itt\rm(s)')
grid on

end