function [filename, gps_week, gps_second] = download_gps_almanac(c)
% 下载GPS历书，自动创建.\temp\almanac文件夹
% c为本地时间，[year, month, day, hour, minute, second]

[gps_week, gps_second] = gps_time(c);

if gps_second<61440
    w = sprintf('%04d',gps_week);
    s = '061440';
elseif gps_second<147456
    w = sprintf('%04d',gps_week);
    s = '147456';
elseif gps_second<233472
    w = sprintf('%04d',gps_week);
    s = '233472';
elseif gps_second<319488
    w = sprintf('%04d',gps_week);
    s = '319488';
elseif gps_second<405504
    w = sprintf('%04d',gps_week);
    s = '405504';
elseif gps_second<589824
    w = sprintf('%04d',gps_week);
    s = '589824';
else
    w = sprintf('%04d',gps_week+1);
    s = '061440';
end
filename = ['.\temp\almanac\',w,'_',s,'.txt'];

% 创建目录
if ~exist('.\temp', 'dir')
    mkdir .\temp
end
if ~exist('.\temp\almanac', 'dir')
    mkdir .\temp\almanac
end

% 如果文件不存在，下载历书
if ~exist(filename, 'file')
    websave(filename, ['http://celestrak.com/GPS/almanac/Yuma/',num2str(c(1)),'/almanac.yuma.week',w,'.',s,'.txt']);
end

end