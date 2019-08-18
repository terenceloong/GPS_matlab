function plot_gnss_file(file_path)
% 画前0.1s的数据（400000个点）

n = 4e5; %0.1s

fileID = fopen(file_path, 'r');
    data = fread(fileID, [2,n], 'int16'); %两行向量
    t = (1:n)/4e6;
    figure
    plot(t, data(1,:)) %实部
    hold on
    plot(t, data(2,:)) %虚部
    xlabel('\itt\rm(s)')
fclose(fileID);

end