function plot_svdata(data, svList, title_str)
% 画一组数据，每列对应一颗卫星，添加卫星编号图例
% svList为卫星编号列表，title_str为标题字符串
% 无数据的值为NaN
% 如果某一列全为NaN，该数据不画
% 列号对应固定的颜色，目前只能画11条线，再多需增加颜色

%----颜色表
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
                  0,     0,     0;];

figure
hold on
grid on
legend_str = []; %图例字符串数组，string类型
for k=1:length(svList)
    if sum(~isnan(data(:,k)))~=0
        plot(data(:,k), 'LineWidth',1, 'Color',colorTable(k,:))
        eval('legend_str = [legend_str; string(num2str(svList(k)))];')
    end
end
legend(legend_str)
title(title_str)

end