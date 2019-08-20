% GPS双天线信号处理程序
% 进行位置解算和姿态解算(姿态解算现在在外面单独做)
% 加载IMU数据，在IMU时间点进行卫星信息测量
% 可以实现接收机时钟闭环和开环两种模式(注释掉时钟反馈校正)
% 卫星信息输出增加载噪比(为了评估测量噪声，可以给最小二乘和滤波器加权值)
% 标*的程序段表示可以单独运行，标#的程序段表示不要修改

clear
clc

%% IMU数据 (#)
imu_data = IMU_parse();
close gcf
close gcf
close gcf

%% 选择文件 (#)
default_path = fileread('.\temp\dataPath.txt'); %数据文件所在默认路径
[file, path] = uigetfile([default_path,'\*.dat'], '选择GPS数据文件'); %文件选择对话框，限制为.dat文件
if file==0
    disp('Invalid file!');
    return
end
if strcmp(file(1:4),'B210')==0
    error('File error!');
end
file_path = [path, file];
file_path_A = [file_path(1:end-5),'1.dat'];
file_path_B = [file_path(1:end-5),'2.dat'];
plot_gnss_file(file_path_A); %显示前0.1s数据
plot_gnss_file(file_path_B);
drawnow

%% 计时开始 (#)
tic

%% 创建日志文件 (#)
fclose('all'); %关闭之前打开的所有文件
logID_A = fopen('.\temp\logA.txt', 'w'); %创建日志文件（时间顺序的日志）
logID_B = fopen('.\temp\logB.txt', 'w');

%% 运行时间
msToProcess = 60*5*1000; %处理总时间
sample_offset = 0*4e6; %抛弃前多少个采样点
sampleFreq = 4e6; %接收机采样频率

%% 参考位置
p0 = [45.730952, 126.624970, 212]; %2A楼顶

%% 数据缓存 (#)
buffBlkNum = 40;                     %采样数据缓存块数量（要保证捕获时存储恰好从头开始）
buffBlkSize = 4000;                  %一个块的采样点数（1ms）
buffSize = buffBlkSize * buffBlkNum; %采样数据缓存大小
buff_A = zeros(2,buffSize);          %采样数据缓存，第一行I，第二行Q
buff_B = zeros(2,buffSize);
buffBlkPoint = 0;                    %数据该往第几块存，从0开始
buffHead = 0;                        %最新数据的序号，buffBlkSize的倍数

%% 获取文件时间 (#)
tf = sscanf(file_path((end-22):(end-8)), '%4d%02d%02d_%02d%02d%02d')'; %数据文件开始采样时间（日期时间数组）
[tw, ts] = gps_time(tf); %tw：GPS周数，ts：GPS周内秒数
ta = [ts,0,0] + sample2dt(sample_offset, sampleFreq); %初始化接收机时间，[s,ms,us]
ta = time_carry(round(ta,2)); %取整

%% 根据历书获取当前可能见到的卫星（*）
% svList = [10;15;20;21;24]; %列向量，为了看方便
svList = gps_constellation(tf, p0);
svN = length(svList);

%% 为每颗可能见到的卫星分配跟踪通道 (#)
channels_A = repmat(GPS_L1_CA_channel_struct(), svN,1);
channels_B = repmat(GPS_L1_CA_channel_struct(), svN,1);
for k=1:svN
    channels_A(k).PRN = svList(k);
    channels_A(k).state = 0; %状态未激活
    channels_B(k).PRN = svList(k);
    channels_B(k).state = 0; %状态未激活
end

%% 预置星历 (#)
ephemeris_file = ['.\temp\ephemeris\',file_path((end-22):(end-8)),'.mat'];
if exist(ephemeris_file, 'file')
    load(ephemeris_file); %星历存在，加载星历文件，星历变量名为ephemeris，星历为列；电离层参数变量名为ion
else
    ephemeris = NaN(26,32); %星历不存在，设置空的星历
    ion = NaN(1,8); %空的电离层参数
end
for k=1:svN
    PRN = svList(k);
    channels_A(k).ephemeris = ephemeris(:,PRN); %为通道的星历赋值
    channels_B(k).ephemeris = ephemeris(:,PRN);
    if ~isnan(ephemeris(1,PRN)) %如果存在某颗卫星的星历，打印日志
        fprintf(logID_A, '%2d: Load ephemeris.\r\n', PRN);
        fprintf(logID_B, '%2d: Load ephemeris.\r\n', PRN);
    end
end

%% 创建跟踪结果存储空间 (#)
% 分配了msToProcess行，每跟踪一次输出一次结果，最后删除多余的行
trackResults_A = repmat(trackResult_struct(msToProcess), svN,1);
trackResults_B = repmat(trackResult_struct(msToProcess), svN,1);
for k=1:svN
    trackResults_A(k).PRN = svList(k);
    trackResults_B(k).PRN = svList(k);
end

%% 接收机状态
receiverState = 0; %接收机状态，0表示未初始化，时间还不对，1表示时间已经校正
deltaFreq = 0; %时钟差，理解为百分比，如果差1e-9，生成1500e6Hz的波会差1.5Hz
dtpos = 10; %定位时间间隔，ms
% tp = [ta(1),0,0]; %tp为下次定位时间
% tp(2) = (floor(ta(2)/dtpos)+1) * dtpos; %加到下个整目标时间
% tp = time_carry(tp); %进位
imu_index = find(imu_data(:,1)>(ta(1)+ta(2)/1e3+ta(3)/1e6), 1); %大于当前接收机时间的imu时间索引
tp = sec2smu(imu_data(imu_index,1)); %取其时间，作为下一定位时间点

%% 创建接收机输出存储空间
% 分配msToProcess/dtpos行，每到时间点输出一次，最后根据接收机状态删除多余的行
nRow = msToProcess/dtpos + 1000; %多分配一些，因为IMU采样频率可能比标称的快
output_ta        = zeros(nRow,2);     %第一列为时间（s），第二列为接收机状态
output_pos       = zeros(nRow,8);     %定位，[位置、速度、钟差、钟频差]
output_sv_A      = zeros(svN,9,nRow); %卫星信息(天线A)，[位置、伪距、速度、伪距率、 ...]，第9列为载噪比
output_sv_B      = zeros(svN,9,nRow); %卫星信息(天线B)，[位置、伪距、速度、伪距率、 ...]
output_df        = zeros(nRow,1);     %修正用的钟频差（滤波后的钟频差）
output_dphase    = NaN(nRow,svN);     %未修正的相位差
no = 1; %指向当前存储行

%% 打开文件，创建进度条 (#)
% 文件A
fileID_A = fopen(file_path_A, 'r');
fseek(fileID_A, round(sample_offset*4), 'bof');
if int64(ftell(fileID_A))~=int64(sample_offset*4)
    error('Sample offset error!');
end
% 文件B
fileID_B = fopen(file_path_B, 'r');
fseek(fileID_B, round(sample_offset*4), 'bof');
if int64(ftell(fileID_B))~=int64(sample_offset*4)
    error('Sample offset error!');
end
% 进度条
f = waitbar(0, ['0s/',num2str(msToProcess/1000),'s']);

%% 信号处理
for t=1:msToProcess %名义上的时间，以采样点数计算
    %% 更新进度条 (#)
    if mod(t,1000)==0 %1s步进
        waitbar(t/msToProcess, f, [num2str(t/1000),'s/',num2str(msToProcess/1000),'s']);
    end
    
    %% 读数据 (#)
    buff_A(:,buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = double(fread(fileID_A, [2,buffBlkSize], 'int16')); %天线A
    buff_B(:,buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = double(fread(fileID_B, [2,buffBlkSize], 'int16')); %天线B
    buffBlkPoint = buffBlkPoint + 1;
    buffHead = buffBlkPoint * buffBlkSize;
    if buffBlkPoint==buffBlkNum
        buffBlkPoint = 0; %缓存从头开始
    end
    
	%% 更新接收机时间
    % 当前最后一个采样的接收机时间
    sampleFreq_real = sampleFreq * (1+deltaFreq); %真实的采样频率
    ta = time_carry(ta + sample2dt(buffBlkSize, sampleFreq_real));
    
    %% 捕获 (#)
    % 每1s的采样点搜索一次
    if mod(t,1000)==0
        for k=1:svN %搜索所有可能见到的卫星
            %====天线A
            if channels_A(k).state==0 %如果通道未激活，捕获尝试激活
                [acqResult, peakRatio] = GPS_L1_CA_acq_one(svList(k), buff_A(:,(end-2*8000+1):end)); %2ms数据捕获
                if ~isempty(acqResult) %成功捕获
                    channels_A(k) = GPS_L1_CA_channel_init(channels_A(k), acqResult, t*buffBlkSize, sampleFreq); %激活通道
                    fprintf(logID_A, '%2d: Acquired at %ds, peakRatio=%.2f\r\n', svList(k), t/1000, peakRatio); %打印捕获日志
                end
            end
            %====天线B
            if channels_B(k).state==0 %如果通道未激活，捕获尝试激活
                [acqResult, peakRatio] = GPS_L1_CA_acq_one(svList(k), buff_B(:,(end-2*8000+1):end)); %2ms数据捕获
                if ~isempty(acqResult) %成功捕获
                    channels_B(k) = GPS_L1_CA_channel_init(channels_B(k), acqResult, t*buffBlkSize, sampleFreq); %激活通道
                    fprintf(logID_B, '%2d: Acquired at %ds, peakRatio=%.2f\r\n', svList(k), t/1000, peakRatio); %打印捕获日志
                end
            end
        end
    end
    
    %% 跟踪 (#)
    for k=1:svN
        %====天线A
        if channels_A(k).state~=0 %如果通道激活，进行跟踪
            while 1
                % 判断是否有完整的跟踪数据
                if mod(buffHead-channels_A(k).trackDataHead,buffSize)>(buffSize/2)
                    break
                end
                % 存跟踪结果（通道参数）
                n = trackResults_A(k).n;
                trackResults_A(k).dataIndex(n,:)    = channels_A(k).dataIndex;
                trackResults_A(k).ts0(n,:)          = channels_A(k).ts0;
                trackResults_A(k).remCodePhase(n,:) = channels_A(k).remCodePhase;
                trackResults_A(k).codeFreq(n,:)     = channels_A(k).codeFreq;
                trackResults_A(k).remCarrPhase(n,:) = channels_A(k).remCarrPhase;
                trackResults_A(k).carrFreq(n,:)     = channels_A(k).carrFreq;
                % 基带处理
                trackDataHead = channels_A(k).trackDataHead;
                trackDataTail = channels_A(k).trackDataTail;
                if trackDataHead>trackDataTail
                    [channels_A(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track(channels_A(k), sampleFreq_real, buffSize, buff_A(:,trackDataTail:trackDataHead), logID_A);
                else
                    [channels_A(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track(channels_A(k), sampleFreq_real, buffSize, [buff_A(:,trackDataTail:end),buff_A(:,1:trackDataHead)], logID_A);
                end
                % 存跟踪结果（跟踪结果）
                trackResults_A(k).I_Q(n,:)          = I_Q;
                trackResults_A(k).disc(n,:)         = disc;
                trackResults_A(k).bitStartFlag(n,:) = bitStartFlag;
                trackResults_A(k).CN0(n,:)          = channels_A(k).CN0;
                trackResults_A(k).carrAcc(n,:)      = channels_A(k).carrAcc;
                trackResults_A(k).Px(n,:)           = sqrt(diag(channels_A(k).Px)')*3;
                trackResults_A(k).n                 = n + 1;
            end
        end
        %====天线B
        if channels_B(k).state~=0 %如果通道激活，进行跟踪
            while 1
                % 判断是否有完整的跟踪数据
                if mod(buffHead-channels_B(k).trackDataHead,buffSize)>(buffSize/2)
                    break
                end
                % 存跟踪结果（通道参数）
                n = trackResults_B(k).n;
                trackResults_B(k).dataIndex(n,:)    = channels_B(k).dataIndex;
                trackResults_B(k).ts0(n,:)          = channels_B(k).ts0;
                trackResults_B(k).remCodePhase(n,:) = channels_B(k).remCodePhase;
                trackResults_B(k).codeFreq(n,:)     = channels_B(k).codeFreq;
                trackResults_B(k).remCarrPhase(n,:) = channels_B(k).remCarrPhase;
                trackResults_B(k).carrFreq(n,:)     = channels_B(k).carrFreq;
                % 基带处理
                trackDataHead = channels_B(k).trackDataHead;
                trackDataTail = channels_B(k).trackDataTail;
                if trackDataHead>trackDataTail
                    [channels_B(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track(channels_B(k), sampleFreq_real, buffSize, buff_B(:,trackDataTail:trackDataHead), logID_B);
                else
                    [channels_B(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track(channels_B(k), sampleFreq_real, buffSize, [buff_B(:,trackDataTail:end),buff_B(:,1:trackDataHead)], logID_B);
                end
                % 存跟踪结果（跟踪结果）
                trackResults_B(k).I_Q(n,:)          = I_Q;
                trackResults_B(k).disc(n,:)         = disc;
                trackResults_B(k).bitStartFlag(n,:) = bitStartFlag;
                trackResults_B(k).CN0(n,:)          = channels_B(k).CN0;
                trackResults_B(k).carrAcc(n,:)      = channels_B(k).carrAcc;
                trackResults_B(k).Px(n,:)           = sqrt(diag(channels_B(k).Px)')*3;
                trackResults_B(k).n                 = n + 1;
            end
        end
    end
    
    %% 检查是否到达定位时间
    dtp = (ta(1)-tp(1)) + (ta(2)-tp(2))/1e3 + (ta(3)-tp(3))/1e6; %当前采样时间与定位时间之差，>=0时表示当前采样时间已经到达或超过定位时间
    
    %% 定位时间已到达
    if dtp>=0
        %% 1.计算卫星位置、速度，测量伪距、伪距率，输出其他需要的信息
        sv_A = NaN(svN,8); %天线A
        sv_B = NaN(svN,8); %天线B
        CN0_A = NaN(svN,1);
        CN0_B = NaN(svN,1);
        for k=1:svN
            if channels_A(k).state==2 %检查所有通道状态，对跟踪到的通道计算卫星信息，[位置、伪距、速度、伪距率]
                dn = mod(buffHead-channels_A(k).trackDataTail+1, buffSize) - 1; %trackDataTail恰好超前buffHead一个时，dn=-1
                dtc = dn / sampleFreq_real; %当前采样时间与跟踪点的时间差
                carrFreq = channels_A(k).carrFreq + 1575.42e6*deltaFreq; %修正后的载波频率
                codeFreq = (carrFreq/1575.42e6+1)*1.023e6; %通过载波频率计算的码频率
                codePhase = channels_A(k).remCodePhase + (dtc-dtp)*codeFreq; %定位点码相位
                ts0 = [floor(channels_A(k).ts0/1e3), mod(channels_A(k).ts0,1e3), 0] + [0, floor(codePhase/1023), mod(codePhase/1023,1)*1e3]; %定位点的码发射时间
                [sv_A(k,:),~] = sv_ecef(channels_A(k).ephemeris, tp, ts0); %根据星历计算卫星[位置、伪距、速度]
                sv_A(k,8) = -carrFreq/1575.42e6*299792458;%载波频率转化为速度
                CN0_A(k) = channels_A(k).CN0; %载噪比
            end
            if channels_B(k).state==2 %检查所有通道状态，对跟踪到的通道计算卫星信息，[位置、伪距、速度、伪距率]
                dn = mod(buffHead-channels_B(k).trackDataTail+1, buffSize) - 1; %trackDataTail恰好超前buffHead一个时，dn=-1
                dtc = dn / sampleFreq_real; %当前采样时间与跟踪点的时间差
                carrFreq = channels_B(k).carrFreq + 1575.42e6*deltaFreq; %修正后的载波频率
                codeFreq = (carrFreq/1575.42e6+1)*1.023e6; %通过载波频率计算的码频率
                codePhase = channels_B(k).remCodePhase + (dtc-dtp)*codeFreq; %定位点码相位
                ts0 = [floor(channels_B(k).ts0/1e3), mod(channels_B(k).ts0,1e3), 0] + [0, floor(codePhase/1023), mod(codePhase/1023,1)*1e3]; %定位点的码发射时间
                [sv_B(k,:),~] = sv_ecef(channels_B(k).ephemeris, tp, ts0); %根据星历计算卫星[位置、伪距、速度]
                sv_B(k,8) = -carrFreq/1575.42e6*299792458;%载波频率转化为速度
                CN0_B(k) = channels_B(k).CN0; %载噪比
            end
        end
        
        %% 2.定位
        % 只使用A天线
        sv_visible = sv_A(~isnan(sv_A(:,1)),:); %提取可见卫星
        pos = pos_solve(sv_visible); %定位，如果不够4颗卫星返回8个NaN
        
        %% 3.计算相位差
        % A相位 - B相位
        dphase = NaN(svN,1); %列向量
        for k=1:svN
            if channels_A(k).state==2 && channels_B(k).state==2 %两个天线都跟踪到该颗卫星
                % 天线A
                dn = mod(buffHead-channels_A(k).trackDataTail+1, buffSize) - 1;
                dtc = dn / sampleFreq_real;
                dt = dtc - dtp;
                phase_A = channels_A(k).remCarrPhase + channels_A(k).carrFreq*dt + 0.5*channels_A(k).carrAcc*dt^2; %载波相位
                % 天线B
                dn = mod(buffHead-channels_B(k).trackDataTail+1, buffSize) - 1;
                dtc = dn / sampleFreq_real;
                dt = dtc - dtp;
                phase_B = channels_B(k).remCarrPhase + channels_B(k).carrFreq*dt + 0.5*channels_B(k).carrAcc*dt^2; %载波相位
                % 相位差
                if channels_A(k).inverseFlag*channels_B(k).inverseFlag==1 %两个天线相位翻转相同
                    dphase(k) = mod((channels_A(k).carrCirc+phase_A)-(channels_B(k).carrCirc+phase_B)    +500,1000) - 500;
                else %两个天线相位翻转不同
                    dphase(k) = mod((channels_A(k).carrCirc+phase_A)-(channels_B(k).carrCirc+phase_B)+0.5+500,1000) - 500;
                end
            end
        end
        
        %% 4.时钟反馈修正
%         if receiverState==1 && ~isnan(pos(7))
%             deltaFreq = deltaFreq + 10*pos(8)*dtpos/1000; %钟频差累加
%             ta = ta - sec2smu(10*pos(7)*dtpos/1000); %时钟修正（可以不用进位，在下次更新时进位）
%         end
        
        %% 5.存储输出
        output_ta(no,1)     = tp(1) + tp(2)/1e3 + tp(3)/1e6;
        output_ta(no,2)     = receiverState;
        output_pos(no,:)    = pos;
        output_sv_A(:,:,no) = [sv_A, CN0_A];
        output_sv_B(:,:,no) = [sv_B, CN0_B];
        output_df(no)       = deltaFreq;
        output_dphase(no,:) = dphase';
        
        %% 6.检查初始化
        if receiverState==0 && ~isnan(pos(7))
            if abs(pos(7))>0.1e-3 %钟差大于0.1ms，修正接收机时间
                ta = ta - sec2smu(pos(7)); %时钟修正
                ta = time_carry(ta);
%                 tp(1) = ta(1); %更新下次定位时间
%                 tp(2) = (floor(ta(2)/dtpos)+1) * dtpos;
%                 tp = time_carry(tp);
                imu_index = find(imu_data(:,1)>(ta(1)+ta(2)/1e3+ta(3)/1e6), 1); %更新下次定位时间
            else %钟差小于0.1ms，初始化结束
                receiverState = 1;
            end
        end
        
        %% 7.更新下次定位时间
%         tp = time_carry(tp + [0,dtpos,0]);
        imu_index = imu_index + 1; %imu数据索引加1，指向下一个
        tp = sec2smu(imu_data(imu_index,1));
        no = no + 1; %指向下一存储位置
    end
    
end

%% 关闭文件，关闭进度条 (#)
fclose(fileID_A);
fclose(fileID_B);
fclose(logID_A);
fclose(logID_B);
close(f);

%% 删除空白数据
for k=1:svN
    trackResults_A(k) = trackResult_clean(trackResults_A(k));
    trackResults_B(k) = trackResult_clean(trackResults_B(k));
end
output_ta(no:end,:)         = [];
output_pos(no:end,:)        = [];
output_sv_A(:,:,no:end)     = [];
output_sv_B(:,:,no:end)     = [];
output_df(no:end,:)         = [];
output_dphase(no:end,:)     = [];
% 删除接收机未初始化时的数据
index = find(output_ta(:,2)==0);
output_ta(index,:)          = [];
output_pos(index,:)         = [];
output_sv_A(:,:,index)      = [];
output_sv_B(:,:,index)      = [];
output_df(index,:)          = [];
output_dphase(index,:)      = [];

%% 打印通道日志（*）
clc
disp('<--------antenna A-------->')
print_log('.\temp\logA.txt', svList);
disp('<--------antenna B-------->')
print_log('.\temp\logB.txt', svList);

%% 保存星历 (#)
% 每次运行完都会保存，有新星历自动添加
for k=1:svN
    PRN = channels_A(k).PRN;
    if isnan(ephemeris(1,PRN)) %星历文件中没有
        if ~isnan(channels_A(k).ephemeris(1))
            ephemeris(:,PRN) = channels_A(k).ephemeris; %插入星历
        elseif ~isnan(channels_B(k).ephemeris(1))
            ephemeris(:,PRN) = channels_B(k).ephemeris; %插入星历
        end
    end
end
save(ephemeris_file, 'ephemeris', 'ion');

%% 画跟踪结果（*）
plot_track_double(sampleFreq, msToProcess, svList, trackResults_A, trackResults_B);

%% 画相位差（*）
plot_svdata(output_dphase, svList, '相位差（未修模糊度）');

%% 清除变量（*）
keepVariables = {...
'sampleFreq'; 'msToProcess';
'p0'; 'tf'; 'svList'; 'svN'; %为了画星座图
'channels_A'; 'trackResults_A'; %A天线跟踪信息
'channels_B'; 'trackResults_B'; %B天线跟踪信息
'output_sv_A'; 'output_sv_B'; %卫星信息
'output_ta'; 'output_pos'; 'output_dphase'; 'output_df'; %输出信息
'ephemeris'; 'ion'; %星历
'imu_data'; %IMU数据
};
clearvars('-except', keepVariables{:})

%% 保存结果 (#)
save .\temp\result_double.mat

%% 计时结束 (#)
toc