% 双天线GPS/INS深组合程序

clear
clc

%% 加载IMU数据
imu_data = IMU_parse();
close gcf
close gcf
close gcf

%% 选择文件
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

%% 计时开始
tic

%% 创建日志文件
fclose('all'); %关闭之前打开的所有文件
result_path = fileread('.\temp\resultPath.txt'); %存储结果的路径
logID_A = fopen([result_path,'\logA.txt'], 'w'); %创建日志文件（时间顺序的日志）
logID_B = fopen([result_path,'\logB.txt'], 'w');

%% 运行时间 (!)
msToProcess = 60*10*1000; %处理总时间
sample_offset = 0*4e6; %抛弃前多少个采样点
sampleFreq = 4e6; %接收机采样频率

%% 预设参数 (!)
p0 = [45.730952, 126.624970, 212]; %参考位置
bl = 1.30; %基线长度
br = 0.02; %基线长度范围
tr = [-5,5]; %初始俯仰角范围，deg
circ_limit = 1000; %相位差数值范围
circ_half = circ_limit/2;

%% 滤波器参数 (!)
lat = p0(1);
dt0 = 0.01;
a = 6371000; %地球半径
para.P = diag([[1,1,1]*1 /180*pi, ...     %初始姿态误差，rad
               [1,1,1]*1, ...             %初始速度误差，m/s
               [1/a,secd(lat)/a,1]*5, ... %初始位置误差，[rad,rad,m]
               2e-8 *3e8, ...             %初始钟差距离，m
               3e-9 *3e8, ...             %初始钟频差速度，m/s
               [1,1,1]*0.2 /180*pi, ...   %初始陀螺仪零偏，rad/s
               [1,1,1]*2e-3 *9.8])^2;     %初始加速度计零偏，m/s^2
para.Q = diag([[1,1,1]*0.15 /180*pi, ...
               ... %姿态一步预测不确定度，rad/s（取陀螺仪噪声标准差）
               [1,1,1]*4e-3 *9.8, ...
               ... %速度一步预测不确定度，m/s/s（因为存在零偏，取加速度计噪声标准差的数倍）
               [1/a,secd(lat)/a,1]*4e-3 *9.8 *(dt0/1), ...
               ... %位置一步预测不确定度，m/s（取速度不确定度的积分或半积分）
               0.01e-9 *3e8 *(dt0/1), ...
               ... %钟差距离一步预测不确定度，m/s（取钟频差速度漂移的积分或半积分）
               0.01e-9 *3e8, ...
               ... %钟频差速度漂移，m/s/s（需根据所用晶振和估计曲线精心调节）
               [1,1,1]*0.01 /180*pi, ...
               ... %陀螺仪零偏漂移，rad/s/s（需根据估计曲线精心调节）
               [1,1,1]*0.05e-3 *9.8])^2 * dt0^2;
                   %加速度计零偏漂移，rad/s/s（需根据估计曲线精心调节，给大点）

%% 数据缓存
buffBlkNum = 40;                     %采样数据缓存块数量（要保证捕获时存储恰好从头开始）
buffBlkSize = 4000;                  %一个块的采样点数（1ms）
buffSize = buffBlkSize * buffBlkNum; %采样数据缓存大小
buff_A = zeros(2,buffSize);          %采样数据缓存，第一行I，第二行Q
buff_B = zeros(2,buffSize);
buffBlkPoint = 0;                    %数据该往第几块存，从0开始
buffHead = 0;                        %最新数据的序号，buffBlkSize的倍数

%% 获取文件时间
tf = sscanf(file_path((end-22):(end-8)), '%4d%02d%02d_%02d%02d%02d')'; %数据文件开始采样时间（日期时间数组）
[tw, ts] = gps_time(tf); %tw：GPS周数，ts：GPS周内秒数
ta = [ts,0,0] + sample2dt(sample_offset, sampleFreq); %初始化接收机时间，[s,ms,us]
ta = time_carry(round(ta,2)); %取整

%% 根据历书获取当前可能见到的卫星
svList = gps_constellation(tf, p0); %列向量，为了看方便
svN = length(svList); %通道数量

%% 为每颗可能见到的卫星分配跟踪通道
channels_A = repmat(GPS_L1_CA_channel_struct(), svN,1);
channels_B = repmat(GPS_L1_CA_channel_struct(), svN,1);
for k=1:svN
    channels_A(k).PRN = svList(k);
    channels_A(k).state = 0; %状态未激活
    channels_B(k).PRN = svList(k);
    channels_B(k).state = 0; %状态未激活
end

%% 预置星历
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

%% 创建跟踪结果存储空间
% 分配了msToProcess行，每跟踪一次输出一次结果，最后删除多余的行
trackResults_A = repmat(trackResult_struct(msToProcess), svN,1);
trackResults_B = repmat(trackResult_struct(msToProcess), svN,1);
for k=1:svN
    trackResults_A(k).PRN = svList(k);
    trackResults_B(k).PRN = svList(k);
end

%% 接收机状态
lamda = 299792458 / 1575.42e6; %波长
code_length = 299792458 / 1.023e6; %码长
receiverState = 0; %接收机状态，0表示未初始化，时间还不对，1表示时间已经校正，2表示导航滤波器启动
deltaFreq = 0; %时钟差，理解为百分比，如果差1e-9，生成1500e6Hz的波会差1.5Hz
deltaPath = 0; %两天线路径差
dtpos = 10; %定位时间间隔，ms
imu_index = find(imu_data(:,1)>(ta(1)+ta(2)/1e3+ta(3)/1e6), 1); %大于当前接收机时间的imu时间索引
tp = sec2smu(imu_data(imu_index,1)); %取其时间，作为下一定位时间点
chSign = zeros(svN,1); %双天线通道符号，同相为0，反向为0.5
track_index0 = ones(svN,1); %存储上一定位时刻跟踪索引，用来计算码鉴相器均值
dphase_mask = NaN(svN,1); %用来识别哪些通道已经完成了整周模糊度确定，NaN为未确定，0为已确定
lla = NaN(1,3); %接收机位置
ele = NaN(svN,1); %卫星高度角，deg
azi = NaN(svN,1); %卫星方位角，deg
rhodot0 = NaN(svN,1); %卫星上一定位时刻伪距率，用来求伪距率的加速度

%% 创建接收机输出存储空间
nRow = msToProcess/dtpos + 1000; %多分配一些，因为IMU采样频率可能比标称的快
no = 1; %指向当前存储行
%----接收机的输出
output_ta        = NaN(nRow,2);     %第一列为时间（s），第二列为接收机状态
output_pos       = NaN(nRow,8);     %定位，[位置、速度、钟差、钟频差]
output_sv_A      = NaN(svN,8,nRow); %卫星信息(天线A)，[位置、伪距、速度、伪距率]
output_sv_B      = NaN(svN,8,nRow); %卫星信息(天线B)，[位置、伪距、速度、伪距率]
output_df        = NaN(nRow,1);     %修正用的钟频差
output_dp        = NaN(nRow,1);     %修正用的路径差
output_dphase    = NaN(nRow,svN);   %相位差
output_Rx        = NaN(nRow,4);     %基线测量
%----导航滤波器的输出
output_filter    = NaN(nRow,9);     %滤波器导航结果，[位置、速度、姿态]
output_bias      = NaN(nRow,6);     %IMU零偏估计
output_imu       = NaN(nRow,3);     %保存陀螺仪输出，为了画图，看陀螺仪零偏估计得准不准
output_P         = NaN(nRow,size(para.P,1));    %滤波器P阵
output_svn       = NaN(nRow,3);     %滤波器量测数量，[伪距数量、伪距率数量、相位差数量]

%% 打开文件，创建进度条
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
    %----------更新进度条--------------------------------------------------%
    if mod(t,1000)==0 %1s步进
        waitbar(t/msToProcess, f, [num2str(t/1000),'s/',num2str(msToProcess/1000),'s']);
    end
    
    %----------读数据------------------------------------------------------%
    buff_A(:,buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = double(fread(fileID_A, [2,buffBlkSize], 'int16')); %天线A
    buff_B(:,buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = double(fread(fileID_B, [2,buffBlkSize], 'int16')); %天线B
    buffBlkPoint = buffBlkPoint + 1;
    buffHead = buffBlkPoint * buffBlkSize;
    if buffBlkPoint==buffBlkNum
        buffBlkPoint = 0; %缓存从头开始
    end
    
    %----------更新接收机时间----------------------------------------------%
    % 当前最后一个采样的接收机时间
    sampleFreq_real = sampleFreq * (1+deltaFreq); %真实的采样频率
    ta = time_carry(ta + sample2dt(buffBlkSize, sampleFreq_real));
    
    %----------捕获--------------------------------------------------------%
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
        end %end for k=1:svN
    end %end if mod(t,1000)==0
    
    %----------跟踪--------------------------------------------------------%
    for k=1:svN
        % 记录两天线的比特开始标志（跟踪到比特开始的最前面一段）
        % 未激活/不是比特开始时，标志为0，检测到比特开始，标志非0
        % 因为两天线离得较近，会同时跟踪到一个码
        % 当两天线的比特标志都不为0时，根据当前I路数据的符号判断两天线是否存在180度相位翻转
        % 当信号失锁时，比特是错的，此时判断结果无意义，不影响，因为弱信号不会进行相位差计算
        bitStartFlag_A = 0;
        bitStartFlag_B = 0;
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
                        GPS_L1_CA_track_deep(channels_A(k), sampleFreq_real, buffSize, buff_A(:,trackDataTail:trackDataHead), logID_A);
                else
                    [channels_A(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track_deep(channels_A(k), sampleFreq_real, buffSize, [buff_A(:,trackDataTail:end),buff_A(:,1:trackDataHead)], logID_A);
                end
                % 记录比特开始标志
                bitStartFlag_A = bitStartFlag;
                I_P_A = I_Q(1);
                % 存跟踪结果（跟踪结果）
                trackResults_A(k).I_Q(n,:)          = I_Q;
                trackResults_A(k).disc(n,:)         = disc;
                trackResults_A(k).bitStartFlag(n,:) = bitStartFlag;
                trackResults_A(k).CN0(n,1)          = channels_A(k).CN0;
                trackResults_A(k).CN0(n,2)          = channels_A(k).CN0i;
                trackResults_A(k).carrAcc(n,:)      = channels_A(k).carrAcc;
                trackResults_A(k).strength (n,:)    = channels_A(k).strength;
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
                        GPS_L1_CA_track_deep(channels_B(k), sampleFreq_real, buffSize, buff_B(:,trackDataTail:trackDataHead), logID_B);
                else
                    [channels_B(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track_deep(channels_B(k), sampleFreq_real, buffSize, [buff_B(:,trackDataTail:end),buff_B(:,1:trackDataHead)], logID_B);
                end
                % 记录比特开始标志
                bitStartFlag_B = bitStartFlag;
                I_P_B = I_Q(1);
                % 存跟踪结果（跟踪结果）
                trackResults_B(k).I_Q(n,:)          = I_Q;
                trackResults_B(k).disc(n,:)         = disc;
                trackResults_B(k).bitStartFlag(n,:) = bitStartFlag;
                trackResults_B(k).CN0(n,1)          = channels_B(k).CN0;
                trackResults_B(k).CN0(n,2)          = channels_B(k).CN0i;
                trackResults_B(k).carrAcc(n,:)      = channels_B(k).carrAcc;
                trackResults_B(k).strength(n,:)     = channels_B(k).strength;
                trackResults_B(k).n                 = n + 1;
            end
        end
        %----判断通道符号
        % 处于连续跟踪的通道不改变通道符号，因为通道符号可能有误判
        % 进入强信号状态之前的符号判断一定是正确的
        if isnan(dphase_mask(k))
            if bitStartFlag_A~=0 && bitStartFlag_B~=0
                if I_P_A*I_P_B>=0
                    chSign(k) = 0;
                else
                    chSign(k) = 0.5;
                end
            end
        end
    end %end for k=1:svN
    
    %----------定位-------------------------------------------------------%
    dtp = (ta(1)-tp(1)) + (ta(2)-tp(2))/1e3 + (ta(3)-tp(3))/1e6; %当前采样时间与定位时间之差，>=0时表示当前采样时间已经到达或超过定位时间
    if dtp>=0
        
        % 1.计算卫星位置、速度，测量伪距、伪距率、相位差，计算滤波器用的数据
        sv_A = NaN(svN,8); %卫星信息，[位置、伪距、速度、伪距率]
        sv_B = NaN(svN,8);
        dphase = NaN(svN,1); %相位差，circ
        rho_m         = NaN(svN,1); %进入滤波器的伪距量测，m
        rhodot_m      = NaN(svN,1); %进入滤波器的伪距率量测，m/s
        sigma_rho_A   = NaN(svN,1); %伪距量测标准差
        sigma_phase_A = NaN(svN,1); %A天线载波相位标准差
        sigma_phase_B = NaN(svN,1); %B天线载波相位标准差
        for k=1:svN
            %====天线A
            if channels_A(k).state==2
                dn = mod(buffHead-channels_A(k).trackDataTail+1, buffSize) - 1; %trackDataTail恰好超前buffHead一个时，dn=-1
                dtc = dn / sampleFreq_real; %当前采样时间与跟踪点的时间差
                dt = dtc - dtp; %定位点到跟踪点的时间差
                carrFreq = channels_A(k).carrFreq + 1575.42e6*deltaFreq; %修正后的载波频率
                codeFreq = (carrFreq/1575.42e6+1)*1.023e6; %通过载波频率计算的码频率
                codePhase = channels_A(k).remCodePhase + dt*codeFreq; %定位点码相位
                ts0 = [floor(channels_A(k).ts0/1e3), mod(channels_A(k).ts0,1e3), 0] + [0, floor(codePhase/1023), mod(codePhase/1023,1)*1e3]; %定位点的码发射时间
                [sv_A(k,:),~] = sv_ecef(channels_A(k).ephemeris, tp, ts0); %根据星历计算卫星[位置、伪距、速度]
                sv_A(k,8) = -carrFreq*lamda; %载波频率转化为速度
                sv_A(k,8) = sv_A(k,8) + channels_A(k).ephemeris(9)*299792458; %修卫星钟频差，卫星钟快测的伪距率偏小
                phase_A = channels_A(k).remCarrPhase + channels_A(k).carrNco*dt;
                % 电离层延迟校正
                if ~isnan(ion(1)) %存在电离层参数
                    if ~isnan(ele(k)) && ~isnan(lla(1))
                        tiono = Klobuchar_iono(ion, ele(k), azi(k), lla(1), lla(2), tp(1)+tp(2)/1e3+tp(3)/1e6); %计算电离层延时
                        sv_A(k,4) = sv_A(k,4) - tiono*299792458; %修正伪距
                    end
                end
                %---------------------------------------------------------%
                if channels_A(k).strength==2 %强信号时给伪距率量测
                    rhodot_m(k) = sv_A(k,8);
                    sigma_phase_A(k) = sqrt(channels_A(k).carrStd.D0);
                end
                if channels_A(k).CN0>35 %平均载噪比大于35时给伪距量测
                    % 计算最近更新间隔内码鉴相器输出的均值
                    % 本地码超前，测的伪距偏短，码鉴相器为负值（因为负斜率），修正伪距是减码鉴相器的值
                    rho_m(k) = sv_A(k,4) - mean(trackResults_A(k).disc(track_index0(k):(trackResults_A(k).n-1),1))*code_length;
                    sigma_rho_A(k) = sqrt(channels_A(k).codeStd.D0)/3.2 * code_length;
                end
                %---------------------------------------------------------%
            end
            track_index0(k) = trackResults_A(k).n; %更新跟踪索引
            %====天线B
            if channels_B(k).state==2
                dn = mod(buffHead-channels_B(k).trackDataTail+1, buffSize) - 1; %trackDataTail恰好超前buffHead一个时，dn=-1
                dtc = dn / sampleFreq_real; %当前采样时间与跟踪点的时间差
                dt = dtc - dtp; %定位点到跟踪点的时间差
                carrFreq = channels_B(k).carrFreq + 1575.42e6*deltaFreq; %修正后的载波频率
                codeFreq = (carrFreq/1575.42e6+1)*1.023e6; %通过载波频率计算的码频率
                codePhase = channels_B(k).remCodePhase + dt*codeFreq; %定位点码相位
                ts0 = [floor(channels_B(k).ts0/1e3), mod(channels_B(k).ts0,1e3), 0] + [0, floor(codePhase/1023), mod(codePhase/1023,1)*1e3]; %定位点的码发射时间
                [sv_B(k,:),~] = sv_ecef(channels_B(k).ephemeris, tp, ts0); %根据星历计算卫星[位置、伪距、速度]
                sv_B(k,8) = -carrFreq*lamda; %载波频率转化为速度
                sv_B(k,8) = sv_B(k,8) + channels_B(k).ephemeris(9)*299792458; %修卫星钟频差，卫星钟快测的伪距率偏小
                phase_B = channels_B(k).remCarrPhase + channels_B(k).carrNco*dt;
                %---------------------------------------------------------%
                if channels_A(k).state==2 && channels_A(k).strength==2 && channels_B(k).strength==2 %A、B都为强信号时计算相位差
                    dphase(k) = mod((channels_A(k).carrCirc+phase_A)-(channels_B(k).carrCirc+phase_B)+chSign(k)+circ_half,circ_limit) - circ_half;
                    dphase(k) = dphase(k) - deltaPath; %路径差修正
                    sigma_phase_B(k) = sqrt(channels_B(k).carrStd.D0);
                end
                %---------------------------------------------------------%
            end
        end
        dphase_m = dphase + dphase_mask; %确定了模糊度的相位差
        
        % 2.直接定位
%         pos = pos_solve(sv_A(~isnan(sv_A(:,1)),:)); %提取可见卫星定位，如果不够4颗卫星返回8个NaN
        pos = pos_solve(sv_A(~isnan(rho_m),:)); %提取信号有一定强度的卫星进行定位
        lla = pos(1:3); %接收机位置
        
        % a.取两个天线接收卫星的并集
        rs = NaN(svN,3);
        vs = NaN(svN,3);
        index = find(~isnan(sv_B(:,1)));
        rs(index,:) = sv_B(index,1:3);
        vs(index,:) = sv_B(index,5:7);
        index = find(~isnan(sv_A(:,1))); %以A为主
        rs(index,:) = sv_A(index,1:3);
        vs(index,:) = sv_A(index,5:7);
        
        % 3.更新导航滤波器
        if receiverState==2
            %----给量测噪声
            sigma_rhodot_A = sigma_phase_A*6 * lamda; %25H带宽情况下频率噪声是输入相位噪声的6倍（仿真得到）
            sigma_dphase = sqrt(sigma_phase_A.^2+sigma_phase_B.^2)/4.5; %25Hz带宽情况下实际相位噪声是输入相位噪声的1/4.5倍（仿真得到）
            sigma_dphase = sigma_dphase + 0.04*(90-ele)/90; %相位噪声加一个与高度角相关的基值
            %----更新导航滤波器
            [NF, rho, rhodot] = NF.update(imu_data(imu_index,2:7), sv_A, [rho_m,rhodot_m,dphase_m], [sigma_rho_A,sigma_rhodot_A,sigma_dphase]);
%             [NF, rho, rhodot] = NF.update(imu_data(imu_index,2:7)+[0,0,0.6,0,0,0], sv_A, [rho_m,rhodot_m,dphase_m*NaN], [sigma_rho_A,sigma_rhodot_A,sigma_dphase]);
            delta_rho_A = rho - sv_A(:,4); %码环伪距误差，计算值减环路值
            delta_rhodot_A = rhodot - sv_A(:,8); %载波环伪距率误差
            lla = NF.pos; %接收机位置
            %----更新A天线通道载波计数，修正整周误差
            rb = [cosd(NF.att(2))*cosd(NF.att(1)), ...
                  cosd(NF.att(2))*sind(NF.att(1)), ...
                 -sind(NF.att(2))] * bl; %地理系下基线矢量
            lat = lla(1) /180*pi;
            lon = lla(2) /180*pi;
            Cen = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
                            -sin(lon),           cos(lon),         0;
                   -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];
            rp = lla2ecef(lla);
            rsp = ones(svN,1)*rp - rs;
            rho = sum(rsp.*rsp, 2).^0.5;
            rspu = rsp ./ (rho*[1,1,1]);
            A = rspu * Cen';
            dphase_cal = A*rb' / lamda; %用姿态算的相位差
            dphase_error = dphase - dphase_cal; %相位差误差
            dN = round(dphase_error); %通道整周修正量，相位差误差取整，如果没有周跳，相位差应在0附近，取整后为0
            dN(isnan(dN)) = 0; %无相位差的通道不修正
            dphase_mask( isnan(dphase)) = NaN;
            dphase_mask(~isnan(dphase)) = 0;
            if NF.cnt>200 %滤波器稳定后修正整周误差
                ki = find(dN~=0)';
                for k=ki
                    channels_A(k).carrCirc = mod(channels_A(k).carrCirc-dN(k), circ_limit);
                end
            end
            %----计算卫星高度角、方位角、卫星运动引起的伪距率、载波加速度
            ele = asind(A(:,3));
            azi = atan2d(-A(:,2),-A(:,1));
            rhodot1 = sum(-vs.*rspu, 2);
            carrAcc = -(rhodot1-rhodot0) / (dtpos/1000) / lamda;
            rhodot0 = rhodot1;
            %----计算B天线环路误差
            % 通过基线矢量计算的B天线位置，误差与姿态误差有关，不大
            % 通过角速度和基线矢量计算的B天线速度与陀螺仪噪声有关，0.4deg/s、1m基线长度，对应速度误差0.007m/s
            rb = rb * Cen; %ecef下基线矢量
            rp = rp + rb; %ecef下B天线位置矢量
            rsp = ones(svN,1)*rp - rs;
            rho = sum(rsp.*rsp, 2).^0.5;
            rspu = rsp ./ (rho*[1,1,1]);
            Cnb = angle2dcm(NF.att(1)/180*pi, NF.att(2)/180*pi, NF.att(3)/180*pi);
            vb = cross((imu_data(imu_index,2:4)-NF.bias(1:3))/180*pi, [bl,0,0]) * Cnb * Cen; %ecef下杆臂速度矢量
            vp = NF.vel*Cen + vb; %ecef下B天线速度矢量
            vsp = ones(svN,1)*vp - vs;
            rhodot = sum(vsp.*rspu, 2);
            delta_rho_B = rho - sv_B(:,4); %码环伪距误差，计算值减环路值
            delta_rhodot_B = rhodot - sv_B(:,8); %载波环伪距率误差
            %----修正通道（通道需要自行进入到状态2）
            for k=1:svN
                %====天线A
                if channels_A(k).state==2
                    % 不论信号强弱，全修码相位（码超前，测量的伪距偏短，delta_rho为正）
                    channels_A(k).remCodePhase = channels_A(k).remCodePhase - delta_rho_A(k)/code_length;
                    % 弱信号修载波频率，顺便更新载波环积分器上下界（本地载波快，测量的伪距率偏小，delta_rhodot为正）
                    if channels_A(k).strength~=2
                        channels_A(k).PLL.Int = channels_A(k).PLL.Int - delta_rhodot_A(k)/lamda;
                        channels_A(k).PLL.upper = channels_A(k).PLL.Int + 1;
                        channels_A(k).PLL.lower = channels_A(k).PLL.Int - 1;
                    else
                    % 强信号只更新载波环积分器上下界
                        carrFreq0 = channels_A(k).PLL.Int - delta_rhodot_A(k)/lamda; %修正后的载波频率，不用来更新载波环积分器，只用来计算积分器上下界
                        channels_A(k).PLL.upper = carrFreq0 + 1;
                        channels_A(k).PLL.lower = carrFreq0 - 1;
                    end
                    % 更新通道载波加速度（如果没有载波加速度前馈，二阶环跟踪频率斜坡，频率有静差）
                    if channels_A(k).trackStage=='D'
                        channels_A(k).carrAcc = carrAcc(k);
                    end
                    % 使刚进入2状态的通道进入深组合跟踪模式
                    if channels_A(k).trackStage=='T'
                        channels_A(k).trackStage = 'D';
                    end
                end
                %====天线B
                if channels_B(k).state==2
                    % 不论信号强弱，全修码相位（码超前，测量的伪距偏短，delta_rho为正）
                    channels_B(k).remCodePhase = channels_B(k).remCodePhase - delta_rho_B(k)/code_length;
                    % 弱信号修载波频率，顺便更新载波环积分器上下界（本地载波快，测量的伪距率偏小，delta_rhodot为正）
                    if channels_B(k).strength~=2
                        channels_B(k).PLL.Int = channels_B(k).PLL.Int - delta_rhodot_B(k)/lamda;
                        channels_B(k).PLL.upper = channels_B(k).PLL.Int + 1;
                        channels_B(k).PLL.lower = channels_B(k).PLL.Int - 1;
                    else
                    % 强信号只更新载波环积分器上下界
                        carrFreq0 = channels_B(k).PLL.Int - delta_rhodot_B(k)/lamda; %修正后的载波频率，不用来更新载波环积分器，只用来计算积分器上下界
                        channels_B(k).PLL.upper = carrFreq0 + 1;
                        channels_B(k).PLL.lower = carrFreq0 - 1;
                    end
                    % 更新通道载波加速度（如果没有载波加速度前馈，二阶环跟踪频率斜坡，频率有静差）
                    if channels_B(k).trackStage=='D'
                        channels_B(k).carrAcc = carrAcc(k);
                    end
                    % 使刚进入2状态的通道进入深组合跟踪模式
                    if channels_B(k).trackStage=='T'
                        channels_B(k).trackStage = 'D';
                    end
                end
            end
            %----存储滤波器输出
            output_filter(no,:) = [NF.pos, NF.vel, NF.att];
            output_bias(no,:)   = NF.bias;
            output_imu(no,:)    = imu_data(imu_index,2:4);
            output_P(no,:)      = sqrt(diag(NF.Px)');
            output_svn(no,1)    = sum(~isnan(rho_m));
            output_svn(no,2)    = sum(~isnan(rhodot_m));
            output_svn(no,3)    = sum(~isnan(dphase_m));
        end
        
        % 4.直接测姿
        if receiverState~=2 %正常模式进行模糊度搜索，修正通道整周误差
            lat = lla(1) /180*pi;
            lon = lla(2) /180*pi;
            Cen = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
                            -sin(lon),           cos(lon),         0;
                   -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];
            rp = lla2ecef(lla);
            rsp = ones(svN,1)*rp - rs;
            rho = sum(rsp.*rsp, 2).^0.5;
            rspu = rsp ./ (rho*[1,1,1]);
            A = rspu * Cen';
            ele = asind(A(:,3)); %卫星高度角
            azi = atan2d(-A(:,2),-A(:,1)); %卫星方位角
            An = [A/lamda, ones(svN,1)]; %单位矢量转换成载波周数，最后添加全为1的一列
            if sum(~isnan(dphase_m))<4 %有效相位差数量小于4，不能定姿
                if sum(~isnan(dphase))>=5 %模糊度搜索
                    % 从A、dphase中取
                    index = find(~isnan(dphase)); %有相位差的索引号
                    Ac = A(index,:);
                    pc = dphase(index);
                    pc = mod(pc,1); %取小数部分
                    Rx = IAR(Ac, pc, lamda, bl+[-br,br], tr);
                    dphase_mask(isnan(dphase)) = NaN;
                    dphase_mask(~isnan(dphase)) = 0;
                    dN = round(dphase-An*Rx);
                    dN(isnan(dN)) = 0; %无有效相位差的通道不修正整周模糊度
                else %无法进行模糊度搜索
                    Rx = NaN(4,1);
                    dphase_mask(isnan(dphase)) = NaN;
                    dN = zeros(svN,1);
                end
            else %有效相位差数量大于等于4，直接定姿
                % 从An、dphase_m中取
                index = find(~isnan(dphase_m));
                Ac = An(index,:);
                pc = dphase_m(index);
                W = diag(Ac(:,3).^3); %高度角越高权值越大
                Rx = (Ac'*W*Ac) \ (Ac'*W*pc); %加权最小二乘
                dphase_mask(isnan(dphase)) = NaN;
                dphase_mask(~isnan(dphase)) = 0;
                dN = round(dphase-An*Rx);
                dN(isnan(dN)) = 0;
            end
            % 修正存在整周误差的通道
            ki = find(dN~=0)'; %索引，行向量
            for k=ki
                channels_A(k).carrCirc = mod(channels_A(k).carrCirc-dN(k), circ_limit);
            end
        else %深组合模式不进行模糊度搜索，不修正通道整周误差，在导航滤波器程序段修正
            if sum(~isnan(dphase_m))<4
                Rx = NaN(4,1);
            else
                An = [A/lamda, ones(svN,1)]; %A在导航滤波器里有算
                index = find(~isnan(dphase_m));
                Ac = An(index,:);
                pc = dphase_m(index);
%                 W = diag(Ac(:,3).^3); %高度角越高权值越大
                W = diag(1./sigma_dphase(index)')^2;
                Rx = (Ac'*W*Ac) \ (Ac'*W*pc); %加权最小二乘
            end
        end
        bl_length = norm(Rx(1:3));
        psi = atan2d(Rx(2),Rx(1));
        theta = -asind(Rx(3)/bl_length);
        
        % 5.时钟反馈修正
        if receiverState==1
            if ~isnan(pos(7))
                deltaFreq = deltaFreq + 10*pos(8)*dtpos/1000; %钟频差累加
                ta = ta - sec2smu(10*pos(7)*dtpos/1000); %时钟修正（可以不用进位，在下次更新时进位）
            end
            if ~isnan(Rx(4))
                deltaPath = deltaPath + 10*Rx(4)*dtpos/1000; %路径差累加
            end
        elseif receiverState==2
            deltaFreq = deltaFreq + NF.dtv;
            ta = ta - sec2smu(NF.dtr);
        end
        
        % 6.存储输出
        output_ta(no,1)     = tp(1) + tp(2)/1e3 + tp(3)/1e6;
        output_ta(no,2)     = receiverState;
        output_sv_A(:,:,no) = sv_A;
        output_sv_B(:,:,no) = sv_B;
        output_df(no)       = deltaFreq;
        output_dp(no)       = deltaPath;
        output_dphase(no,:) = dphase_m';
        output_pos(no,:)    = pos;
        output_Rx(no,:)     = [psi, theta, bl_length, Rx(4)];
        
        % 7.检查初始化
        if receiverState==0 && ~isnan(pos(7))
            if abs(pos(7))>0.1e-3 %钟差大于0.1ms，修正接收机时间
                ta = ta - sec2smu(pos(7)); %时钟修正
                ta = time_carry(ta);
                imu_index = find(imu_data(:,1)>(ta(1)+ta(2)/1e3+ta(3)/1e6), 1); %更新下次定位时间
            else %钟差小于0.1ms，初始化结束
                receiverState = 1;
            end
        end
        
        % 8.初始化导航滤波器
        if receiverState==1
            if abs(pos(8))<1e-10 && abs(Rx(4))<0.005 %钟频差收敛、路径差收敛
                %----初始化滤波器
                NF = navFilter_deep_dphase(pos(1:3), [0,0,0], [psi,theta,0], dtpos/1000, lamda, bl, para);
                %----卫星运动带来的伪距率，用来算载波加速度
                [~, rhodot0, ~, ~, ~] = cal_rho_rhodot(rs, vs, pos(1:3), [0,0,0]); %载体速度为0
                %----计算各通道的载波频率，确定载波环积分器上下界
                [~, rhodot , ~, ~, ~] = cal_rho_rhodot(rs, vs, pos(1:3), pos(4:6));
                carrFreq0 = -(rhodot/299792458 + deltaFreq) * 1575.42e6; %接收机钟频快，测的载波频率偏小
                %----更新通道跟踪模式
                for k=1:svN
                    if channels_A(k).state==2
                        channels_A(k).trackStage = 'D';
                        channels_A(k).PLL.upper = carrFreq0(k) + 1;
                        channels_A(k).PLL.lower = carrFreq0(k) - 1;
                    end
                    if channels_B(k).state==2
                        channels_B(k).trackStage = 'D';
                        channels_B(k).PLL.upper = carrFreq0(k) + 1;
                        channels_B(k).PLL.lower = carrFreq0(k) - 1;
                    end
                end
                %----更新接收机状态
                receiverState = 2;
            end
        end
        
        % 9.更新下次定位时间
        imu_index = imu_index + 1; %imu数据索引加1，指向下一个
        tp = sec2smu(imu_data(imu_index,1));
        no = no + 1; %指向下一存储位置
        
    end %end if dtp>=0
end

%% 关闭文件，关闭进度条
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
output_dp(no:end,:)         = [];
output_dphase(no:end,:)     = [];
output_Rx(no:end,:)         = [];
output_filter(no:end,:)     = [];
output_imu(no:end,:)        = [];
output_bias(no:end,:)       = [];
output_P(no:end,:)          = [];
output_svn(no:end,:)        = [];
% 删除接收机未初始化时的数据
index = find(output_ta(:,2)==0);
output_ta(index,:)          = [];
output_pos(index,:)         = [];
output_sv_A(:,:,index)      = [];
output_sv_B(:,:,index)      = [];
output_df(index,:)          = [];
output_dp(index,:)          = [];
output_dphase(index,:)      = [];
output_Rx(index,:)          = [];
output_filter(index,:)      = [];
output_imu(index,:)         = [];
output_bias(index,:)        = [];
output_P(index,:)           = [];
output_svn(index,:)         = [];

%% 打印通道日志
clc
disp('<--------antenna A-------->')
print_log([result_path,'\logA.txt'], svList);
disp('<--------antenna B-------->')
print_log([result_path,'\logB.txt'], svList);

%% 保存星历
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

%% 画跟踪结果
plot_track_double(sampleFreq, msToProcess, svList, trackResults_A, trackResults_B);

%% 画相位差
plot_svdata(output_dphase, svList, '相位差');

%% 清除变量
keepVariables = { ...
'sampleFreq'; 'msToProcess';
'p0'; 'tf'; 'svList'; 'svN'; %为了画星座图
'channels_A'; 'trackResults_A'; %A天线跟踪信息
'channels_B'; 'trackResults_B'; %B天线跟踪信息
'output_sv_A'; 'output_sv_B'; %卫星信息
'ephemeris'; 'ion'; %星历
'imu_data'; %IMU数据
'output_ta'; 'output_pos'; 'output_df'; 'output_dp'; 'output_dphase'; 'output_Rx';
'output_filter'; 'output_bias'; 'output_imu'; 'output_P'; 'output_svn';
'file';
};
clearvars('-except', keepVariables{:})

%% 保存结果
% 以时间命名
t0 = clock;
time_str = sprintf('%4d%02d%02d_%02d%02d%02d', t0(1),t0(2),t0(3),t0(4),t0(5),floor(t0(6)));
result_path = fileread('.\temp\resultPath.txt');
save([result_path,'\',time_str,'__deeply_double__',file(1:end-8),'.mat'])

%% 计时结束
toc

%% 函数
function ch = GPS_L1_CA_channel_struct()
% 声明通道结构体所有场
    ch.PRN              = []; %卫星编号
    ch.state            = []; %通道状态（数字）
    ch.trackStage       = []; %跟踪阶段（字符）
    ch.msgStage         = []; %电文解析阶段（字符）
    ch.strength         = []; %信号强度（数字）
    ch.cnt_t            = []; %跟踪时用的计数器
    ch.cnt_m            = []; %电文解析时用的计数器
    ch.stableCnt        = []; %信号稳定计数器
    ch.code             = []; %伪码
    ch.timeIntMs        = []; %积分时间，ms
    ch.trackDataTail    = []; %跟踪开始点在数据缓存中的位置
    ch.blkSize          = []; %跟踪数据段采样点个数
    ch.trackDataHead    = []; %跟踪结束点在数据缓存中的位置
    ch.dataIndex        = []; %跟踪开始点在文件中的位置
    ch.ts0              = []; %跟踪开始点所在码周期的理论发射时间，ms
    ch.carrNco          = []; %载波发生器频率
    ch.codeNco          = []; %码发生器频率
    ch.carrAcc          = []; %载波频率加速度
    ch.carrFreq         = []; %载波频率测量
    ch.codeFreq         = []; %码频率测量
    ch.remCarrPhase     = []; %跟踪开始点的载波相位
    ch.remCodePhase     = []; %跟踪开始点的码相位
    ch.carrCirc         = []; %记录载波经过的整周数，0~999
    ch.I_P0             = []; %上次跟踪的I_P
    ch.Q_P0             = []; %上次跟踪的Q_P
    ch.FLL              = []; %锁频环（结构体）
    ch.PLL              = []; %锁相环（结构体）
    ch.DLL              = []; %延迟锁定环（结构体）
    ch.bitSyncTable     = []; %比特同步统计表
    ch.bitBuff          = []; %比特缓存
    ch.frameBuff        = []; %帧缓存
    ch.frameBuffPoint   = []; %帧缓存指针
    ch.ephemeris        = []; %星历
    ch.codeStd          = []; %计算码鉴相器误差标准差结构体
    ch.carrStd          = []; %计算载波鉴相器误差标准差结构体
    ch.NWmean           = []; %计算NBP/WBP均值结构体
    ch.CN0              = []; %平均载噪比
    ch.CN0i             = []; %瞬时载噪比
end

function ch = GPS_L1_CA_channel_init(ch, acqResult, n, sampleFreq)
% 通道结构体初始化，执行完后通道被激活
% n表示已经经过了多少个采样点
% 需要预先赋PRN
    code = GPS_L1_CA_generate(ch.PRN); %C/A码

    % ch.PRN 卫星号不变
    ch.state = 1; %激活通道
    ch.trackStage = 'F'; %频率牵引
    ch.msgStage = 'I'; %空闲
    ch.strength = 0; %信号失锁
    ch.cnt_t = 0;
    ch.cnt_m = 0;
    ch.stableCnt = 0;
    ch.code = [code(end),code,code(1)]'; %列向量，为了求积分时用矢量相乘加速
    ch.timeIntMs = 1;
    ch.trackDataTail = sampleFreq*0.001 - acqResult(1) + 2;
    ch.blkSize = sampleFreq*0.001;
    ch.trackDataHead = ch.trackDataTail + ch.blkSize - 1;
    ch.dataIndex = ch.trackDataTail + n;
    ch.ts0 = NaN;
    ch.carrNco = acqResult(2);
    ch.codeNco = 1.023e6 + ch.carrNco/1540;
    ch.carrAcc = 0;
    ch.carrFreq = ch.carrNco;
    ch.codeFreq = ch.codeNco;
    ch.remCarrPhase = 0;
    ch.remCodePhase = 0;
    ch.carrCirc = 0;
    ch.I_P0 = NaN;
    ch.Q_P0 = NaN;

    ch.FLL.K = 40;
    ch.FLL.Int = ch.carrNco;

    [K1, K2] = orderTwoLoopCoef(25, 0.707, 1);
    ch.PLL.K1 = K1;
    ch.PLL.K2 = K2;
    ch.PLL.Int = 0;
    ch.PLL.upper = 0; %积分器上界
    ch.PLL.lower = 0; %积分器下界

    [K1, K2] = orderTwoLoopCoef(2, 0.707, 1);
    ch.DLL.K1 = K1;
    ch.DLL.K2 = K2;
    ch.DLL.Int = ch.codeNco;

    ch.bitSyncTable = zeros(1,20);
    ch.bitBuff = zeros(2,20); %第一行I_P，第二行Q_P
    ch.frameBuff = zeros(1,1502);
    ch.frameBuffPoint = 0;
    % ch.ephemeris 星历不变

    % 计算码鉴相器误差标准差结构体
    ch.codeStd.buff = zeros(1,200);
    ch.codeStd.buffSize = length(ch.codeStd.buff);
    ch.codeStd.buffPoint = 0;
    ch.codeStd.E0 = 0;
    ch.codeStd.D0 = 0;

    % 计算载波鉴相器误差标准差结构体
    ch.carrStd.buff = zeros(1,200);
    ch.carrStd.buffSize = length(ch.carrStd.buff);
    ch.carrStd.buffPoint = 0;
    ch.carrStd.E0 = 0;
    ch.carrStd.D0 = 0;

    % 计算NBP/WBP均值结构体
    ch.NWmean.buff = zeros(1,50); %50个点求均值
    ch.NWmean.buffSize = length(ch.NWmean.buff);
    ch.NWmean.buffPoint = 0;
    ch.NWmean.E0 = 0;
    ch.CN0 = 0;
    ch.CN0i = 0;
end

function trackResult = trackResult_struct(m)
% 跟踪结果结构体
    trackResult.PRN = 0;
    trackResult.n = 1; %指向当前存储的行号
    trackResult.dataIndex     = zeros(m,1); %码周期开始采样点在原始数据文件中的位置
    trackResult.ts0           = zeros(m,1); %码周期理论发射时间，ms
    trackResult.remCodePhase  = zeros(m,1); %码周期开始采样点的码相位，码片
    trackResult.codeFreq      = zeros(m,1); %码频率
    trackResult.remCarrPhase  = zeros(m,1); %码周期开始采样点的载波相位，周
    trackResult.carrFreq      = zeros(m,1); %载波频率
    trackResult.I_Q           = zeros(m,6); %[I_P,I_E,I_L,Q_P,Q_E,Q_L]
    trackResult.disc          = zeros(m,5); %[codeError,std, carrError,std, freqError]，鉴相器
    trackResult.bitStartFlag  = zeros(m,1); %比特开始标志
    trackResult.CN0           = zeros(m,2); %载噪比（平均和瞬时）
    trackResult.carrAcc       = zeros(m,1); %载波加速度
    trackResult.strength      = zeros(m,1); %信号强度
end

function trackResult = trackResult_clean(trackResult)
% 清理跟踪结果中的空白空间
    n = trackResult.n;
    trackResult.dataIndex(n:end,:)    = [];
    trackResult.ts0(n:end,:)          = [];
    trackResult.remCodePhase(n:end,:) = [];
    trackResult.codeFreq(n:end,:)     = [];
    trackResult.remCarrPhase(n:end,:) = [];
    trackResult.carrFreq(n:end,:)     = [];
    trackResult.I_Q(n:end,:)          = [];
    trackResult.disc(n:end,:)         = [];
    trackResult.bitStartFlag(n:end,:) = [];
    trackResult.CN0(n:end,:)          = [];
    trackResult.carrAcc(n:end,:)      = [];
    trackResult.strength(n:end,:)     = [];
end
