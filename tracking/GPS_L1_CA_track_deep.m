function [ch, I_Q, disc, bitStartFlag] = GPS_L1_CA_track_deep(ch, sampleFreq, buffSize, rawSignal, logID)

bitStartFlag = 0;

%% 提取通道信息（跟踪算法用到的控制参数）
trackStage     = ch.trackStage;
msgStage       = ch.msgStage;
cnt_t          = ch.cnt_t;
cnt_m          = ch.cnt_m;
code           = ch.code;
timeIntMs      = ch.timeIntMs;
blkSize        = ch.blkSize;
carrNco        = ch.carrNco;
codeNco        = ch.codeNco;
remCarrPhase   = ch.remCarrPhase;
remCodePhase   = ch.remCodePhase;
carrCirc       = ch.carrCirc;
I_P0           = ch.I_P0;
Q_P0           = ch.Q_P0;
FLL            = ch.FLL;
PLL            = ch.PLL;
DLL            = ch.DLL;
bitSyncTable   = ch.bitSyncTable;
bitBuff        = ch.bitBuff;
frameBuff      = ch.frameBuff;
frameBuffPoint = ch.frameBuffPoint;

ch.dataIndex = ch.dataIndex + blkSize; %下个数据段开始的数据索引
ch.ts0       = ch.ts0 + timeIntMs; %当前码周期结束的时间，ms

timeInt = timeIntMs * 0.001; %积分时间，s
codeInt = timeIntMs * 1023; %积分时间内码片个数
pointInt = 20 / timeIntMs; %一个电文比特中有多少个积分点

%% 基本处理
% 时间序列
t = (0:blkSize-1) / sampleFreq;
te = blkSize / sampleFreq;

% 生成本地载波
theta = (remCarrPhase + carrNco*t) * 2; %变频率载波，乘2因为后面是以pi为单位求三角函数
carr_cos = cospi(theta); %本地载波
carr_sin = sinpi(theta);
theta_next = remCarrPhase + carrNco*te;
remCarrPhase = mod(theta_next, 1); %剩余载波相位，周
carrCirc = mod(floor(carrCirc+theta_next), 1000); %载波经过的整周数

% 生成本地码
tcode = remCodePhase + codeNco*t + 2; %加2保证求滞后码时大于1
earlyCode  = code(floor(tcode+0.5)); %超前码
promptCode = code(floor(tcode));     %即时码
lateCode   = code(floor(tcode-0.5)); %滞后码
remCodePhase = remCodePhase + codeNco*te - codeInt; %剩余载波相位，周

% 原始数据乘载波
% 将复数相乘拆开，为了避免取实部、虚部时调用函数耗时
iBasebandSignal = rawSignal(1,:).*carr_cos + rawSignal(2,:).*carr_sin; %乘负载波
qBasebandSignal = rawSignal(2,:).*carr_cos - rawSignal(1,:).*carr_sin;

% 六路积分
% 用行向量与列向量相乘代替sum(X.*X)，加速
I_E = iBasebandSignal * earlyCode;
Q_E = qBasebandSignal * earlyCode;
I_P = iBasebandSignal * promptCode;
Q_P = qBasebandSignal * promptCode;
I_L = iBasebandSignal * lateCode;
Q_L = qBasebandSignal * lateCode;

% 码鉴相器
S_E = sqrt(I_E^2+Q_E^2);
S_L = sqrt(I_L^2+Q_L^2);
codeError = 0.5 * (S_E-S_L)/(S_E+S_L); %单位：码片
[ch.codeStd, codeSigma] = std_rec(ch.codeStd ,codeError); %计算码鉴相器误差标准差

% 载波鉴相器
carrError = atan(Q_P/I_P) / (2*pi); %单位：周
[ch.carrStd, carrSigma] = std_rec(ch.carrStd ,carrError); %计算载波鉴相器误差标准差

% 鉴频器
if ~isnan(I_P0)
    yc = I_P0*I_P + Q_P0*Q_P;
    ys = I_P0*Q_P - Q_P0*I_P;
    freqError = atan(ys/yc)/timeInt / (2*pi); %单位：Hz
else
    freqError = 0;
end

%% 跟踪算法
switch trackStage
    case 'F' %<<====频率牵引
        %----FLL
        FLL.Int = FLL.Int + FLL.K*freqError*timeInt; %锁频环积分器
        carrNco = FLL.Int;
        carrFreq = FLL.Int;
        % 500ms后转到传统跟踪
        cnt_t = cnt_t + 1;
        if cnt_t==500
            cnt_t = 0; %计数器清零
            PLL.Int = FLL.Int; %初始化锁相环积分器
            trackStage = 'T';
            fprintf(logID, '%2d: Start traditional tracking at %.8fs\r\n', ...
                    ch.PRN, ch.dataIndex/sampleFreq);
        end
        %----DLL
        DLL.Int = DLL.Int + DLL.K2*codeError*timeInt; %延迟锁定环积分器
        codeNco = DLL.Int + DLL.K1*codeError;
        codeFreq = DLL.Int;
        
	case 'T' %<<====传统跟踪
        %----PLL
        PLL.Int = PLL.Int + PLL.K2*carrError*timeInt; %锁相环积分器
        carrNco = PLL.Int + PLL.K1*carrError;
        carrFreq = PLL.Int;
        % 500ms后进行比特同步
        if msgStage=='I'
            cnt_t = cnt_t + 1;
            if cnt_t==500
                cnt_t = 0; %计数器清零
                msgStage = 'B';
                fprintf(logID, '%2d: Start bit synchronization at %.8fs\r\n', ...
                        ch.PRN, ch.dataIndex/sampleFreq);
            end
        end
        %----DLL
        DLL.Int = DLL.Int + DLL.K2*codeError*timeInt; %延迟锁定环积分器
        codeNco = DLL.Int + DLL.K1*codeError;
        codeFreq = DLL.Int;
        
	case 'D' %<<====深组合码环和锁相环都开环
        %----PLL
        if ch.strength==2 %强信号时使用二阶环（比例积分控制）得到载波频率测量
            PLL.Int = PLL.Int + PLL.K2*carrError*timeInt + ch.carrAcc*timeInt; %锁相环积分器
            if PLL.Int>PLL.upper %积分饱和
                PLL.Int = PLL.upper;
            elseif PLL.Int<PLL.lower
                PLL.Int = PLL.lower;
            end
            carrNco = PLL.Int + PLL.K1*carrError;
            carrFreq = PLL.Int;
        else %弱信号时使用频率辅助的一阶环（比例控制）跟踪载波相位
            carrNco = PLL.Int + PLL.K1*carrError;
%             carrNco = PLL.Int + 8*carrError;
            carrFreq = PLL.Int;
        end
        %----DLL
        codeNco = 1.023e6 + carrFreq/1540; %码频率直接由载波频率计算
        codeFreq = 1.023e6 + carrFreq/1540;
        
	otherwise
end

%% 电文解析算法
switch msgStage %I, B, W, H, C, E
    case 'I' %<<====空闲
        
    case 'B' %<<====比特同步
        % 必须是1ms积分时间，进行2s，有100个比特
        % 比特同步后可以实现更长的积分时间
        % 比特同步后可以进行载噪比计算
        cnt_m = cnt_m + 1;
        if (I_P0*I_P)<0 %发现电平翻转
            index = mod(cnt_m-1,20) + 1;
            bitSyncTable(index) = bitSyncTable(index) + 1; %统计表中的对应位加1
        end
        if cnt_m==2000 %2s后检验统计表
            if max(bitSyncTable)>10 && (sum(bitSyncTable)-max(bitSyncTable))<=2 %确定电平翻转位置（电平翻转大都发生在一个点上）
                [~,cnt_m] = max(bitSyncTable);
                bitSyncTable = zeros(1,20); %比特同步统计表清零
                cnt_m = -cnt_m + 1;
                if cnt_m==0
                    msgStage = 'H'; %进入寻找帧头模式
                    fprintf(logID, '%2d: Start find head at %.8fs\r\n', ...
                            ch.PRN, ch.dataIndex/sampleFreq);
                else
                    msgStage = 'W'; %进入等待cnt_m==0模式
                end
            else
                ch.state = 0; %比特同步失败，关闭通道
                fprintf(logID, '%2d: ***Bit synchronization fails at %.8fs\r\n', ...
                        ch.PRN, ch.dataIndex/sampleFreq);
            end
        end
        
    case 'W' %<<====等待cnt_m==0
        cnt_m = cnt_m + 1;
        if cnt_m==0
            msgStage = 'H'; %进入寻找帧头模式
            fprintf(logID, '%2d: Start find head at %.8fs\r\n', ...
                    ch.PRN, ch.dataIndex/sampleFreq);
        end
        
    otherwise %<<====已经完成比特同步
        cnt_m = cnt_m + 1;
        bitBuff(1,cnt_m) = I_P; %往比特缓存中存数
        bitBuff(2,cnt_m) = Q_P; %往比特缓存中存数
        if cnt_m==1 %标记当前跟踪的数据段为比特开始位置
            bitStartFlag = double(msgStage);
        end
        if cnt_m==pointInt %跟踪完一个比特
            cnt_m = 0; %计数器清零
            %-------------------------------------------------------------%
            %====计算平均载噪比
            Ps = bitBuff(1,1:pointInt).^2 + bitBuff(2,1:pointInt).^2; %每个点的功率
            WBP = sum(Ps); %宽带功率，所有点的功率求和（先平方再求和）
            Is = sum(bitBuff(1,1:pointInt)); %合成I
            Qs = sum(bitBuff(2,1:pointInt)); %合成Q
            NBP = Is^2 + Qs^2; %窄带功率，合成IQ的功率，信号越好，窄带功率越大（先求和再平方）
            if ch.CN0==0 %初始化求均值结构体（通道刚激活时CN0为0）
                ch.NWmean.buff = ones(1,ch.NWmean.buffSize)*(NBP/WBP);
                ch.NWmean.E0 = NBP/WBP;
            end
            [ch.NWmean, NWm] = mean_rec(ch.NWmean, NBP/WBP); %计算Z的均值
            S = (NWm-1) / (pointInt-NWm) / timeInt;
            if S>10
                CN0 = 10*log10(S); %载噪比
            else
                CN0 = 10; %计算载噪比的最小值为10，刚初始化时为0
            end
            ch.CN0 = CN0;
            %====没进深组合时失锁要放弃通道
            if CN0<30
                if trackStage=='T'
                    ch.state = 0;
                    fprintf(logID, '%2d: ***Abandon at %.8fs\r\n', ...
                            ch.PRN, ch.dataIndex/sampleFreq);
                end
            end
            %====计算瞬时载噪比
            Z = NBP / WBP;
            S = (Z-1) / (pointInt-Z) / timeInt;
            if S>10
                CN0i = 10*log10(S); %瞬时载噪比
            else
                CN0i = 10;
            end
            ch.CN0i = CN0i;
            %====信号穿越判断
            if length(unique(sign(bitBuff(1,1:pointInt))))==1 %一个比特内的所有点符号都相同
                through_flag = 0; %无穿越
            else
                through_flag = 1; %有穿越
            end
            %====更新信号稳定计数器
            stableCnt = ch.stableCnt;
            if CN0i>=30 && through_flag==0 %载噪比大于阈值并且无穿越时，计数器加1，计数器最大值为50（1s）
                if (stableCnt+1)>50
                    stableCnt = 50; %表示信号已经稳定了1s
                else
                    stableCnt = stableCnt + 1;
                end
            else
                stableCnt = 0; %检测到载噪比小于阈值或者在一个比特内存在积分点穿越，计数器清零
            end
            ch.stableCnt = stableCnt;
            %====判断信号强度
            strength = ch.strength;
            if strength==2 %强信号
                if stableCnt==0
                    ch.strength = 1;
                    fprintf(logID, '%2d: ***Weak signal at %.8fs\r\n', ...
                            ch.PRN, ch.dataIndex/sampleFreq);
                end
            elseif strength==1 %弱信号
                if CN0<30
                    ch.strength = 0;
                    fprintf(logID, '%2d: ***Lose lock at %.8fs\r\n', ...
                            ch.PRN, ch.dataIndex/sampleFreq);
                end
                if stableCnt==50
                    ch.strength = 2;
                    fprintf(logID, '%2d: Strong signal at %.8fs\r\n', ...
                            ch.PRN, ch.dataIndex/sampleFreq);
                end
            elseif strength==0 %失锁
                if CN0>=32
                    ch.strength = 1;
                    fprintf(logID, '%2d: Weak signal at %.8fs\r\n', ...
                            ch.PRN, ch.dataIndex/sampleFreq);
                end
            end
            %-------------------------------------------------------------%
            bit = sum(bitBuff(1,1:pointInt)) > 0; %判断比特值，0/1
            frameBuffPoint = frameBuffPoint + 1;
            frameBuff(frameBuffPoint) = (double(bit) - 0.5) * 2; %存储比特值，±1
            switch msgStage
                case 'H' %<<====寻找帧头
                    if frameBuffPoint>=10 %至少有10个比特，前两个用来校验
                        if abs(sum(frameBuff(frameBuffPoint+(-7:0)).*[1,-1,-1,-1,1,-1,1,1]))==8 %检测到疑似帧头
                            frameBuff(1:10) = frameBuff(frameBuffPoint+(-9:0)); %将帧头提前
                            frameBuffPoint = 10;
                            msgStage = 'C'; %进入校验帧头模式
                        end
                    end
                    if frameBuffPoint==1502 %防止Bug，一般到不了这里，30s还没找到帧头早就被判定为失锁了
                        frameBuffPoint = 0;
                    end
                case 'C' %<<====校验帧头
                    if frameBuffPoint==310 %存储了一个子帧，2+300+8
                        if GPS_L1_CA_check(frameBuff(1:32))==1 && GPS_L1_CA_check(frameBuff(31:62))==1 && ... %校验通过
                            abs(sum(frameBuff(303:310).*[1,-1,-1,-1,1,-1,1,1]))==8
                            % 获取电文时间
                            % frameBuff(32)为上一字的最后一位，校验时控制电平翻转，为1表示翻转，为0表示不翻转，参见ICD-GPS最后几页
                            bits = -frameBuff(32) * frameBuff(33:49); %电平翻转，31~47比特
                            bits = dec2bin(bits>0)'; %±1数组转化为01字符串
                            TOW = bin2dec(bits); %01字符串转换为十进制数
                            ch.ts0 = (TOW*6+0.16)*1000; %ms，0.16=8/50
                            if ~isnan(ch.ephemeris(1))
                                ch.state = 2; %更新状态（知道码发射时间，而且有星历）
                            end
                            msgStage = 'E'; %进入解析星历模式
                            fprintf(logID, '%2d: Start parse ephemeris at %.8fs\r\n', ...
                                    ch.PRN, ch.dataIndex/sampleFreq);
                        else %校验未通过
                            for k=11:310 %检查其他比特中有没有帧头
                                if abs(sum(frameBuff(k+(-7:0)).*[1,-1,-1,-1,1,-1,1,1]))==8 %检测到疑似帧头
                                    frameBuff(1:320-k) = frameBuff(k-9:310); %将帧头和后面的比特提前，320-k = 310-(k-9)+1
                                    frameBuffPoint = 320-k;
                                    break
                                end
                            end
                            if frameBuffPoint==310 %没检测到疑似帧头
                                frameBuff(1:9) = frameBuff(302:310); %将未检测的比特提前
                                frameBuffPoint = 9;
                                msgStage = 'H'; %再次寻找帧头
                            end
                        end
                    end
                case 'E' %<<====解析星历
                    if frameBuffPoint==1502 %跟踪完5帧
                        ephemeris = GPS_L1_CA_ephemeris(frameBuff); %解析星历
                        if isempty(ephemeris) %星历错误
                            fprintf(logID, '%2d: ***Ephemeris error at %.8fs\r\n', ...
                                    ch.PRN, ch.dataIndex/sampleFreq);
                            %-------------------------------------------------------------%
%                             bits = -frameBuff(62) * frameBuff; %电平翻转
%                             bits = dec2bin(bits>0)'; %±1数组转化为01字符串
%                             fprintf(logID, ['%2d: ',bits(1:2),'\r\n'], ch.PRN);
%                             for k=1:50 %将错误的电文输出，查找电文错误原因
%                                 fprintf(logID, ['%2d: ',bits((k-1)*30+2+(1:30)),'\r\n'], ch.PRN);
%                             end
                            %-------------------------------------------------------------%
                            frameBuffPoint = 0;
                            msgStage = 'H'; %重新寻找帧头
                            fprintf(logID, '%2d: Start find head at %.8fs\r\n', ...
                                    ch.PRN, ch.dataIndex/sampleFreq);
                        else
                            if ephemeris(2)~=ephemeris(3) %星历改变
                                fprintf(logID, '%2d: ***Ephemeris changes at %.8fs, IODC=%d, IODE=%d\r\n', ...
                                        ch.PRN, ch.dataIndex/sampleFreq, ephemeris(2), ephemeris(3));
                            else
                                ch.ephemeris = ephemeris; %更新星历
                                ch.state = 2; %更新状态（知道码发射时间，而且有星历）
                                fprintf(logID, '%2d: Ephemeris is parsed at %.8fs\r\n', ...
                                        ch.PRN, ch.dataIndex/sampleFreq);
                            end
                            frameBuff(1:2) = frameBuff(1501:1502); %将最后两个比特提前
                            frameBuffPoint = 2;
                        end
                    end
                otherwise
            end
        end
        
end

%% 更新下一数据块位置
trackDataTail = ch.trackDataHead + 1;
if trackDataTail>buffSize
    trackDataTail = 1;
end
blkSize = ceil((codeInt-remCodePhase)/codeNco*sampleFreq);
trackDataHead = trackDataTail + blkSize - 1;
if trackDataHead>buffSize
    trackDataHead = trackDataHead - buffSize;
end
ch.trackDataTail = trackDataTail;
ch.blkSize       = blkSize;
ch.trackDataHead = trackDataHead;

%% 更新通道信息
ch.trackStage     = trackStage;
ch.msgStage       = msgStage;
ch.cnt_t          = cnt_t;
ch.cnt_m          = cnt_m;
ch.code           = code;
ch.timeIntMs      = timeIntMs;
ch.carrNco        = carrNco;
ch.codeNco        = codeNco;
ch.carrFreq       = carrFreq;
ch.codeFreq       = codeFreq;
ch.remCarrPhase   = remCarrPhase;
ch.remCodePhase   = remCodePhase;
ch.carrCirc       = carrCirc;
ch.I_P0           = I_P;
ch.Q_P0           = Q_P;
ch.FLL            = FLL;
ch.PLL            = PLL;
ch.DLL            = DLL;
ch.bitSyncTable   = bitSyncTable;
ch.bitBuff        = bitBuff;
ch.frameBuff      = frameBuff;
ch.frameBuffPoint = frameBuffPoint;

%% 输出
I_Q = [I_P, I_E, I_L, Q_P, Q_E, Q_L];
disc = [codeError, codeSigma, carrError, carrSigma, freqError];

end