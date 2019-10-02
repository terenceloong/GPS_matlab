for k=1:svN
    %====天线A
    t = trackResults_A(k).dataIndex/sampleFreq;
    I_Q = trackResults_A(k).I_Q(:,1);
    
    I_Q0 = I_Q;
    I_Q0(trackResults_A(k).strength~=0) = NaN;
    
    I_Q1 = I_Q;
    I_Q1(trackResults_A(k).strength~=1) = NaN;
    
    I_Q2 = I_Q;
    I_Q2(trackResults_A(k).strength~=2) = NaN;
    
    figure
    subplot(2,1,1)
    plot(t, I_Q2) %强信号
    hold on
    plot(t, I_Q1) %弱信号
    plot(t, I_Q0) %失锁
    set(gca, 'XLim',[0,msToProcess/1000])
    
    %====天线A
    t = trackResults_B(k).dataIndex/sampleFreq;
    I_Q = trackResults_B(k).I_Q(:,1);
    
    I_Q0 = I_Q;
    I_Q0(trackResults_B(k).strength~=0) = NaN;
    
    I_Q1 = I_Q;
    I_Q1(trackResults_B(k).strength~=1) = NaN;
    
    I_Q2 = I_Q;
    I_Q2(trackResults_B(k).strength~=2) = NaN;
    
    subplot(2,1,2)
    plot(t, I_Q2) %强信号
    hold on
    plot(t, I_Q1) %弱信号
    plot(t, I_Q0) %失锁
    set(gca, 'XLim',[0,msToProcess/1000])
end

clearvars t I_Q I_Q0 I_Q1 I_Q2