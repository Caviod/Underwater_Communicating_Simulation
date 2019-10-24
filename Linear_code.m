clc;clear all;close all;
%%%%%%%%生成初始序列%%%%%%%
SignNum = 100; %信息长度
genmat=[1 0 0 0 1 1 1;0 1 0 0 1 1 0;0 0 1 0 1 0 1;0 0 0 1 0 1 1];
Length = SignNum ; %原始码序列长度
OrigiSeq=(sign(randn(1,SignNum))+1)/2; %产生原始的0或1信息序列
%%%%%%%%hamming编码%%%%%
ConCode=encode(OrigiSeq,7,4,'linear',genmat); %线性分组码编码后序列
%%%%%BPSK调制%%%%%%%%
BPSKCode =pskmod(ConCode,2);
%%%%%%%定义信噪比%%%%%%%%%%%%
EbN0 =-6:6; %EbN0-dB
snr=10.^(EbN0./10); %转化为功率信噪比


%%%%%%%%%%%%%信噪比循环接收统计误码率%%%%%%%%%%%%
for k = 1:length(EbN0)
    for b=1:100;
        %%%%%%%%%%%%%%加性高斯白噪声信道%%%%%%%%%%%
        RecCode = awgn(BPSKCode,EbN0(k),'measured');
        %%%%%%%%%%%%%%BPSK解调%%%%%%%%%%%%%%%
        
        % BPSKdecode = BPSKDecode(RecCode, Length*7/4); ? %BPSK解调
        
        BPSKdecode=pskdemod(RecCode,2); %%%%%%%上面的式子可以直接用自带pskdemod函数来解调，而不用调用函数%%%%%%%%%
        
        %% 纠错部分
        
        OrigiSeq2 =reshape(OrigiSeq,(SignNum)/4,4); %%原始bit信息重装成 25*4矩阵
        SignNums=length(BPSKdecode);
        BPSKdecode1=reshape((BPSKdecode)',7,(SignNums)/7); %%重装成 25*7矩阵
        BPSKdecode1=(BPSKdecode1)';
        
        
        BPSKdecode=jiucuo(OrigiSeq2,BPSKdecode1); %%%%%调用纠错
        
        %%%% ?a:注意reshape函数的使用，按列读取，按列存储
        %%%% ?b:BPSKdecode=jiucuo(OrigiSeq2,BPSKdecode1)中其实只需一个变量BPSKdecode1即可
        
        
        
        [g,h]=size(BPSKdecode); % g=1,h=175
        %%%%%%%%%%%%%% ? ? ? ?linear码译码%%%%%%%%%%%%%%
        Decoder=decode(BPSKdecode,7,4,'linear',genmat);
        [m,n]=size(Decoder); % m=1,n=100
        
        
        %%%%%%%统计误码率%%%%%%%%
        
        [o,p]=size(OrigiSeq); %%% ?o=1,p=100
        error(b) = sum(abs(Decoder-OrigiSeq))/Length;
    end
    %%%%%%%%%%%%%%计算平均误码率%%%%%%%%%%%%%
    errorout(k)=mean(error);
    OrigiSeq1=OrigiSeq(1:20);%提取前20个仿真结果
    Decoder1=Decoder(1:20);
%     figure
    t=1:20;
    subplot(211);stairs(t,OrigiSeq1,'r');axis([1 20 -0.5 1.5]);title(['编码前序列',num2str(EbN0(k)),'dB']);
    subplot(212);stairs(t,Decoder1,'b');axis([1 20 -0.5 1.5]);title(['译码后序列',num2str(EbN0(k)),'dB']);
    F = getframe;
    movie(F,1,0.5)
end
close all;


%%%%%%%%%%%%%计算理论误码率%%%%%%%%%%%%%
ber_theory=0.5*erfc(sqrt(snr)); %计算理论误比特率


ber_theory;

%%%%%%%%%%%%%绘制误码率曲线%%%%%%%%%%%%%%
figure; 
semilogy (EbN0,errorout,'b*-',EbN0,ber_theory,'rd-');
hold on
xlabel('Eb/N0(dB)');
ylabel('误比特率');
legend('加 linear code','理论值');
title('BPSK+linear误比特率曲线');
grid on; %%%%添加网格线






