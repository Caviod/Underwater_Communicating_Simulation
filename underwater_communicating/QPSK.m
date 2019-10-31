clc;clear all;close all;
%假定接收端已经实现载波同步，位同步（盲信号解调重点要解决的问题：载波同步（costas环（未见到相关代码）），位同步（Gardner算法（未见相关代码）），帧同步）
% carrier frequency  for modulation and demodulation  
fc=5e6;
%QPSK transmitter
data=5000  ;   %码数率为5MHZ     %原码个数
rand_data=randn(1,5000);
for  i=1:data
    if rand_data(i)>=0.5
        rand_data(i)=1;
    else
        rand_data(i)=0;
    end
end
%seriel to parallel        %同时单极性码转为双极性码
for  i=1:data
    if rem(i,2)==1
        if  rand_data(i)==1
            I(i)=1;
            I(i+1)=1;
        else
            I(i)=-1;
            I(i+1)=-1;
        end
    else
        if rand_data(i)==1
            Q(i-1)=1;
            Q(i)=1;
        else
            Q(i-1)=-1;
            Q(i)=-1;
        end
    end
end
% zero insertion   ，此过程称为成形。成形的意思就是实现由消息到波形的转换，以便发射，脉冲成形应该是在基带调制之后。
zero=5;         %sampling  rate  25M HZ  ,明白了，zero为过采样率。它等于 采样率fs/码速率。
for  i=1:zero*data     % 采样点数目=过采样率*原码数目
    if rem(i,zero)==1
        Izero(i)=I(fix((i-1)/zero)+1);
        Qzero(i)=Q(fix((i-1)/zero)+1);
    else
        Izero(i)=0;
        Qzero(i)=0;
    end
end
%pulse shape filter， 接着，将进行低通滤波，因为 随着传输速率的增大，基带脉冲的频谱将变宽
%如果不滤波（如升余弦滤波）进行低通滤波，后面加载频的时候可能会出现困难。
%平方根升余弦滤波器
% psf=rcosfir(rf,n_t,rate,fs,'sqrt')   rate:过采样率，rf:滚降因子，n_t:滤波器阶数，fs:采样率
%用在调制或发送之前，用在解调或接受之后，用来降低过采样符号流带宽并不引发ISI（码间串扰）
  
NT=50;
N=2*zero*NT;    % =500
fs=25e6;
rf=0.1;
psf=rcosfir(rf,NT,zero,fs,'sqrt');% psf大小为500
Ipulse=conv(Izero,psf);
Qpulse=conv(Qzero,psf);
%为什么数字信号传输也要过采样，成形滤波？
%答：过采样的数字信号处理起来对低通滤波器的要求相对较低，如果不过采样，滤波的时候滤波器需要很陡峭，指标会很严格
%成形滤波的作用是保证采样点不失真。如果没有它，那信号在经过带限信道后，眼图张不开，ISI非常严重。成形滤波的位置在基带调制之后。
%因为经成形滤波后，信号的信息已经有所损失，这也是为避免ISI付出的代价。换句话说，成形滤波的位置在载波调制之前，仅挨着载波调制。
%即：（发送端）插值（采样）-成形-滤波（LPF)-加载频(载波调制)-加噪声至（接收端）乘本振-低通-定时抽取-判决。

%modulation
for i=1:zero*data+N   %采样点数目改变 （因为卷积的缘故）
    t(i)=(i-1)/(fs);  %这里因为假设载频与码速率大小相等，所以用载频fc乘以过采样率=采样率。
    Imod(i)=Ipulse(i)*sqrt(2)*cos(2*pi*fc*t(i));
    Qmod(i)=Qpulse(i)*(-sqrt(2)*sin(2*pi*fc*t(i)));
end
sum=Imod+Qmod;
%QPSK  receiver
%demodulation
   for i=1:zero*data+N
       Idem(i)=sum(i)*sqrt(2)*cos(2*pi*fc*t(i));
       Qdem(i)=sum(i)*(-sqrt(2)*sin(2*pi*fc*t(i)));
   end
   %matched  filter
   mtf=rcosfir(rf,NT,zero,fs,'sqrt');
   Imat=conv(Idem,mtf);
   Qmat=conv(Qdem,mtf);
   %data selection
   for  i=1:zero*data
       Isel(i)=Imat(i+N);
       Qsel(i)=Qmat(i+N);
   end
   %sampler        %提取码元  
   for i=1:data
       Isam(i)=Isel((i-1)*zero+1);
       Qsam(i)=Qsel((i-1)*zero+1);
   end
   %decision  threshold
   threshold=0.2;
   for  i=1:data
       if Isam(i)>=threshold
           Ifinal(i)=1;
       else
           Ifinal(i)=-1;
       end
       if Qsam(i)>=threshold
           Qfinal(i)=1;
       else
           Qfinal(i)=-1;
       end
   end
   %parallel to serial
   for i=1:data
       if rem (i,2)==1
           if Ifinal(i)==1
               final(i)=1;
           else
               final(i)=0;
           end
       else
           if  Qfinal(i)==1
               final(i)=1;
           else
               final(i)=0;
           end
       end
   end
   % 绘图
   figure(1)
   plot(20*log(abs(fft(rand_data))));
   axis([0  data  -40  100]);
   grid on;
   title('spectrum  of input binary data');
   figure(2)
   subplot(221);
   plot(20*log(abs(fft(I))));
   axis([0 data -40 140]);
   grid  on;
   title('spectrum of I-channel data');
   subplot(222);
   plot(20*log(abs(fft(Q))));
   axis([0  data   -40  140]);
   grid  on;
   title('spectrum of Q-channel data');
   subplot(223);
   plot(20*log(abs(fft(Izero))));
   axis([0 zero*data  -20  140]);
   grid  on;
   title('spectrum of I-channel after zero insertion');
   subplot(224);
   plot(20*log(abs(fft(Qzero))));
   axis([0  zero*data   -20 140]);
   grid  on;
   title('spectrum of Q-channel after zero insertion');
   figure(3);
   subplot(221);
   plot(psf);
   axis([200    300     -0.2    0.6]);
   title('time domain response of pulse shaping filter');
   grid  on;
   subplot(222);
   plot(20*log(abs(fft(psf))));
   axis([0  N   -350 50]);
   grid on;
   title('transfer  function  of pulse  shaping filter');
   subplot(223);
   plot(20*log(abs(fft(Ipulse))));
   axis([0  zero*data+N  -250 150]);
   grid on;
   title('spectrum of I-channel after  impulse shaping filter');
   subplot(224);
   plot(20*log(abs(fft(Qpulse))));
   axis([0  zero*data+N -250  150]);
   grid  on;
   title('spectrum of Q-channel  after pluse shaping  filter');
   figure(4)
   subplot(211);
   plot(20*log(abs(fft(Imod))));
   axis([0  zero*data+N  -250 150]);
   grid  on ;
   title('spectrum of I-channel  after modulation');
   subplot(212);
   plot(20*log(abs(fft(Qmod))));
   axis([0  zero*data+N  -250 150]);
   grid  on;
   title('spectrum  of  Q-channel after modulation');
   figure(5)
   subplot(221);
   plot(20*log(abs(fft(Idem))));
   axis([0 zero*data  -200  150]);
   grid on;
   title('spectrum  of I-channel after  demodulation');
   subplot(222);
   plot(20*log(abs(fft(Qdem))));
   axis([0  zero*data+N  -200  150 ]);
   grid  on;
   title('spectrum of Q-channel after demodulation');
   subplot(223);
   plot(20*log(abs(fft(Imat))));
   axis([0  zero*data  -400  200]);
   grid  on;
   title('spectrum  of I-channel  after  matched filter');
   subplot(224);
   plot(20*log(abs(fft(Qmat))));
   axis([0  zero*data  -400  200]);
   grid  on;
   title('spectrum of  Q-channel after matched filter');
   figure(6)
   subplot(221);
   plot(20*log(abs(fft(Isam))));
   axis([0 data  -40  150]);
   grid  on;
   title('spectrum of I-channel after sampler');
   subplot(222);
   plot(20*log(abs(fft(Qsam))));
   axis([0  data -40  150 ]);
   grid  on;
   title('spectrum of Q-channel after  sampler');
   subplot(223);
   plot(20*log(abs(fft(Ifinal))));
   axis([0 data  -40  150]);
   grid on;
   title('spectrum of  I-channel after  decision threshold');
   subplot(224);
   plot(20*log(abs(fft(Qfinal))));
   axis([0 data  -40  150]);
   grid on;
   title('spectrum of  Q-channel after  decision threshold');
   figure(7)
   plot(Isel,Qsel);
   axis([-1.6 1.6  -1.6  1.6]);
   grid  on;
   title('constellation  of  matched  filter  output');
   figure(8)
   plot(Isam,Qsam,'X');
   axis([-1.2  1.2   -1.2  1.2]);
   grid on;
   title('constellation  of  sampler');
   figure(9)
   plot(20*log(abs(fft(final))));
   axis([0  data  0  100]);
   grid  on;
   title('aspectrum  of  final  received  binary  data');
