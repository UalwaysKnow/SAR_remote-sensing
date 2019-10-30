clc
clear all
close all
%基带chirp信号
%参数设置
Tp=20*10^(-6);%脉冲宽度
B=150*10^6;%发射信号宽度
fs=200*10^6;%采样率
kr=B/Tp;%调频系数
c=3*10^8;
H=3000;%飞行高度
lambda=0.05;%载波波长
A=0.025;%方位向波束宽度
Rmin=15000;%初始采样距离
Rmax=15312;%可能计算错误
PRF=1000;
PRT=1/1000;
v=100;%飞机飞行速度
%点目标参数设置
R_1=15200;%点目标1和航迹的垂直距离
x_1=0;%点目标沿航向坐标
R_2=15200;%点目标1和航迹的垂直距离
x_2=3;%点目标沿航向坐标
R_3=15250;%点目标1和航迹的垂直距离
x_3=0;%点目标沿航向坐标
Ptarget=[x_1,R_1
         x_2,R_2
         x_3,R_3];
%距离向轴和方位向轴
x0=-195;
Nr=(Rmax+Tp*c/2-Rmin)*2*fs/c;
t_r=[2*Rmin/c+(0:(Nr-1))/fs];
Rng = t_r*c/2;
Na = -ceil(2*x0/v/PRT);
t_a = x0/v + (0:(Na-1))*PRT;
Azi = t_a*v;
%%
Srmn=zeros(Na,Nr);%回波数据
for i=1:3
    xT=Ptarget(i,1);
    R_n=Ptarget(i,2);
    Ls=A*R_n;
    for n=1:Na
        Rn=sqrt(R_n^2+(xT-Azi(n))^2);
        sr=rectpuls((Azi(n)-xT),Ls)*rectpuls((t_r-2*Rn/c-Tp/2),Tp) ...
            .*exp(1j*pi*kr*(t_r-2*Rn/c-Tp/2).^2).*exp(-1j*4*pi*Rn/lambda);
        Srmn(n,:) = Srmn(n,:) + sr;
    end
end
figure(1);
imagesc(Rng/1000,Azi,real(Srmn));
colormap gray;
xlabel('距离向 /km');
ylabel('方位向 /m');
title('3个点回波数据');
%%
%距离压缩
Srff = fft(Srmn.').';
s_rej = rectpuls(t_r-2*Rmin/c-Tp/2).*exp(1j*pi*kr*(t_r-2*Rmin/c-Tp/2).^2);
S_rej = conj(fft(s_rej));
Srff0 = Srff.*S_rej;
Srmn0 = ifft(Srff0.').';%距离压缩后的时域数据

figure(2);
imagesc(Rng/1000,Azi,abs(Srmn0));
colormap gray;
xlabel('距离向 /km');
ylabel('方位向 /m');
title('距离压缩后的数据');
%%
%频率轴坐标设置
fd_r = linspace(-1*fs/2 , 1*fs/2,Nr);
fd_a = linspace(-1*PRF/2 , 1*PRF/2,Na);

Srmn0_T = fftshift(Srmn0).';
%方位向fft
Srff_T = zeros(Nr,Na);
for n=1:Nr
    Srff_T(n,:) = fft(Srmn0_T(n,:));
end
Srmn1 = Srff_T.';%RCMC前的距离多普勒域数据

figure(3);
imagesc(Rng/1000,fd_a,fftshift(abs(Srmn1)));
colormap gray;
xlabel('距离向 /km');
ylabel('方位向 /Hz');
title('RCMC前的RD域数据');
%%
%距离徙动矫正
ffff = fftshift(fft(fftshift(Srff0)));%RCMC前的二维频域的数据
%ffff0 = zeros(Na,Nr);%RCMC后的二维频域数据
Fa = fd_a'*ones(1,Nr);
Fr = ones(Na,1)*fd_r;
%for i=1:3
%    xT=Ptarget(i,1);
%    R_n=Ptarget(i,2);
%    Ls=A*R_n;
fdr = 2 * v^2 ./ (lambda *[Rmin+(0:Nr-1)*(Rmax-Rmin)/Nr]);
fdr = ones(Na,1)*fdr;
Fc = c/lambda;
RMC=exp(-(1j*pi/Fc^2./fdr.*(Fa.*Fr).^2-1j*pi/Fc./fdr.*Fa.^2.*Fr));
ffff0 = ffff.*RMC;
%    ffff0 = ffff0 + ff;
%end
Srmn2 = fftshift(ifft((ffff0.'))).';%RCMC后的RD域数据

figure(4);
imagesc(Rng/1000,fd_a,fftshift(abs(Srmn2)));
colormap gray;
xlabel('距离向 /km');
ylabel('方位向 /Hz');
title('RCMC后的RD域数据');
%%
%测试RCMC后时域的变化
Srmn2_T = fftshift(Srmn2).';
%方位向ifft
Srmn777_T = zeros(Nr,Na);
for n=1:Nr
    Srmn777_T(n,:) = ifft(Srmn2_T(n,:));
end
Srmn777 = fftshift(Srmn777_T).';%最终数据

figure(5);
imagesc(Rng/1000,Azi,fftshift(abs(Srmn777)));
colormap gray;
xlabel('距离向 /km');
ylabel('方位向 /m');
title('RCMC后时域数据');
%%
%%Azimuth compression
%Lsar = A*R_n;
%Tsar = Lsar / v;
%sn=linspace((-Lsar/2)/v,(3+Lsar/2)/v,Na);
%Ka = -2 * v^2 / (lambda * R_n); % 多普勒调频率
%Refa = exp(j * pi * Ka * sn.^2) .* (abs(sn) < Tsar/2); % 方位压缩参考函数
%REFA = fftshift(fft(fftshift(Refa)));
%Sa_RCMC0 = Srmn2.*(conj(REFA).'*ones(1, Nr)); % 对距离徙动校正后的距离方位压缩结果
%Sa_RCMC = fftshift(ifft(fftshift(Sa_RCMC0)));

%p0_2fu=exp(j*pi/Ka*Fa.^2);%方位向压缩因子
%s2rcac_fut=Srmn2.*conj(p0_2fu);%方位压缩
%s2rcac_fuf=fftshift(fft(fftshift(s2rcac_fut)));%距离方位压缩后的二维频谱
%s2rcac_ut=fftshift(ifft(fftshift(s2rcac_fut).')).';%方位向IFFT

%fa_azimuth_MF = fd_a;
%R0_RCMC = (c/2).*t_r;
%Ka = 2*v^2./(lambda.* R0_RCMC);  
%Ka_1 = 1./Ka;
%Haz = exp( -1j*pi.*(((fa_azimuth_MF).').^2*Ka_1) );
%S_rd_c = Srmn2.*Haz; 
%s_ac = ifft(S_rd_c,[],1); 

Srmn3 = zeros(Na,Nr);
fdr = 2*v^2/(lambda*R_n);
for n=1:Na
    fa = -2*v*Azi(n)/(lambda*sqrt(R_n+(Azi(n))^2));
    Ha = exp(-1j*pi*fa^2/fdr);
    Srmn3(n,:) = Srmn2(n,:).*conj(Ha);
end
Srmn4 = fftshift(ifft(fftshift(Srmn3)));
figure(6);
imagesc(Rng/1000,Azi,abs(Srmn4));
colormap gray;
xlabel('距离向 /km');
ylabel('方位向 /m');
title('最终成像数据');





