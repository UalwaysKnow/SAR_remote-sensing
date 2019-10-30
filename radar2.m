clc
close all
clear all
%基带chirp信号
%参数设置
Tp=20*10^(-6);%脉冲宽度
B=150*10^6;%发射信号宽度
fs=200*10^6;%采样率
kr=B/Tp;%调频系数
t=-10*10^(-6):1/fs:10*10^(-6);
s=rectpuls(t,Tp).*exp(j*pi*kr*(t.^2));
figure(1);
subplot(311);
plot(t,real(s));
subplot(312);
S=fft(s);
plot(t,fftshift(abs(S)));
subplot(313);
plot(t,phase(S));
%%
%单/多点回波信号参数设置
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
Srmn=zeros(Na,Nr);
xT=Ptarget(1,1);
R_n=Ptarget(1,2);
Ls=A*R_n;
for n=1:Na
    Rn=sqrt(R_n^2+(xT-Azi(n))^2);
    sr=rectpuls((Azi(n)-xT),Ls)*rectpuls((t_r-2*Rn/c-Tp/2),Tp) ...
        .*exp(1j*pi*kr*(t_r-2*Rn/c-Tp/2).^2).*exp(-1j*4*pi*Rn/lambda);
    Srmn(n,:) = Srmn(n,:) + sr;
end
figure(2);
imagesc(Rng/1000,Azi,real(Srmn));
colormap gray;
xlabel('距离向 /km');
ylabel('方位向 /m');
title('单点回波数据')
%%
Srmn=zeros(Na,Nr);
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

figure(3);
imagesc(Rng/1000,Azi,real(Srmn));
colormap gray;
xlabel('距离向 /km');
ylabel('方位向 /m');
title('3个点回波数据')