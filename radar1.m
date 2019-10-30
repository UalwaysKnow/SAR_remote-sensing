clear all;
clc;
X=20:1:45;
A=zeros(1,26);
B=zeros(1,26);
C=zeros(1,26);
D=zeros(1,26);
f=9.35*10^9;%载频9.35GHz
lambda=3*10^8/f;
B=120*10^6;%信号带宽120MHz
Pav=4.534;%平均发射功率4.534W
PRF=1429;%PRF=1429Hz
Tp=4*10^(-6);%脉冲带宽4us
B=120*10^6;%带宽120MHz
theta_a=1.8*pi/180;%波数方位向宽度1.8
theta_r=8.5*pi/180;%波数距离向宽度8.5
Ae=lambda^2/(theta_a*theta_r);
G=4*pi*Ae/lambda^2;
pr=3*10^8/(2*B);
N=4096;%采样数4096
Va=200;%飞机速度200m/s
H=5*10^3;%飞机高度5km
L=10^(-5/10);%系统损耗-5dB
Fn=10^(5/10);%等效噪声系数5dB
K=1.3806*10^(-23);T=290;
j=1;
for i=20:45
    alpha=i*pi/180;%下视角范围20-45
    R=H/cos(alpha);
    NESZ=((4*pi)^2*R^3*Va*K*T*Fn)/(Pav*G*Ae*lambda*pr*L);
    A(j)=NESZ;
    Ls=R*theta_a;
    Ts=Ls/200;
    B(j)=Ts;
    pgr=pr/sin(alpha);
    C(j)=pgr;
    Rn=H*tan(alpha-theta_r/2);
    Rf=H*tan(alpha+theta_r/2);
    D(j)=Rf-Rn;
    j=j+1;
end
figure(1);
plot(X,A);
xlabel('下视角角度');
ylabel('NERZ/dB');
title('NERZ变化曲线');
grid on;
figure(2);
plot(X,B);
xlabel('下视角角度');
ylabel('合成孔径时间/s');
title('合成孔径时间变化曲线')
grid on;
figure(3);
plot(X,C);
xlabel('下视角角度');
ylabel('地距分辨率/m');
title('地距分辨率变化曲线')
grid on;
figure(4);
plot(X,D);
xlabel('下视角角度');
ylabel('测绘带宽度/m');
title('测绘带宽度变化曲线')
grid on;
