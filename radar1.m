clear all;
clc;
X=20:1:45;
A=zeros(1,26);
B=zeros(1,26);
C=zeros(1,26);
D=zeros(1,26);
f=9.35*10^9;%��Ƶ9.35GHz
lambda=3*10^8/f;
B=120*10^6;%�źŴ���120MHz
Pav=4.534;%ƽ�����书��4.534W
PRF=1429;%PRF=1429Hz
Tp=4*10^(-6);%�������4us
B=120*10^6;%����120MHz
theta_a=1.8*pi/180;%������λ����1.8
theta_r=8.5*pi/180;%������������8.5
Ae=lambda^2/(theta_a*theta_r);
G=4*pi*Ae/lambda^2;
pr=3*10^8/(2*B);
N=4096;%������4096
Va=200;%�ɻ��ٶ�200m/s
H=5*10^3;%�ɻ��߶�5km
L=10^(-5/10);%ϵͳ���-5dB
Fn=10^(5/10);%��Ч����ϵ��5dB
K=1.3806*10^(-23);T=290;
j=1;
for i=20:45
    alpha=i*pi/180;%���ӽǷ�Χ20-45
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
xlabel('���ӽǽǶ�');
ylabel('NERZ/dB');
title('NERZ�仯����');
grid on;
figure(2);
plot(X,B);
xlabel('���ӽǽǶ�');
ylabel('�ϳɿ׾�ʱ��/s');
title('�ϳɿ׾�ʱ��仯����')
grid on;
figure(3);
plot(X,C);
xlabel('���ӽǽǶ�');
ylabel('�ؾ�ֱ���/m');
title('�ؾ�ֱ��ʱ仯����')
grid on;
figure(4);
plot(X,D);
xlabel('���ӽǽǶ�');
ylabel('�������/m');
title('������ȱ仯����')
grid on;
