clc
close all
clear all
%����chirp�ź�
%��������
Tp=20*10^(-6);%������
B=150*10^6;%�����źſ��
fs=200*10^6;%������
kr=B/Tp;%��Ƶϵ��
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
%��/���ز��źŲ�������
c=3*10^8;
H=3000;%���и߶�
lambda=0.05;%�ز�����
A=0.025;%��λ�������
Rmin=15000;%��ʼ��������
Rmax=15312;%���ܼ������
PRF=1000;
PRT=1/1000;
v=100;%�ɻ������ٶ�
%��Ŀ���������
R_1=15200;%��Ŀ��1�ͺ����Ĵ�ֱ����
x_1=0;%��Ŀ���غ�������
R_2=15200;%��Ŀ��1�ͺ����Ĵ�ֱ����
x_2=3;%��Ŀ���غ�������
R_3=15250;%��Ŀ��1�ͺ����Ĵ�ֱ����
x_3=0;%��Ŀ���غ�������
Ptarget=[x_1,R_1
         x_2,R_2
         x_3,R_3];
%��������ͷ�λ����
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
xlabel('������ /km');
ylabel('��λ�� /m');
title('����ز�����')
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
xlabel('������ /km');
ylabel('��λ�� /m');
title('3����ز�����')