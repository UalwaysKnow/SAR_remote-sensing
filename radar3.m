clc
clear all
close all
%����chirp�ź�
%��������
Tp=20*10^(-6);%������
B=150*10^6;%�����źſ��
fs=200*10^6;%������
kr=B/Tp;%��Ƶϵ��
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
Srmn=zeros(Na,Nr);%�ز�����
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
xlabel('������ /km');
ylabel('��λ�� /m');
title('3����ز�����');
%%
%����ѹ��
Srff = fft(Srmn.').';
s_rej = rectpuls(t_r-2*Rmin/c-Tp/2).*exp(1j*pi*kr*(t_r-2*Rmin/c-Tp/2).^2);
S_rej = conj(fft(s_rej));
Srff0 = Srff.*S_rej;
Srmn0 = ifft(Srff0.').';%����ѹ�����ʱ������

figure(2);
imagesc(Rng/1000,Azi,abs(Srmn0));
colormap gray;
xlabel('������ /km');
ylabel('��λ�� /m');
title('����ѹ���������');
%%
%Ƶ������������
fd_r = linspace(-1*fs/2 , 1*fs/2,Nr);
fd_a = linspace(-1*PRF/2 , 1*PRF/2,Na);

Srmn0_T = fftshift(Srmn0).';
%��λ��fft
Srff_T = zeros(Nr,Na);
for n=1:Nr
    Srff_T(n,:) = fft(Srmn0_T(n,:));
end
Srmn1 = Srff_T.';%RCMCǰ�ľ��������������

figure(3);
imagesc(Rng/1000,fd_a,fftshift(abs(Srmn1)));
colormap gray;
xlabel('������ /km');
ylabel('��λ�� /Hz');
title('RCMCǰ��RD������');
%%
%�����㶯����
ffff = fftshift(fft(fftshift(Srff0)));%RCMCǰ�Ķ�άƵ�������
%ffff0 = zeros(Na,Nr);%RCMC��Ķ�άƵ������
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
Srmn2 = fftshift(ifft((ffff0.'))).';%RCMC���RD������

figure(4);
imagesc(Rng/1000,fd_a,fftshift(abs(Srmn2)));
colormap gray;
xlabel('������ /km');
ylabel('��λ�� /Hz');
title('RCMC���RD������');
%%
%����RCMC��ʱ��ı仯
Srmn2_T = fftshift(Srmn2).';
%��λ��ifft
Srmn777_T = zeros(Nr,Na);
for n=1:Nr
    Srmn777_T(n,:) = ifft(Srmn2_T(n,:));
end
Srmn777 = fftshift(Srmn777_T).';%��������

figure(5);
imagesc(Rng/1000,Azi,fftshift(abs(Srmn777)));
colormap gray;
xlabel('������ /km');
ylabel('��λ�� /m');
title('RCMC��ʱ������');
%%
%%Azimuth compression
%Lsar = A*R_n;
%Tsar = Lsar / v;
%sn=linspace((-Lsar/2)/v,(3+Lsar/2)/v,Na);
%Ka = -2 * v^2 / (lambda * R_n); % �����յ�Ƶ��
%Refa = exp(j * pi * Ka * sn.^2) .* (abs(sn) < Tsar/2); % ��λѹ���ο�����
%REFA = fftshift(fft(fftshift(Refa)));
%Sa_RCMC0 = Srmn2.*(conj(REFA).'*ones(1, Nr)); % �Ծ����㶯У����ľ��뷽λѹ�����
%Sa_RCMC = fftshift(ifft(fftshift(Sa_RCMC0)));

%p0_2fu=exp(j*pi/Ka*Fa.^2);%��λ��ѹ������
%s2rcac_fut=Srmn2.*conj(p0_2fu);%��λѹ��
%s2rcac_fuf=fftshift(fft(fftshift(s2rcac_fut)));%���뷽λѹ����Ķ�άƵ��
%s2rcac_ut=fftshift(ifft(fftshift(s2rcac_fut).')).';%��λ��IFFT

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
xlabel('������ /km');
ylabel('��λ�� /m');
title('���ճ�������');





