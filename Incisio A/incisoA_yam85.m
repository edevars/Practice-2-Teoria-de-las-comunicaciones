clc;
clear all;
close all;

%Frecuencia de muestreo 
fs=8000;
%Frecuencia de la señal portadora
fc=250;
%tiempo de muestreo
ts=1/fs;

%limite de tiempo de la señal
to=0.15;
%intervalo de tiempo
t=0:ts:to;
%Indice de figuras
k=0;

%Vector de frecuencias
w=linspace(-fs/2,fs/2,20000)*2*pi;

%Señal moduladora
m_t=double(0<=t & t<=to/3)-2.*double(to/3<t & t<=2.*to/3);
%Señal portadora
c_t=cos(2*pi*fc*t);

k=k+1;
figure(k)
%Grafica de m(t)
subplot(2,1,1);
plot(t,m_t,'color',[0 .8 .8],'linewidth',2)
title(' Grafica de $$ m(t) $$','interpreter','latex','fontsize',21)
axis([-.02 to+.01 -2.2 1.2])
grid on

%Proceso para calcular el espectro de magnitud de m(t) => |M(w)|

M=fftshift(fft(m_t,20000))*ts;

%Grafica del espectro de magnitud de m(t) => |M(w)|
subplot(2,1,2);
plot(w/(2*pi),abs(M),'color',[0.01 .46 .99],'linewidth',1.75)
title(' Grafica de $$ |M(w)| $$','interpreter','latex','fontsize',21)
grid on
grid minor
xlim([-500 500])


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%u=0.85
k=k+1;
figure(k)
mp=abs(min(m_t));
u=0.85;
A=mp/u;

%Grafica de la señal moduldada y(t)_am con u=0.85
subplot(2,1,1)
yam=(A+m_t).*c_t;
plot(t,yam,'linewidth',2,'color',[.24 .63 .33])
title(' Grafica de la senial modulada $$   y(t)_{AM}$$, $$ \mu=0.85 $$','interpreter','latex','fontsize',21)
grid on

Yam=0;
n=0;
for tt=t
    n=n+1;
    Yam=Yam+yam(n).*exp(-1j.*w.*tt)*ts;
end

%Espectro de magnitud de  |Y(t)_am| con u=0.85
subplot(2,1,2);
plot(w/(2*pi),abs(Yam),'linewidth',2,'color',[.46 .93 0]);
title(' Grafica del espectro de magnitud   $$  |Y(t)_{AM}| $$, $$ \mu=0.85 $$','interpreter','latex','fontsize',21)
grid on

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


%Demodulación coherente
k=k+1;
figure(k)
%Señal demodulada
subplot(2,1,1)
r_t= yam.*cos(2.*pi.*fc.*t);
plot(t,r_t);
title(' Grafica de la senial demodulada  $$  r(t) $$','interpreter','latex','fontsize',21)
grid on

w=(-1000:1000)*2*pi;

R=0;
n=0;
for tt=t
    n=n+1;
    R=R+r_t(n).*exp(-1j.*w.*tt)*ts;
end

%Grafica del espectro de magnitud de m(t) => |M(w)|
subplot(2,1,2);
plot(w/(2*pi),abs(R),'color',[0.01 .46 .99],'linewidth',1.75)
title(' Espectro de la senial demodulada $$ |R(w)| $$','interpreter','latex','fontsize',21)
grid on

%Creacion del filtro
FPB=zeros(1,length(w));

I=find(abs(w)<=200*2*pi);

FPB(I)=1;

k=k+1;
figure(k)
subplot(2,1,1)
plot(w/(2*pi),abs(FPB),'color','m','linewidth',1.75);
title(' Filtro pasa bajas $FPB$','interpreter','latex','fontsize',21)
grid on;

%Señal filtrada
M1=R.*FPB;
subplot(2,1,2)
plot(w/(2*pi),abs(M1),'color','c','linewidth',1.75);
title(' Senial Filtrada','interpreter','latex','fontsize',21)
grid on;

%Transformada inversa de Fourier
k=k+1;
figure(k)
dw=w(2)-w(1);
m1=0;
n=0;

for ww=w
    n=n+1;
    m1=m1+M1(n)*exp(1j*ww*t)*dw/(2*pi);
end

plot(t,real(m1));
title(' Senial Recuperada','interpreter','latex','fontsize',21)
grid on;