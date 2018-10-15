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
to=0.1;
%intervalo de tiempo
t=-to:ts:to;

%Señal moduladora
m_t=sinc(100*t);
%Señal portadora
c_t=cos(2*pi*fc*t).*double(abs(t)<=to);

figure(1)
%Grafica de m(t)
subplot(2,1,1);
plot(t,m_t,'color',[0 .8 .8],'linewidth',2)
title(' Grafica de $$ m(t) $$','interpreter','latex','fontsize',21)

grid on

%Proceso para calcular el espectro de magnitud de m(t) => |M(w)|

w=-1000:1000;
M=0;
n=0;
for tt=t
    n=n+1;
    M=M+m_t(n).*exp(-1j.*w.*tt)*ts;
end

%Grafica del espectro de magnitud de m(t) => |M(w)|
subplot(2,1,2);
plot(w,abs(M),'color',[0.01 .46 .99],'linewidth',1.75)
title(' Grafica de $$ |M(w)| $$','interpreter','latex','fontsize',21)
grid on

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

figure(2)
%Calculo de y_tdsb_sc
y_tdsb_sc=m_t.*c_t;

%Grafica de y(t)_DSB-SC
subplot(2,1,1);%, 81%, 82%
plot(t,y_tdsb_sc,'color','m','linewidth',1.75)
title(' Grafica de $$ y(t)_{DSB-SC} $$','interpreter','latex','fontsize',21)
grid on

%Proceso para calcular el espectro de magnitud de  y(t)_DSB-SC
w=-2500:2500;
Y_dsb_sc=0;
n=0;
for tt=t
    n=n+1;
    Y_dsb_sc=Y_dsb_sc+y_tdsb_sc(n).*exp(-1j.*w.*tt)*ts;
end

%Grafica del espectro de magnitud de  y(t)_DSB-SC
subplot(2,1,2);
plot(w,abs(Y_dsb_sc),'color',[.9 .71 .15],'linewidth',1.75);
title(' Grafica de $$ |Y(t)_{DSB-SC} |$$','interpreter','latex','fontsize',21)
grid on

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Demodulación coherente
figure(3);
%Señal demodulada
subplot(2,1,1)
r_t= y_tdsb_sc.*cos(2.*pi.*fc.*t);
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

figure(4)
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
figure(5)
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