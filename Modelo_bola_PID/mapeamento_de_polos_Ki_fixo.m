%Mapeamento de pólos

%Autor: Tiago da Costa Ferreira

%Código criado para realizar o mapeamento de polos do Sistema barra-bola
%%
clc
clear all
close all
s=tf('s');
d=.5:.5:2;    %Define a curva sigma=d no plano complexo sigma x jw no qual os pólos ficarão localizados à esquerda
theta=deg2rad(15:15:75);    %Define o ângulo medido com o eixo imaginário que define o cone onde os pólos do sistema devem residir
c=2.85;  %Define o ganho Ki que será fixo
c2=0.1;
w=0.1:0.1:100;
lambda=9.677;   %Constante calculada para definir o sistema de segunda ordem analisado neste programa: G(s) = lambda/s^2

%%
%mapeando para valores no plano complexo tal que as raízes fiquem dentro da
%área à esquerda da curva w=d
Ki=c;
for i=1:length(d)
    G=lambda./(j.*w-d(i)).^2;
    F=-(j.*w-d(i))./G;
    Kd=((real(F)-c).*w+imag(F)*d(i))./(-w.*(d(i)^2+(w).^2));
    Kp=(imag(F).*(d(i)^2-w.^2)+2*(real(F)-c)*d(i).*w)./(-w.*(d(i)^2+w.^2));
    plot(Kd,Kp,'Color','r');
    hold on
end


%%
%mapeando para valores no plano complexo tal que as raízes fiquem dentro do
%cone definido pelo angulo theta medido com o eixo vertical
Ki=c;
for i=1:length(theta)
    G=lambda./(j.*w.*exp(j.*theta(i))).^2;
    E=-(j.*w.*exp(j.*theta(i)))./G;
    Kd=(imag(E).*sin(theta(i))+(real(E)-c).*cos(theta(i)))./(-w.^2.*cos(theta(i)));
    Kp=((c-real(E)).*sin(2.*theta(i))+imag(E).*cos(2*theta(i)))./(w.*cos(theta(i)));
    plot(Kd,Kp,'Color','b');
    grid on
end
%legend('d=1','d=2','d=3','d=4','d=5','\theta=15','\theta=30','\theta=45','\theta=60','\theta=75')
ylabel('Kp')
xlabel('Kd')
title('Mapeamento de polos')
xlim([0 3])
ylim([1 5])
%%
%testa raízes escolhidas 
Ki=c;
Kp=2.5;
Kd=0.875;
%Kp=(0.6447+0.8495)/2;
%Kd=(1.854+2.793)/2;
Tn=Kd/100;
Gc=Kp+Ki/s+Kd*s/(1+Tn*s);
G=lambda/s^2;
figure()
step(feedback(Gc*G,1))
[y t]=step(feedback(Gc*G,1));
figure()
step(Gc/(1+Gc*G),t)
ylim([-5 5])
title('Esforço de controle PID')
stepinfo(feedback(Gc*G,1))
figure()
pzmap(feedback(Gc*G,1))
