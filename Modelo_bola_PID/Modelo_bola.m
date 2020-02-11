%Sistema bola-barra e sintonia Ziegler-Nichols
%Autor: Tiago da Costa Ferreira

%Código criado para criar o modelo do sistema Bola-Barra e avaliar a
%sintonia por Ziegler-Nichols

clc 
clear all
close all


%dados gerais
g=9.7838163; %gravidade em BH

%dados da bola
mb=100e-3;
rb=5.5e-2/2;
Jb=(2*mb*rb^2)/5; %momento de inércia

%dados da barra
Lh = 70e-2;
mh = 0.146;
Jh=(mh*Lh^2)/12; %momento de inércia

%função de transferência
lambda = (mb*g)/(Jb/rb+mb);
s=tf('s');
G=(lambda/s^2)

%Espaço de estados
A=[0 0;1 0];
B=[lambda; 0];
C=[0 1];
D=0;
SS = ss(A,B,C,D);

%%
%Sintonia PID por ZN
Kcr=1;
Pcr=2.02;

%P por ZN

Kp1=0.5*Kcr;
Ti1=0;
Td1=0;
Gc1=Kp1;

%PI por ZN

Kp2=0.45*Kcr;
Ti2=Pcr/1.2;
Td2=0;
Gc2=Kp2*(1+1/(s*Ti2)+Td2*s);

%PID por ZN

Kp3=0.6*Kcr;
Ti3=Pcr/2;
Td3=Pcr/8;
Gc3=Kp3*(1+1/(s*Ti3)+((Td3*s)/((Td3/1000)*s+1)));


%%
%step(G)
figure()
step(G) %sistema em malha aberta
figure()
rlocus(feedback(G,1))   %lugar das raizes do sistema com realimentação unitária
figure()
step(feedback(G,1)) %resposta ao degrau sistema com realimentação unitária
xlim([0 10])
figure()
step(feedback(G*Gc1,1)) %resposta ao degrau com controlador P
xlim([0 10])
[y1 t1]=step(feedback(Gc1*G,1));
figure()
step(Gc1/(1+Gc1*G),t1)
title('Esforço de controle P')
xlim([0 10])
stepinfo(feedback(G*Gc1,1))
figure()
step(feedback(G*Gc2,1)) %resposta ao degrau com controlador PI
[y2 t2]=step(feedback(Gc2*G,1));
figure()
step(Gc2/(1+Gc2*G),t2)
title('Esforço de controle PI')
stepinfo(feedback(G*Gc2,1))
figure()
step(feedback(G*Gc3,1)) %resposta ao degrau com controlador PID
[y3 t3]=step(feedback(Gc3*G,1));
figure()
step(Gc3/(1+Gc3*G),t3)
ylim([-1 1])
title('Esforço de controle PID')

stepinfo(feedback(G*Gc3,1))
figure()
pzmap(feedback(G*Gc3,1))
