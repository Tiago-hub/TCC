%Modelagem sistema bola barra
%retorna o modelo bola barra como fun��o de transfer�ncia G e espa�o de
%estados SS
%Autor: Tiago da Costa Ferreira


function [G SS] = Modelo_bola()
clear all
close all


%dados gerais
g=9.7838163; %gravidade em BH

%dados da bola
mb=100e-3;
rb=5.5e-2/2;
Jb=(2*mb*rb^2)/5; %momento de in�rcia

%dados da barra
Lh = 70e-2;
mh = 0.146;
Jh=(mh*Lh^2)/12; %momento de in�rcia

%fun��o de transfer�ncia
lambda = (mb*g)/(Jb/rb+mb);
s=tf('s');
G=(lambda/s^2);

%Espa�o de estados
x0=[0;0]; %Estado x1 -> velocidade ; x2 -> posi��o
A=[0 0;1 0];
B=[lambda; 0];
C=[0 1];
D=0;
SS = ss(A,B,C,D);