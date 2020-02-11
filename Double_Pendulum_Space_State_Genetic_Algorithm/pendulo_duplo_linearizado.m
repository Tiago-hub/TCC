%Essa função retorna as funçoes de transferencia e o modelo em espaço de estados do pendulo duplo quando linearizado
%Deseja-se descrever um sistema de pendulo duplo onde as equaçoes diferenciais foram obtidas a partir do modelo 
%Lagrangiano
%As variaveis de estado são: x1->angulo superior (angulo 1), x2=x1p->velocidade angular 1, x2p->aceleração angular 1
% x3->angulo inferior (angulo 2), x4=x3p->velocidade angular 2, x4p->aceleração angular 2
%As demais variaveis são constantes físicas como tamanho do braço do pendulo, localização do centro de gravidade, atrito e aceleração da gravidade
function [Q SS]=pendulo_duplo_linearizado();
syms x1 x2 x3 x4 x1p x2p x3p x4p w1 w2 l1 l2 m1 m2 g b

constantes=[9.7838163 0.25 10e-2 10e-2/2 0.25 10e-2 10e-2/2 0.01]; %vetor com os valores constantes na seguinte ordem: [g m1 w1 l1 m2 w2 l2 b]

eq1=x2p==-(b*x2+x4p*w1*l2*m2*cos(x1-x3)+(m1*l1+m2*w1)*g*sin(x1)+x4^2*w1*l2*m2*sin(x1-x3))/(m1*l1^2+m2*w1^2);
eq2=x4p==-(x2p*w1*cos(x1-x3)-x2^2*w1*sin(x1-x3)+g*sin(x3)+b*x4)/l2;

f1=x2;
f3=x4;

temp1=-(b*x2+x4p*w1*l2*m2*cos(x1-x3)+(m1*l1+m2*w1)*g*sin(x1)+x4^2*w1*l2*m2*sin(x1-x3))/(m1*l1^2+m2*w1^2);
temp2=-(x2p*w1*cos(x1-x3)-x2^2*w1*sin(x1-x3)+g*sin(x3)+b*x4)/l2;
eq1_desacoplada=subs(eq1,x4p,temp2);
f2=solve(eq1_desacoplada,x2p);
eq2_desacoplada=subs(eq2,x2p,temp1);
f4=solve(eq2_desacoplada,x4p);

jacobiano=[diff(f1,x1) diff(f1,x2) diff(f1,x3) diff(f1,x4);
           diff(f2,x1) diff(f2,x2) diff(f2,x3) diff(f2,x4);
           diff(f3,x1) diff(f3,x2) diff(f3,x3) diff(f3,x4);
           diff(f4,x1) diff(f4,x2) diff(f4,x3) diff(f4,x4)];
 
jacobiano_avaliado_zero=subs(jacobiano,[x1 x2 x3 x4],[0 0 0 0]);
A=subs(jacobiano_avaliado_zero,[g m1 w1 l1 m2 w2 l2 b],constantes);
A=eval(A);
B=[1 0;0 0;0 1;0 0];
C=[1 0 0 0;0 0 1 0];
D=zeros(2,2);

SS = ss(A,B,C,D);
[num den]=ss2tf(A,B,C,D,1);
Q11=tf(num(1,:), den);
Q12=tf(num(2,:),den);
[num den]=ss2tf(A,B,C,D,2);
Q21=tf(num(1,:), den);
Q22=tf(num(2,:),den);
Q=[Q11 Q12;Q21 Q22];
end



