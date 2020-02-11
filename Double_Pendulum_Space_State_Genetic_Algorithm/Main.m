
clc 
clear all
close all
tic
iteracoes=100;
numero_individuos=200;
chance_mutacao=1;
chance_reproducao=80;
[Matriz_tf,SS]=pendulo_duplo_linearizado();
planta=Matriz_tf(1,1);
limites_ganhos=[.1 1e2];
Q=[1 0 0 0; 0 0 0 0;0 0 1 0;0 0 0 0];
R=[1 0; 0 1];
n=[0 0;0 0;0 0;0 0];
u=1;
for geracoes=1:iteracoes
    if geracoes~=1
        clear all
        close all
        clc
        load('Filhos.mat');
    end
 
[filho,pop_inicial]=Algoritmo_Genetico(numero_individuos,chance_mutacao,chance_reproducao,limites_ganhos,SS,Q,R,u,n);

save('Filhos.mat');

end
toc
%%
load('Filhos.mat');

%testando a população:
melhor=avalia_individuos(u,SS,filho,R,n,1e9);
best=ss(SS.A-SS.B*filho(melhor).K,SS.B,SS.C,SS.D);
t=0:0.01:10;
u1=0.2;
u2=0.1;
u=[u1.*ones(length(t),1) u2.*ones(length(t),1)];
[y,t,x]=lsim(best,u,t);
plot(t,x(:,1),t,x(:,2),t,x(:,3),t,x(:,4))
melhorQ=filho(melhor).Q
melhorK=filho(melhor).K
polos=eig(best.A)
toc
