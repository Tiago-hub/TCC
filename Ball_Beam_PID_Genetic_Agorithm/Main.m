
clc 
clear all
close all
tic
iteracoes=50;
[planta,SS]=Modelo_bola();
numero_individuos=20;
chance_mutacao=1;
chance_reproducao=80;
limites_ganhos=[0.01 2];
Q=[1 0; 0 1];
R=1;
n=0;
u=1;
for geracoes=1:iteracoes
    if geracoes~=1
        clear all
        close all
        clc
        load('Filhos.mat');
    end
 
[filho,pop_inicial]=Algoritmo_Genetico(numero_individuos,chance_mutacao,chance_reproducao,limites_ganhos,planta,Q,R,u,n);

save('Filhos.mat');

end
toc
%%
load('Filhos.mat');

%testando a população:
Tf=single(0.01);

for i=1:numero_individuos
    controlador=pid(filho(i,1),filho(i,2),filho(i,3),Tf);
    MF=feedback(controlador*planta,1);
    [y,t]=step(MF);
    y=single(y);
    t=single(t);
    yp=diff(y)/(t(2)-t(1)); %calculo da derivada da posição y
    y(length(y))=[];
    t(length(t))=[];
    for j=1:length(y)
        Y(j,i)=y(j);
        T(j,i)=t(j);
        Yp(j,i)=yp(j);
    end
end


Desempenho=zeros(1,numero_individuos);
Limite_Sup_Desempenho=1e9;

for i=1:numero_individuos
   Desempenho(1,i)=LQR_calc([Y(:,i) Yp(:,i)],u,Q,R,n,t);
   if Desempenho(1,i) >= Limite_Sup_Desempenho
       Desempenho(1,i)=Limite_Sup_Desempenho;
   end
end
end
[M,I]=min(Desempenho);
    controlador=pid(filho(I,1),filho(I,2),filho(I,3),Tf);
    MF=feedback(controlador*planta,1);
    [y,t]=step(MF);
    esforco=controlador/(1+planta*controlador);
    [y2,t2]=step(esforco);
    erro=1-y;
    plot(t,y)
    figure()
    plot(t,erro)
    figure()
    plot(t2,y2)
filho(I,:)
toc
