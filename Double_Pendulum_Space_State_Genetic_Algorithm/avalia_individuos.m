function melhor_individuo=avalia_individuos(u,planta,populacao,R,n,Limite_Sup_Desempenho)
    numero_individuos=length(populacao)
    %testando a população:
    Tf=single(0.01);

    for i=1:numero_individuos
        Qteste=(populacao(i).Q);
        [K,S,E]=lqr(planta,Qteste,R,n);
        filho(i).K=K;
        MF=ss(planta.A-planta.B*K,planta.B,planta.C,planta.D);
        temp=size(planta.B);
        t=0:abs(max(real(eig(planta.A))))/100:6*abs(max(real(eig(planta.A))));
        U=u.*ones(length(t),temp(2));
        [y,t,x]=lsim(MF,U,t);
    end

    Desempenho=zeros(1,numero_individuos);  %aloca memoria para vetor de desempenho
    %U=u.*ones(size(t));
    %avaliação da população
    for i=1:numero_individuos
       Qteste=(populacao(i).Q);
       Desempenho(1,i)=LQR_calc(x,U,Qteste,R,n,t);
       if Desempenho(1,i) >= Limite_Sup_Desempenho
           Desempenho(1,i)=Limite_Sup_Desempenho;
       end
    end
    [M melhor_individuo]=min(Desempenho)
end
