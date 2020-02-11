function [filho,VetorQ] = Algoritmo_Genetico(numero_individuos,chance_mutacao,chance_reproducao,limites_ganhos,planta,Q,R,u,n)
    
    a=single(limites_ganhos(1)); %menor valor possível
    b=single(limites_ganhos(2)); %maior valor possível
    rng('shuffle'); %seed do rng
    for k=1:numero_individuos
        for i=1:length(planta.A)
            for j=1:length(planta.A)
                if ((i==j) & (j==1))||((i==j)&(j==3))
                    VetorQ(k).Q(i,j)=a+(b-a)*rand();
                else
                    VetorQ(k).Q(i,j)=0;
                end
            end
        end
    end
    Limite_Sup_Desempenho=1e9;
    
    %testando a população:
    Tf=single(0.01);

    for i=1:numero_individuos
        Qteste=cell2mat(struct2cell(VetorQ(i)));
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
        Qteste=cell2mat(struct2cell(VetorQ(i)));
       Desempenho(1,i)=LQR_calc(x,U,Qteste,R,n,t);
       if Desempenho(1,i) >= Limite_Sup_Desempenho
           Desempenho(1,i)=Limite_Sup_Desempenho;
       end
    end

    %classificação dos indivíduos
    Base_Desempenho=max(Desempenho);    %Define a base do vetor de desempenho a ser normalizado
    Classificapu=(Desempenho./Base_Desempenho); %Como desejo minimizar desempenho, numeros baixos de classifica sao bons desempenhos
    Classifica=sort(Classificapu);
    
    %Seleção da nova geração
    %ira selecionar metade do tamanho do numero de individuos para
    %reproducao, a metade caso seja impar é definida como o numero central
    if mod(numero_individuos,2)==0 
        tam_reproducao=(numero_individuos/2+1); %caso numero de individuos par
    else
        tam_reproducao=(numero_individuos/2+1/2);   %caso numero de individuos impar
    end
    
    for i=1:tam_reproducao
        %if i<=tam_reproducao/4;
         %   Roleta(i)=round(1+(length(Classifica)/2-1)*rand());   %prioriza selecao dos mais aptos reservando 1/4 da populacao selecionada para estar entre a metade melhor
        %else
            Roleta(i) = round(1+(length(Classifica)-1)*rand());
        %end
    end


    for i=1:tam_reproducao
        Nova_Populacao(i)=Classifica(Roleta(i));
    end

    for i=1:tam_reproducao
        for j=1:length(Classificapu)
            if Nova_Populacao(i)==Classificapu(j)
                
                Novo_K(i)=VetorQ(j);
                %break;
            end
        end
    end

    %Cruzamento dos selecionados

    for i=1:numero_individuos
        P1=cell2mat(struct2cell(Novo_K(round(1+(length(Novo_K)-1)*rand()))));
        P2=cell2mat(struct2cell(Novo_K(round(1+(length(Novo_K)-1)*rand()))));
        for j=1:length(P1)
            for k=1:length(P1)
                if ((k==j) & (j==1))||((k==j)&(j==3))
                    if 100*rand()>=chance_reproducao
                        if rand()>=0.5
                            cruzamento(j,k)=P1(j,k);
                        else
                            cruzamento(j,k)=P2(j,k);
                        end
                    else
                        cruzamento(j,k)=P1(j,k);
                    end
                    if 100*rand()<=chance_mutacao    %Chance de mutacao
                        cruzamento(j,k) = a + (b-a).*rand(); %mutacao
                    end
                else
                    cruzamento(j,k)=0;
                end
            end
        end
        filho(i).Q=cruzamento;
    end

end