function [filho,K] = Algoritmo_Genetico(numero_individuos,chance_mutacao,chance_reproducao,limites_ganhos,planta,Q,R,u,n)
    
    a=single(limites_ganhos(1)); %menor valor possível
    b=single(limites_ganhos(2)); %maior valor possível
    rng('shuffle'); %seed do rng
    K = a + (b-a).*rand(numero_individuos,3);    %matriz de ganhos kp ki e kd com numeros aleatorios limitados no intervalo [a,b]
    K=single(K);
    Limite_Sup_Desempenho=1e9;
    
    %testando a população:
    Tf=single(0.01);

    for i=1:numero_individuos
        controlador=pid(K(i,1),K(i,2),K(i,3),Tf);
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

    Desempenho=zeros(1,numero_individuos);  %aloca memoria para vetor de desempenho
    %U=u.*ones(size(t));
    %avaliação da população
    for i=1:numero_individuos
       Desempenho(1,i)=LQR_calc([Y(:,i) Yp(:,i)],u,Q,R,n,t);
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
                Novo_K(i,:)=K(j,:);
                %break;
            end
        end
    end

    %Cruzamento dos selecionados

    for i=1:numero_individuos
        P1=Novo_K(round(1+(length(Novo_K)-1)*rand()),:);
        P2=Novo_K(round(1+(length(Novo_K)-1)*rand()),:);
        for j=1:3
            if 100*rand()>=chance_reproducao
                if rand()>=0.5
                    filho(i,j)=P1(j);
                else
                    filho(i,j)=P2(j);
                end
            else
                filho(i,j)=P1(j);
            end
            if 100*rand()<=chance_mutacao    %Chance de mutacao
                filho(i,j) = a + (b-a).*rand(); %mutacao
            end
        end
    end

end