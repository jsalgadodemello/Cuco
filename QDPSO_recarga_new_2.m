function[] = QDPSO()
clc,clear all
global D;

seed =8980809898;
rand('state', seed);

%%%%%%%%%%%%%%%%%%% Parametro da Funçao %%%%%%%%%%%%%%%%

xmin=0;
xmax=100;
nGeracoes=1000;
M=40; %populacao
D=20;  %dimensao

%%%%%%%%%%%%%%%% Parametro do QDPSO %%%%%%%%%%%%%%%%%%

 g = 1.45;

% Poço potencial outro

%  fator1 =  1.4427;

% Poço oscilador harmônico

%  fator2 = 2.0967;



%Incializar Populacao

for i = 1:M
    for j =1:D
        
        x(i,j)= xmin + rand*(xmax-xmin);
    end
    
    pib(i,:) = x(i,:);
    
end

% loop

Convergencia = [];

for loop = 1:nGeracoes
    disp ('loop')
    disp (loop)
        for i = 1:M
        
        %calculo das fitness
        
        fxi=fitness(x(i,:));
       
        fxp(i)=fitness(pib(i,:));
        
        if fxi < fxp(i) 
            pib(i,:) = x(i,:); %pib=melhor do i em todas as iteraçoes; x=individuo i atual
        end
    end
    
    %calcula o p global best
    [fpgb ipgb]= min (fxp);
    pgb = pib(ipgb,:); %best
    
    %calcula novas posicoes
    %Populacao gracao t
    
    for i = 1:M
        
        for j =1:D
            phi1 = rand;
            phi2 = rand;
            pd(j)=(phi1*pib(i,j) + phi2*pgb(j))/(phi1+phi2);
            u = rand;
            Ld = (1/g)*abs(x(i,j)-pd(j));         
            
            if rand > 0.5
                
           %%%%%% Poço Delta potencial %%%%%%%%%
           
              x(i,j) = pd(j)- ((Ld)*(log(1/u)));
           
           %%%% Nova versão com fator1
                
%              x(i,j) = pd(j)- ((Ld)*(fator1)*(log(1/u)));
         
            %%%% Poço Oscilador Harmônico com fator2
                
%              x(i,j) = pd(j)-((Ld)*(fator2)*(sqrt((log(1/u)))));

          
            %%%% Tratamento de fronteira %%%%%
            
%                 % colado na fronteira
%                                 if (x(i,j) > xmax) 
%                                     x(i,j) = xmax;
%                                 end
%                                 if (x(i,j) < xmin ) 
%                                     x(i,j) = xmin;
%                                 end
% %                 
                % rebate aleatorio
                %                 
                %                 if (x(i,j) > xmax ) 
                %                     x(i,j) =  xmax + rand*(xmax-xmin);
                %                 end
                %                 if (x(i,j) < xmin ) 
                %                     x(i,j) = xmin + rand*(xmax-xmin);;
                %                 end
                %                 
                
                % rebate deterministico
%                 %                 
%                                 if (x(i,j) > xmax ) 
%                                     x(i,j) =  xmax + 0.5*(xmax-xmin);
%                                 end
%                                 if (x(i,j) < xmin ) 
%                                     x(i,j) = xmin + 0.5*(xmax-xmin);;
%                                 end
%                 
                
                
            else
    


          %%%%%% Poço Delta potencial %%%%%%%%%
           
              x(i,j) = pd(j)+ ((Ld)*(log(1/u)));
           
           %%%% Nova versão com fator1
%                 
%               x(i,j) = pd(j)+ ((Ld)*(fator1)*(log(1/u)));
         
            %%%% Poço Oscilador Harmônico com fator2
                
%              x(i,j) = pd(j)+((Ld)*(fator2)*(sqrt((log(1/u)))));
%           

      %%%% Tratamento de fronteira %%%%%
      
                % colado na fronteira
%                                 if (x(i,j) > xmax) 
%                                     x(i,j) = xmax;
%                                 end
%                                 if (x(i,j) < xmin ) 
%                                     x(i,j) = xmin;
%                                 end
%                 
                % rebate aleatorio
                %                 
                %                 if (x(i,j) > xmax ) 
                %                     x(i,j) =  xmax + rand*(xmax-xmin);
                %                 end
                %                 if (x(i,j) < xmin ) 
                %                     x(i,j) = xmin + rand*(xmax-xmin);;
                %                 end
                %                 
                
                % rebate deterministico
                %                 
%                                 if (x(i,j) > xmax ) 
%                                     x(i,j) =  xmax + 0.5*(xmax-xmin);
%                                 end
%                                 if (x(i,j) < xmin ) 
%                                     x(i,j) = xmin + 0.5*(xmax-xmin);;
%                                 end
%                 
                
            end
        end 
        
    boro = 1/(fpgb);
    Convergencia = [Convergencia; [loop, boro]];
   
  end % novo loop
    
   
     disp ('boro')
    disp (boro)
%     disp('EC')
%     disp(pgb)
end% programa

end

%---------------------------------------------------------------------
% PLOT FIGURA (EVALUATION x FITNESS)
%----------------------------------------------------------------------

% function [] = figura ();
% figure, hold on
% plot (Convergencia(:,1),Convergencia(:,2))
% 
% hxlabel = xlabel('Generations');
% hylabel = ylabel('Fitness');
% htitle  = title('Convergencia');
% 
% end


function [Av_F] = fitness(vec);
global D;

    
 fit = 0;

    
    %**********************************************************************
    % Gravar os vinte valores em KUATO_SAIDA.dat
    %**********************************************************************
 
%    [vord vind] = sort(x);
    [vord vind] = sort(vec);
  
    %disp (vind)
    
    % Início Alan 27/10/2010 - Modificação para fazer com 10 e 10
    for ix = 1:10
        %vetor1(ix) = x(ix);
        vetor1(ix) = vec(ix);
    end
    for ix = 11:20
        vetor2(ix-10) = vec(ix);
%        vetor2(ix) = x(ix);
    end
%    [vord1 vind1] = sort(vetor1);
    [vord1 vind1] = sort(vetor1);
    [vord2 vind2] = sort(vetor2);

    
    for ix = 1:10

        vetor3(ix) = vind1(ix);
    end
    for ix = 11:20
        vetor3(ix) = vind2(ix-10)+ 10;
    end
    % Fim Alan 27/10/2010

    fid = fopen('KUATO_SAIDA.dat' ,'wt');
    for ix = 1:20
%        cnt = fprintf(fid, '%12.8f\n', vind(ix)); Faz com 20
        cnt = fprintf(fid, '%12.8f\n', vetor3(ix)); % Faz com 10 e 10
    end
    sta = fclose(fid);
    
    
%    %>>>>>>>>>>>>>>>..26-08-2013
%     fid = fopen('KUATO_SAIDAacumulada2.dat' ,'at'); 
%    vetor3(21)=1111111111111111111111111111;
%     for ix = 1:21
%         %cnt = fprintf(fid, '%12.8f\n', vind(ix)); Faz com 20
%         cnt = fprintf(fid, '%12.8f\n', vetor3(ix)); % Faz com 10 e 10
%     end
%     
% 
%     sta = fclose(fid);
%    
    
    
    
    %**********************************************************************    
    % Executar o KUATO_FORMIGA.exe
    %**********************************************************************
    
    [sta, wlst] = system ('KUATO_FORMIGA.exe');
    if (sta ~= 0)
        m_Fitness = +Inf; disp (wlst); pause(5), return    
    end
    
    %**********************************************************************    
    % Simular KUATO_FORMIGA.exe
    %**********************************************************************

    % fid = fopen('KUATO_APTIDAO.dat' ,'wt');
    % m_Fitness = rand; cnt = fprintf(fid, '%12.8f\n', m_Fitness);
    % m_Pico    = rand; cnt = fprintf(fid, '%12.8f\n', m_Pico);
    % m_Boro    = rand; cnt = fprintf(fid, '%12.8f\n', m_Boro);
    % sta = fclose(fid);
    
    %**********************************************************************    
    % Ler arquivo KUATO_APTIDAO.dat com Fitness, Pico e Boro
    %**********************************************************************    

    [fid, mensagem] = fopen('KUATO_APTIDAO.dat','rt');
    tline = fgets(fid); m_Fitness = str2num(tline);
    tline = fgets(fid); m_Pico    = str2num(tline);
    tline = fgets(fid); m_Boro    = str2num(tline);
    sta = fclose(fid);
%     fit=m_Fitness;
      if m_Pico ~= 10
          m_Pico = m_Pico /1000;
      end
      fit = m_Pico;
      if fit < 1.395
          fit = 1 / m_Boro;
      end
  
%     disp('fittttttt')
%     disp(fit)
%     
%     disp('pico')
%     disp(m_Pico)
%     
%     disp('boro')
%     disp(m_Boro)
%     
%     
%     if m_Pico<1.395
%     disp('if m_Pico<1.395')
%     %disp(m_Pico)
%     disp(m_Boro)
%     %disp(fit)
%     end
    if m_Boro > 1400 & m_Pico<1.395
   %disp('if m_Pico<1.395')
  
   disp('m_Pico')
    disp(m_Pico)
    disp('m_Boro')
    disp(m_Boro)
    disp('fit')
    disp(fit)
    disp('vetor3')
    disp(vetor3)
    disp('vind1')
    disp(vind1)
    disp('vind2')
    disp(vind2)
    pause 
    end
        Av_F =  fit;
  


return

end