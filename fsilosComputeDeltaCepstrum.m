function [deltaCepstrum] = fsilosComputeDeltaCepstrum(ceps,deltaWindow)

ceps = ceps(:,1:13);        % Joyner: Se pasa la matriz de Cepstrum entera, pero sólo se necesitan los primeros 13 coeficientes

[filas,columnas]=size(ceps);
deltaCepstrum = ceps;

% JC: En esta funcion se calcula la diferencia entre 2 tramas antes y dos
% tramas después. La corrida en frío del código está escrita en una hoja (ver mis apuntes)

%Cada fila es un conjunto de 13 coef ceps de datos de la trama de la
%palabra a tratar. Hay que trabajar con ellos con la Superfórmula para que
%los nuevos datos de trama representen los DELTACEPS con unos nuevos 13
%coeficientes. Esos nuevos coeficientes se pondrán en una matriz para ser
%añadidos más adelante a los .ascci junto a los CEPS iniciales.

%Superformula.


for t=1:filas   %t es el instante temporal. La fila amos... la columnas son 13 y son los coeficientes.
   
   numerador = zeros(1,columnas); 
   denominador = 0;
      
   for x=1:deltaWindow                  % CALCULO DEL NUMERADOR. El for es el sumatorio.
       
        if (t-x)<1 && (t+x)<filas       % Caso de que estemos en el borde de la primera trama!!
            numerador = numerador + x*(ceps(t+x,:)-ceps(1,:));
        
        elseif (t-x)>1 && (t+x)>filas   % Caso de que estemos en el borde de la ultima trama!!
            numerador = numerador + x*(ceps(filas,:)-ceps(t-x,:));
        
        elseif (t-x)<1 && (t+x)>filas   % Caso de deltaWindowEnorme, creo que nunca entra porque da ERROR MATLAB.
            numerador = numerador + x*(ceps(filas,:)-ceps(1,:));
            disp('DeltaWindow más grande que las tramas de la palabra!!!!');
        
        else  % Caso correcto y normal.
            numerador = numerador + x*(ceps(t+x,:)-ceps(t-x,:));
        end 
        
   end     % Final del numerador.

   for x=1:deltaWindow                  % CALCULO DEL DENOMINADOR. El for es el sumatorio.
       denominador = denominador + x*x;
   end
   denominador = 2*denominador;         % Final del denominador.

   % ESCRIBO LOS DELTACEPSTRUM.
   deltaCepstrum(t,:)=numerador/denominador;
   % disp('Fila escrita.');
   
end
