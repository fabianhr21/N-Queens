
clear all;clc;

muestras = 100;                 %Tamaño de población
reinas = 8;                     %Número de reinas
gps = zeros(2,reinas);           %Matriz de ubicación; fila 2 es x y fila 3 es y
aptitud = ones(1,muestras);    %Aptitud por muestra
size = 8;                      %Tamaño del tablero size*size

tic
%Generación de población inicial
for i = 1:muestras
    gps(1,1,i) = 0;
    gps(2,1,i) = 0;
    for j = 1:reinas
        x = randi(size);
        y = randi(size);
        for k = 1:reinas
            while x == gps(1,k,i) && y == gps(2,k,i)
                x = randi(size);
                y = randi(size);   
            end
        end    
        gps(1,j,i) = x;
        gps(2,j,i) = y;
    end  
    aptitud(i) = apt(i,reinas,gps);     %Funcion aptitud
end

%%% Algoritmo genético %%%
for i = 1:100000
    %   Selección
    [~,ind] = mink(aptitud,2);          %Par de individuos óptimos para reproducirse 
    padre1 = gps(:,:,ind(1));
    padre2 = gps(:,:,ind(2));
    %   Reproducción y mutacion
    cruz = randi(reinas);
    hijo1 = [padre1(:,[1:cruz]), padre2(:,[cruz+1:reinas])];    %Hijo1
    hijo1 = mutacion(hijo1,reinas,size);                             %Mutacion hijo 1
    hijo2 = [padre2(:,[1:cruz]), padre1(:,[cruz+1:reinas])];    %Hijo2
    hijo2 = mutacion(hijo2,reinas,size);                             %Mutacion Hijo 2
    % Reemplazo 
    [~,indm] = maxk(aptitud,2);
    gps(:,:,indm(1)) = hijo1;
    gps(:,:,indm(2)) = hijo2;
    aptitud(indm(1)) = apt(indm(1),reinas,gps);
    aptitud(indm(2)) = apt(indm(2),reinas,gps);
    if min(aptitud) == 0
        break
    end
    
end

[a,b] = min(aptitud);
x = sprintf('El indice %d tiene aptitud %d en %d iteraciones.',b,a,i);
z = sprintf('Número de reinas: %d, tamaño del tablero: %d * %d.',reinas,size,size);
disp(x)
disp(z)
plot(gps(1,:,b)-1/2,gps(2,:,b)-1/2,'d')
title(x)
toc
grid on
xticks(0:size)
yticks(0:size)


%%%%%%%%%%%%%%%%%%%% Funciones trucutru %%%%%%%%%%%%%%%%%%%%%%
function apt = apt(i,reinas,gps)
    apt = 0;
    for j = 1:reinas
        for k = 1:reinas
            if j ~= k
                if abs(gps(1,k,i)-gps(1,j,i)) == abs(gps(2,k,i)-gps(2,j,i)) || abs(gps(1,k,i)-gps(1,j,i)) == 0 || abs(gps(2,k,i)-gps(2,j,i)) == 0
                    apt = apt + 1;
                end
            end           
        end  
    end
end

function hijo = mutacion(hijo,reinas,size)
    p_mutacion = rand();
    if p_mutacion <= 0.8
        g_mutado = randi(reinas);
        x = randi(size);
        y = randi(size);
        for j = 1:reinas       
            while x == hijo(1,j) && y == hijo(2,j)
            x = randi(size);
            y = randi(size);   
            end
        end
        hijo(1,g_mutado)= x;
        hijo(2,g_mutado)=y;
    end
end