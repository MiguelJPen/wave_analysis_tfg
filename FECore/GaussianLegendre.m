% Returns Gauss Legendre integration points and weights for 2 and three point integration.
% 
% Created:       27 August, 2017
% Last Modified: 23 Nov, 2020
% Author: Abdullah Waseem

% https://es.wikipedia.org/wiki/Cuadratura_de_Gauss
% Nos permite aproximar el resultado de una integral definida de una función
% de grado 2*nnp - 1 o menos.

% nnp: Número de puntos para la curvatura de Gauss-Legendre para integrar
% en el dominio [-1, 1]
% glz: Vector donde se almacenan los puntos.
% glw: Vector donde se almacenan los pesos.

% glw contiene la información donde se ve el valor de los puntos glw del
% polinomio de Legendre.

if nnp == 1
	
	glz(1,1) = 0;
	glw(1,1) = 2;
    
    ngp = 1;
	
elseif nnp == 2
    
    glz(1,1) = -1/sqrt(3);
    glz(2,1) =  1/sqrt(3);
    
    glw(1,1) =  1;
    glw(2,1) =  1;

    ngp = 2;

elseif nnp == 3
        
    glz(1,1) = -sqrt(3/5);
    glz(2,1) =  0;
    glz(3,1) =  sqrt(3/5);
    
    glw(1,1) =  5/9;
    glw(2,1) =  8/9;
    glw(3,1) =  5/9;

    ngp = 3;
    
elseif nnp == 4

    glz(1,1) = -sqrt(3/7 + 2/7 * sqrt(6/5));
    glz(2,1) = -sqrt(3/7 - 2/7 * sqrt(6/5));
    glz(3,1) = 0;
    glz(4,1) = sqrt(3/7 - 2/7 * sqrt(6/5));
    glz(5,1) = sqrt(3/7 + 2/7 * sqrt(6/5));

    glw(1,1) = (322 - 13 * sqrt(70))/900;
    glw(2,1) = (322 + 13 * sqrt(70))/900;
    glw(3,1) = 128/225;
    glw(4,1) = (322 + 13 * sqrt(70))/900;
    glw(5,1) = (322 - 13 * sqrt(70))/900;

    ngp = 5;
end