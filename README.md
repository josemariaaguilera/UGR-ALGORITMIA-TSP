## Algorítmo para resolver el problema del vendedor (TSP)

###Instalación y ejecución
el programa está escrito en c++ por lo que tendrás que tener el compilador g++ para poder compilarlo y ejecutar el programa.


Primero tenemos que clonar el repositorio en nuestro dispositivo, se puede hacer usando la orden:   
``git clone https://github.com/josemariaaguilera/UGR-ALGORITMIA-TSP.git ``

####Windows
Para compilar el programa:  
``g++ -std=c++11 -O3 tsp2.1.cpp -o tsp.exe``  
Para ejecutar el programa:  
``./tsp.exe ./datos/ulysses<n>.tsp``

####Linux
Para compilar el programa:  
``g++ -std=c++11 -O3 tsp2.1.cpp -o tsp``  
Para ejecutar el programa:  
``./tsp ./datos/ulysses<n>.tsp``

###Problema
Antes de nada, hay que describir el problema a resolver.
Dado un numero de ciudades con su posición, trata de calcular el camino mas corto para ir todas las ciudades, volviendo a la de partida, sin pasar dos veces por una ciudad.

###Descripción del algoritmo
- genera una matriz con todas las distancias entre todas las ciudades.

- ya que es un ciclo, da igual en que ciudad se empiece, asi que añade la
primera ciudad a la solución.

- genera una solucion no optima escogiendo como siguiente ciudad la mas
cercana a la ultima añadida. esto nos servirá para eliminar todas las opciones
que sean peores que esta.

- calcula las distancias mas cortas entre cada ciudad y cualquier otra. esto
servirá para crear una estimación optimista para completar el camino. si
esta estimación es peor que la mejor solucion que tenemos hasta el momento
(cotasuperior) este camino y todos sus derivados se podran descartar.

- una vez hecho esto, hasta que las posibilidades no se acaben:
  - comprueba que la estimación optimista del camino no es peor que la mejor
  solucion que tenemos hasta el momento, sino la descarta y a todas sus variantes.

  - genera todas las variantes posibles (sin repetir ciudad) del camino que
  evaluamos.

  - calcula la estimación optimista para el camino que evaluamos (cotalocal).

  - calcula un valor haciendo la media ponderada entre la estimación optimista y
  la mejor solucion hasta el momento (el camino en el que este valor sea menor
  sera el siguiente en explorarse).

  - si el camino ya paso por todas las ciudades, comprueba si es mejor que el
  camino mas corto que hemos descubierto hasta el momento, si es asi, este se
  convierte en el mejor camino.

###Datos de interes
Para comprender un poco mejor el código, te podría ser de ayuda ver la imagen que esta en la carpeta /diagrama, pero antes, es recomendable leer la descripción del algoritmo.

###Autores
- José María Aguilera Barea (@josemariaaguilera)
- Matilde Cabrera Gonzalez (@mati3)
