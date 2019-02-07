// g++ -std=c++11 -O3 tsp2.1.cpp -o tsp
// ./tsp ./datos/ulysses<n>.tsp
/**
* el programa lee de un fichero una lista de ciudades con sus coordenadas y
* calcula el camino mas corto para ir todas las ciudades, volviendo a la de
* partida, sin pasar dos veces por una ciudad. para ello el algoritmo:
*
* - genera una matriz con todas las distancias entre todas las ciudades.
*
* - ya que es un ciclo, da igual en que ciudad se empiece, asi que añade la
* primera ciudad a la solución.
*
* - genera una solucion no optima escogiendo como siguiente ciudad la mas
* cercana a la ultima añadida. esto nos servirá para eliminar todas las opciones
* que sean peores que esta.
*
* - calcula las distancias mas cortas entre cada ciudad y cualquier otra. esto
* servirá para crear una estimación optimista para completar el camino. si
* esta estimación es peor que la mejor solucion que tenemos hasta el momento
* (cotasuperior) este camino y todos sus derivados se podran descartar.
*
* - una vez hecho esto, hasta que las posibilidades no se acaben:
*
* - comprueba que la estimación optimista del camino no es peor que la mejor
* solucion que tenemos hasta el momento, sino la descarta y a todas sus variantes.
*
* - genera todas las variantes posibles (sin repetir ciudad) del camino que
* evaluamos.
*
* - calcula la estimación optimista para el camino que evaluamos (cotalocal).
*
* - calcula un valor haciendo la media ponderada entre la estimación optimista y
* la mejor solucion hasta el momento (el camino en el que este valor sea menor
* sera el siguiente en explorarse).
*
* - si el camino ya paso por todas las ciudades, comprueba si es mejor que el
* camino mas corto que hemos descubierto hasta el momento, si es asi, este se
* convierte en el mejor camino.
*/

#include <cstdlib>
#include <iostream>
#include <list>
#include <vector>
#include <chrono>
#include <queue>
#include <map>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

///contendra el numero de ciudades.
int n_ciudades ;
///contendra la mejor solución.
vector<int> solCiudades;
///contendra la posicion de las ciudades identificadas por un entero.
map<int, pair<double, double>> M;
///contendra las distancias mas pequeñas entre una ciudad y su mas cercana.
vector<double> menores;

/**
* La estructura nodo sirve para representar un nodo en el arbol de busqueda.
* consta de:
* cotaLocal => la distancia que se ha recorrido recorriendo las ciudades
* determinadas en la variable solucion.
* estimador => la estimación de la distancia que se prevee como minimo tendra
* si se sigue profundizando hasta la solucion final por este nodo.
* solucion => almacena el identificador de las ciudades que se han recorrido.
*/
struct Nodo{
	double cotalocal=0 ;
	int estimador=0 ;
	vector<int> solucion ;
};

/**
* operador < de nodos que determina que un nodo es menor que otro comparando sus estimadores.
*/
bool operator<(const Nodo & n1, const Nodo & n2){
	return (n1.estimador > n2.estimador);
}

////////////////////////////////////////////////////////////////////////////////
///clase para obtener datos de la linea de comandos proporcionada por el profesor
class CLParser{
	public:
		CLParser(int argc_, char * argv_[],bool switches_on_=false);
		~CLParser(){}
		string get_arg(int i);
		string get_arg(string s);
	private:
		int argc;
		vector<string> argv;
		bool switches_on;
		map<string,string> switch_map;
};

CLParser::CLParser(int argc_, char * argv_[],bool switches_on_){
	argc = argc_;
	argv.resize(argc);
	copy(argv_, argv_ + argc, argv.begin());
	switches_on = switches_on_;
	if(switches_on){
		vector<string>::iterator it1, it2;
		it1 = argv.begin();
		it2 = it1 + 1;
		while(true){
			if(it1 == argv.end())
				break;
			if(it2 == argv.end())
				break;
			if((*it1)[0] == '-')
				switch_map[*it1]=*(it2);
			it1++;
			it2++;
		}
	}
}

string CLParser::get_arg(int i){
	if(i >= 0 && i < argc)
		return argv[i];
	return "";
}

string CLParser::get_arg(string s){
	if(!switches_on)
		return "";
	if(switch_map.find(s) != switch_map.end())
		return switch_map[s];
	return "";
}
///////////////////////////////////////////////////////////////////////////////

/**
* Lee las ciudades de un fichero y las almacena en el mapa<identificador,<x,y>>.
* @param nombre Nombre del fichero que contiene la lista de ciudades con su posición.
* @param M Mapa donde se van a almacenar las ciudades con su posicion.
*/
void leer_puntos(string& nombre, map<int,pair<double,double>>& M){
	ifstream datos;
	string s;
	pair<double,double> p;
	int n, act;
	datos.open(nombre.c_str());
	if(datos.is_open()){
		datos >> s;
		datos >> n;
		n_ciudades = n;
		int i = 0;
		while(i < n){
			datos >> act >> p.first >> p.second;
			M[act] = p;
			i++;
		}
		datos.close();
	}else{
		cout << "Error de Lectura en " << nombre << endl;
	}
 }

 /**
 * Calcula la distancia entre 2 puntos dados(ciudades).
 * @param d1 coordenadas del punto 1.
 * @param d2 coordenadas del punto 2.
 * @return distancia entre las dos posiciones.
 */
double calcularDistancia(pair<double,double>& d1, pair<double,double>& d2){
	double distX = d1.first - d2.first;
	double distY = d1.second - d2.second;
	return (double)(sqrt((distX * distX) + (distY * distY)));
}


/**
* Elimina todas columnas de la matriz de longitudes entre ciudades (las pone a -1).
* @param m matriz de longitudes entre ciudades.
* @param ciudad ciudad de la cual se quiere eliminar sus medidas de la matriz.
*/
void anularColumnas(vector< vector<double> >& m, int ciudad){
	for(int i = 0; i < n_ciudades; i++){
		m[i][ciudad-1] = -1;
	}
}

/**
* Rellena la matriz de longitudes entre ciudades a partir del mapa de ciudades.
* @param longitudes matriz de longitudes que se va a rellenar.
* @param ciudades mapa de ciudades almacenadas en formato <id,coordenadas>.
*/
void rellenarMatriz(vector< vector<double> >& longitudes, map<int,pair<double,double>>& ciudades){
	for(int i = 0; i < n_ciudades; i++){
		for(int j = 0; j < n_ciudades; j++){
			if(i != j){
				longitudes[i].push_back(calcularDistancia(ciudades[i+1], ciudades[j+1]));
			}else{
				longitudes[i].push_back(-1);
			}
		}
	}
}

/**
* Devuelve la siguiente ciudad que está mas cerca de la ciudad de partida,
* elimina las columnas de la matriz de distancias porque esa ciudad ya se ha
* usado y altualiza la distancia desde la ciudad de partida y la ciudad devuelta.
* @param ciudadPartida ciudad de partida para encontrar la siguiente ciudad mas cercana.
* @param m matriz de longitudes entre ciudades calculada.
* @param contDistancia variable donde se almacena la distancia total de
* la posible solucion.
* @return la siguiente ciudad que está mas cerca de la ciudad de partida.
*/
int Greedy(int ciudadPartida, vector< vector<double> >& m, double& contDistancia){
	int ciudadSiguiente = -1;
	double distanciaMinima = 99999999;
	for(int i = 0; i < n_ciudades; i++){
		if(m[ciudadPartida-1][i] < distanciaMinima && m[ciudadPartida-1][i] != -1){
			ciudadSiguiente = i+1;
			distanciaMinima = m[ciudadPartida-1][i];
		}
	}
	if(ciudadSiguiente != -1){
		contDistancia += distanciaMinima;
		solCiudades.push_back(ciudadSiguiente);
		anularColumnas(m, ciudadPartida);
	}
	return ciudadSiguiente;
}

/**
* calcula una posible buena solucion y devuelve la distancia total de esta.
* para poder eliminar los nodos que la superen y no seguir profundizando por ese
* camino para encontrar la solución.
* @param matriz matriz de distancias entre ciudades ya calculada.
* @return la distancia total de una posible buena solución.
*/
double calcularCamino(vector<vector<double> > matriz){
	solCiudades.clear();
	int inicio = 1, fin = -1;
	int siguiente = inicio;
	double contDistancia=0;
	solCiudades.push_back(1);
	for(int i = 0; i < n_ciudades; i++){
		siguiente = Greedy(siguiente, matriz, contDistancia);
	}
	contDistancia += calcularDistancia(M[inicio], M[solCiudades.back()]);
	return contDistancia;
}

/**
* comprueba si en la lista de ciudades del nodo no hay ciudades repetidas.
* @param y nodo del cual se quiere calcular si no hay ciudades repetidas.
* @return true si no ha ciudades repetidas, false de lo contrario.
*/
bool ciudadesNoRepetidas(Nodo & y){
	bool bol=true;
	int fin=y.solucion.back();
	for(int i=0; i<y.solucion.size()-1 && bol; i++){
		if(y.solucion.at(i)==fin){
			bol=false;
		}
	}
	return bol;
}

/**
* Si el nodo es factible quiere decir que  en la lista de ciudades del nodo no
* hay ciudades repetidas.
* @param y nodo del cual se quiere calcular si no hay ciudades repetidas.
* @return true si no ha ciudades repetidas, false de lo contrario.
*/
bool Factible(Nodo & y){
	return ciudadesNoRepetidas(y);
}

/**
* comprueba si el nodo tiene todas las ciudades utilizadas.
* @param y nodo del cual se quiere comprobar si es solucion.
* @return true si el numero de ciudades del nodo es el numero total de ciudades,
* false en caso contrario.
*/
bool esSolucion(Nodo & nodo){
	return nodo.solucion.size() >= n_ciudades;
}

/**
* devuelve el minimo de 2 números.
* @param num1 numero 1.
* @param num2 numero 2.
* @return minimo de 2 numeros.
*/
int Min(int num1, int num2){
	if(num1 < num2){
		return num1 ;
	}else
		return num2 ;
}

/**
* devuelve la distancia entre la ciudad y la mas proxima.
* @param ciudad id de ciudad elegida.
* @param m vector de distancias entre ciudades.
* @return la distancia entre la ciudad y la mas proxima.
*/
double menorArista(int ciudad, const vector<vector<double> > &m){
	double distancia=9999999;
		for(int i=0; i<n_ciudades; i++){
			if(m.at(ciudad-1).at(i)<distancia && i!=ciudad-1){
				distancia=m.at(ciudad-1).at(i);
			}
		}
	return distancia;
}

/**
* rellena el vector global de distancias minimas entre ciudades.
* @param m matriz de distancias entre ciudades rellena.
*/
void calcularMenoresAristas(const vector<vector<double> > &m){
	for(int i=1; i<=n_ciudades; i++){
		menores.push_back(menorArista(i,m));
	}
}

/**
* devuelve una estimación de la distancia que puede ser tan buena, o mas (pero no posible),
* que la solucion optima.
* @param sol lista de ciudades de la solucion.
* @param m matriz de distancias entre las ciudades.
* @return una estimación de la distancia que puede ser tan buena, o mas (pero no posible),
* que la solucion optima.
*/
double greedyLocal2(const vector<int> &sol, const vector<vector<double> > &m){
	double estimacionOptima= 0;
	for(int j=0 ; j < n_ciudades ; j++){
		if (find(sol.begin(), sol.end(), j+1) == sol.end()){//para cada una que no este en sol

			double distancia=9999999;
				for(int o=0; o<n_ciudades; o++){
				 //comprueba que la fila no este en la solucion
					if (find(sol.begin(), sol.end(), o+1) == sol.end()){
					 // si no es -1 y la distancia de la posicion es menor a la que ya habia
						if(m.at(j).at(o)<distancia && o!=j){
							distancia=m.at(j).at(o);
						}
					}
				}
				if(distancia==9999999){estimacionOptima+=menores.at(j);}else{
				estimacionOptima+=distancia;}
		}

	}
	estimacionOptima+=menores.at(sol.back()-1);//mas la ultima

	return estimacionOptima;
}

/**
* calcula la distancia total a recorrer para visitar las ciudades de un nodo.
* @param solnodo lista de ciudades del nodo.
* @param m matriz de distancias entre las ciudades.
* @return la distancia total a recorrer para visitar las ciudades de un nodo.
*/
double Coste(vector<int> &solnodo, vector< vector<double> > &m){
	double coste=0;
	for(int i=0; i < solnodo.size()-1; i++){
		coste+=m[solnodo.at(i)-1][solnodo.at(i+1)-1];
	}
	return coste;
}

/**
* calcula el camino mas corto para recorrer todas las ciudades, volviendo a la
* ciudad origen, sin pasar dos veces por una misma ciudad.
* @param m matriz de distancias entre las ciudades.
* @return vector de identificadores de las ciudades en orden, para recorrer la
* menor distancia.
*/
vector<int> TSP (vector<vector<double> > & m){
	vector<int> optimo;

	priority_queue<Nodo> LNV ; // lista de nodos vivos.
	Nodo inicio, x ; // nodo inicial, no tiene cota local ni estimador
	inicio.cotalocal=0;
	inicio.estimador=0;
	optimo.push_back(1);
	inicio.solucion=optimo;
	LNV.push(inicio);

	double cotasup=calcularCamino(m);
	calcularMenoresAristas(m);

	while(!LNV.empty()){
		x = LNV.top();
		LNV.pop();

		if(x.cotalocal < cotasup){

			for(int i=1 ; i < n_ciudades; i++){

				Nodo y;
				y.solucion=x.solucion;
				y.solucion.push_back(i+1);


				if(Factible(y)){
					double coste=Coste(y.solucion,m);
					y.cotalocal=greedyLocal2(y.solucion,m)+coste;
					y.estimador=(y.cotalocal+cotasup)*0.5;

					if(esSolucion(y)){
						/*si optimo solo tiene el nodo raiz lo añades, en caso contrario,
						se comprueba el mejor*/
						if(optimo.size()==1){
							optimo= y.solucion ;
							cotasup = coste;
							cotasup +=m[y.solucion.back()-1][y.solucion.at(0)-1];
						}else{
							coste+=m[y.solucion.back()-1][y.solucion.at(0)-1];
							if(coste < cotasup){
								optimo= y.solucion ;
								cotasup = coste ;
							}
						}//si la estimacion optimista es peor que la mejor solucion que tenemos.
					}else if(y.cotalocal < cotasup){ // estrategia de poda
						LNV.push(y);
					}
				}
			}
		}
	}
	return optimo;
}


///////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------

int main(int argc, char * argv[]){
	string filename;
	//tiempos
  double total = 0;
  high_resolution_clock::time_point tantes,tdespues;
  duration<double> tiempo_transcurrido;

	if(argc < 2){
		cout << "Error Formato: ./tsp puntos orden_correcto" << endl;
		exit(1);
	}

	CLParser cmd_line(argc, argv, false);
	filename = cmd_line.get_arg(1);
	leer_puntos(filename,M);

	cout <<"n_ciudades: "<< n_ciudades << endl;

	//rellena la matriz de ciudades
	vector< vector<double> > matriz(n_ciudades);
	rellenarMatriz(matriz, M);

	double contDistancia = 0;
	tantes = high_resolution_clock::now();

	//llamar al algoritmo
	vector<int> sol=TSP(matriz);

	//medir tiempos
	tdespues = high_resolution_clock::now();
	tiempo_transcurrido  = duration_cast<duration<double>>(tdespues-tantes);
	total = tiempo_transcurrido.count();

	cout << "tiempo: " << total << endl;

	//muestra resultado
	int ciudad;
	for(int i=0; i<sol.size(); i++){
		ciudad = sol.at(i);
		cout << ciudad << "\t" << M[ciudad].first << "\t" << M[ciudad].second << endl;
	}
	contDistancia=Coste(sol,matriz);
	contDistancia+=matriz[sol.back()-1][sol.at(0)-1];

	cout << endl << "Distancia total: " << contDistancia << endl;

	return 0;
}
