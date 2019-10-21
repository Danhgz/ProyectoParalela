/*
Daniel Salazar M. B87214
Ricardo Franco R. B83050
*/
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "Aleatorizador.h"
#include <exception>
#include <omp.h>
#include "Generador.h"
#include <math.h>

using namespace std;

/* Métodos para leer los datos de un archivo de texto
*/
double stod_wrapper(string v) throw (invalid_argument, out_of_range) { return std::stod(v); }
int stoi_wrapper(string v)throw (invalid_argument, out_of_range) { return std::stoi(v); }

template < typename T, class F >
vector< vector< T > > carga_valida_datos(ifstream& archivo, F t) throw (invalid_argument, out_of_range);

void importe();
/***********************************************/


/* Métodos */
void enseneVec();
void muestreDistancias();
void conseguir_centroides();
//void genereCosto(int);
//double genereDistanciaCent(int, int);
void agrupe();
void imprime_centroides();
vector< vector <int> > agrupe_KMeans(vector <int>&);
int busque_pesos(vector <int>&);
void KMeans_Centroides();
void kgrupos(vector<vector <int>>&);
void agregue_centroides(vector < vector <int> >&);
void agrupe_valores(vector < vector <int> >&);
vector <double> reagrupe(vector <int>&);
void genereCosto();

/**********************************/

/*Variables*/
vector < vector <double> > datos;									//vector de los vectores de datos
vector < vector <int> > g_centroides;							//cada vector es un grupo de datos y el primero es su centroide respectivo.
vector < vector <double> > centroides;											//conjunto de los indices de aquellos vectores que son centroides.
vector <vector<double>> distancias_centroides;						//distancias generadas entre cada dato y cada centroide
double costo = 0;

int n;																//dimension: ej si es de 2 es x,y
int m;																//cantidad de vectores presentes con datos
int k;																//cantidad de agrupaciones y de centroides
int numeroHilos;
int epsi;
double c_;
/**********************************************/
Generador generator;

int main(int argc, const char* argv[]) {
	//cout << "Escriba la cantidad de agrupamientos que desea " << endl;
	//cin >> k;
	//cout << "Escriba la cantidad de vectores de datos que desea " << endl;
	//cin >> m;
	//cout << "Escriba la dimensión deseada " << endl;
	//cin >> n;
	k = 100;
	m = 10000;
	n = 100;
	epsi = 500;

	vector < vector <int> > nombre(k, *new vector <int>);							//necesario para luego igualar a g_centroides
	cout << "Numero de hilos: " << endl;
	cin >> numeroHilos;
	numeroHilos = numeroHilos * omp_get_num_procs();
	generator.hilos = numeroHilos;
	importe();
	cout << "Sus datos son los siguientes: " << endl;
	//enseneVec();
	cin.ignore();

	/*****************************************Comienza a paralelizar*******************************************/
	double start = omp_get_wtime();
	conseguir_centroides();
	cout << "Termino centroides" << endl;
	kgrupos(nombre);
	double finish = omp_get_wtime();
	/****************************************Termina de paralelizar*********************************************/
	double tiempo = finish - start;
	g_centroides = nombre;
	cout << "----- Distancias de cada vector de datos con respecto a cada centroide: ------" << endl;
	cout << "Los datos se han agrupado de la siguiente manera: " << endl;
	//muestreDistancias();
	cout << "----- Agrupamientos de los vectores ------" << endl;

	generator.imprime_grupos(g_centroides, datos, centroides, n);
	cout << "El tiempo que duro realizando el trabajo paralelizado fue de: " << tiempo << " segundos." << endl;
	cin.ignore();
	return 0;
}

void conseguir_centroides() {
	srand(time(NULL));
	int pos_aleatoria;
	pos_aleatoria = rand() % m;
	centroides.push_back(datos[pos_aleatoria]);					//un valor aleatorio sera el nuevo centroide. Se guarda su indice
	
	double l = k/2;
	double prob = 0;
	double distancia = 0;
	/* initialize random seed: */
	srand(time(NULL));
	/* inicializa el generador de números entre 0 y 1 */
	Aleatorizador::inicializar_generador_random();	// inicializa aleatorizador
	//for que va a iterar 5 veces buscando conseguir varios centroides nuevos
	for (int it = 0; it < 6; it++)
	{
#pragma omp parallel num_threads(numeroHilos) private(distancia,prob)
	{
		genereCosto();
		//doble for que genera una probabilidad y otra aleatoria. Si la primera > aleatoria entonces el punto respectivo = nuevo centroide.
#pragma omp barrier
#pragma omp  for 
			for (int i = 0; i < m; i++) {
				distancia = 0;
				for (int j = 0; j < centroides.size(); j++) {
					distancia += generator.genereDistanciaCent(i, j, datos, n);
				}
				prob = (l * distancia) / costo;														//probabilidad de cada vector de datos
					if (Aleatorizador::random_uniform_real(Aleatorizador::generador) < prob) {		//si la prob es mayor que el valor aleatorio, se vuelve centroide
#pragma omp critical
						{
							if (find(centroides.begin(), centroides.end(), datos[i]) == centroides.end()) {
								centroides.push_back(datos[i]);
							}
						}
					}
#pragma omp critical
					cout << i << endl;
			}
	}
#pragma omp master
		cout << "Iteracion terminada " << endl;
	}
	cout << "Termino etapa 1 donde solo hay centroides aleatorios" << endl;
	agrupe();			//recordar que centroides trabaja con los indices de datos[][]

}

void agrupe() {		//Entra un solo hilo
	vector<vector<int>> grupos_centroides;									// se añadiran los k centroides y luego aquellos otros centroides cercanos. Tiene los indices de centroides 
	srand(time(NULL));
	int c;
	vector <int> centroides_distanciados;									//vector temporal para meter un centroide y luego guardarlo en el vector de grupos_centroides.
	c = rand() % centroides.size();
	centroides_distanciados.push_back(c);	
	double distancia = 0;
	double max = -1;																		//mayor distancia
	int c_maslargo = 0;																//centroide que posee esa mayor distancia
	bool agregado;
	int counter;
	double prob=0;
	Aleatorizador::inicializar_generador_random();
	while (centroides_distanciados.size() < k) {
		agregado = false;
		counter = 0;																						//counter es el contador que recorre centroides
		c = 0;																								//c va a recorrer centroides_distanciados
		while ( (agregado == false) && (counter < centroides.size() ) ) {
			distancia = generator.genereDistanciaCent(counter, c,centroides,n);								//saca distancias trabajando solo entre centroides
			prob = distancia / costo;
			if((Aleatorizador::random_uniform_real(Aleatorizador::generador)) < prob && find(centroides_distanciados.begin(), centroides_distanciados.end(), counter) == centroides_distanciados.end() ) {
				centroides_distanciados.push_back(counter);										//agrega a centroides_dist el indice del centroide en centroides.
				agregado = true;
			}
			if (c + 1 == centroides_distanciados.size()) {
				c = 0;
			}
			else {
				c++;
			}
			counter++;
		}
	}
		grupos_centroides = agrupe_KMeans(centroides_distanciados);					//Actualizacion. Ahora cambian los centroides.
		distancias_centroides.clear();												//se vacía para posteriormente llenar la matriz con sus verdaderos datos correspondientes.


		vector< vector<double> > nuevos_centroides;
		for (int i = 0; i < grupos_centroides.size(); i++) {
			nuevos_centroides.push_back( centroides[busque_pesos(grupos_centroides[i])] );						
		}
		centroides = nuevos_centroides;												//cambio de centroides completado.
}


int busque_pesos(vector<int> &set_centroides) {							//Solo entra 1 hilo
	vector<int> sumas_pesos(set_centroides.size(), 0);
	double distancia = 0;
	double valor = 1000000;
	int c_indice = 0;
	//doble for que va a buscar cual es el centroide que tiene mayor peso muestral. Compara cada vector de datos[][] con cada posible centroide.
#pragma omp parallel for num_threads(numeroHilos) private(distancia,valor,c_indice) shared(sumas_pesos)
	for (int i = 0; i < datos.size(); i++) {
		distancia = 0;
		valor = 1000000;
		for (int it = 0; it < set_centroides.size(); it++) {
			for (int j = 0; j < n; j++)
			{
				distancia += pow((datos[i][j] - centroides[set_centroides[it]] [j]), 2);					// Sumatoria de (x-x) al cuadrado (distancia)
			}
			distancia = sqrt(distancia);	
			if (valor > distancia){
				valor = distancia;
				c_indice = it;
			}
		}
#pragma omp atomic
		sumas_pesos[c_indice]++;
	}

	int contador = 0;
	valor = 0;
	c_indice = 0;
	//for que va a iterar y buscar el punto que tenga mayor peso. Ese se convertirá en un nuevo centroide.
	for (auto it = set_centroides.begin(); it != set_centroides.end(); it++) {
		if (valor < sumas_pesos[contador]) {
			valor = sumas_pesos[contador];
			c_indice = *it;
		}
		contador++;
	}
	return c_indice;												
}

void KMeans_Centroides() {  /*Entran varios hilos. Implementación de KMeans Paralelo*/
	vector<double> vector_tmp;
	double distancia = 0;
	//doble for para sacar distancias de cada vector de datos con respecto a cada centroide definitivo.
//#pragma omp for private(vector_tmp)
	for (int it = 0; it < datos.size(); it++) {
		for (int it2 = 0; it2 < centroides.size(); it2++) {
			for (int j = 0; j < n; j++)
			{
				distancia += pow((datos[it][j] - centroides[it2][j]), 2);			// Sumatoria de (x-x) al cuadrado (distancia)
			}
			distancia = sqrt(distancia);
			vector_tmp.push_back(distancia);
			distancia = 0;
		}
//#pragma omp critical
		distancias_centroides.push_back(vector_tmp);
		vector_tmp.clear();
	}
//#pragma omp barrier
}

void agregue_centroides(vector<vector <int>>&nombre) {
#pragma omp parallel for num_threads(numeroHilos)
	for (int i = 0; i < centroides.size(); i++) {							//agrega cada centroide a su grupo respectivo, y de primero. Lo que se guarda es el indice del centroide al que pertenece.
		nombre[i].push_back(i);												//centroide pos 0 en la pos 0 de nombre, etc.
	}
#pragma omp barrier
}

void agrupe_valores(vector < vector <int> >& nombre){
		//agrupa los valores dentro del vector nombre
		double minimo = 10000000;
		int centroide = 0;
		bool repite = false;
		int c = 0;
		int cantidad = 0;
#pragma omp parallel for num_threads(numeroHilos) private(centroide, minimo, c, cantidad)
		for (int it = 0; it < distancias_centroides.size(); it++) {
			cantidad = 0;
			c = 0;
			minimo = 10000000;
			for (int i = 0; i < k; i++) {
				if (distancias_centroides[it][i] < minimo) {
					minimo = distancias_centroides[it][i];
					centroide = i;																													//Guarda el indice del centroide
				}
			}
#pragma omp critical
				nombre[centroide].push_back(it);												//agrega los demás puntos al grupo respectivo de cercanía con dicho centroide.	
		}
}

void kgrupos(vector < vector <int> >& nombre) {
	vector< vector <double> > centroides_tmp;
	agregue_centroides(nombre);
	KMeans_Centroides();
	agrupe_valores(nombre);																			//vector de grupacion nombre: listo.
	double distancia = 0;
	double valor = 0;
	bool cambio = false;																			//si pasa toda una iteracion y queda en falso, significa que ya estan bien agrupados.
	int cont = 0;
	double costo_nuevo = 0;
	cout << "empieza while: " << endl;
	while (costo - costo_nuevo >= epsi) {
		costo = costo_nuevo;
		cout << "Costo: " << costo << "Costo Nuevo: " << costo_nuevo << endl;
		costo_nuevo = 0;
		//llamado al metodo de reagrupar: genera nuevos centroides sacando medias
		for (int j = 0; j < nombre.size(); j++) {
			centroides_tmp.push_back(reagrupe(nombre[j]));
		}
		centroides = centroides_tmp;																//elimina los centroides actuales y los reemplaza por los nuevos.
		centroides_tmp.clear();
		nombre.clear();																				//borra las agrupaciones anteriores para crear nuevas
		cambio = true;

		//comienza a conseguir nuevas distancias. Si bool cambio se mantiene en true significa que ya no realiza mas cambios y el agrupamiento está completo.
#pragma omp parallel for num_threads(numeroHilos) private(valor,distancia) shared(cambio,nombre) 		
		for (int it = 0; it < datos.size(); it++) {
			valor = 0;
			for (int it2 = 0; it2 < centroides.size(); it2++) {
				for (int j = 0; j < n; j++)
				{
					valor += pow((datos[it][j] - centroides[it2][j]), 2);					// Sumatoria de (x-x) al cuadrado (distancia)
				}
				valor = sqrt(valor);
				if (distancias_centroides[it][it2] != valor) {
					distancias_centroides[it][it2] = valor;
					cambio = false;
				}
			}
		}

		while (cont < k) {
			nombre.push_back(*new vector<int>);
			cont++;
		}
		//muestreDistancias();
		agregue_centroides(nombre);																//guarda de primeros los centroides en sus agrupaciones respectivas
		agrupe_valores(nombre);																	//compara distancias y agrupa en k grupos.

		for (int i = 0; i < nombre.size(); i++) {
			for (int j = 1; j < nombre[i].size(); j++) {									//empieza en 1 porque la pos 0 es un centroide
				costo_nuevo += generator.recalculeCosto(j, i, datos, centroides, m, n);
			}
		}

		cont = 0;
	}
}


vector< vector <int> >  agrupe_KMeans(vector <int>& centroides_kmeans) {		//Solo entra 1 hilo
	vector< vector <int> > grupos_centroides;															//matriz de int que tendra k columnas y centroides.size() filas.
	vector<int> vector_tmp;
	for (int it2 = 0; it2 < centroides_kmeans.size(); it2++) {
		vector_tmp.clear();
		vector_tmp.push_back(centroides_kmeans[it2]);													//guarda de primero en cada k grupo el indice al centroide respectivo.
		grupos_centroides.push_back(vector_tmp);
	}
	vector<double> distancias_tmp;
	double distancia = 0;

	//En este doble for se van a sacar las distancias de todos los centroides actuales con respecto a los centroides elegidos segun kmeans++.
	//Recordar: En distancias_centroides: filas = vectores de datos, columnas = los centroides.
	for (int it = 0; it < centroides.size(); it++) {
		for (int it2 = 0; it2 < centroides_kmeans.size(); it2++) {
			distancia = generator.genereDistanciaCent(it, centroides_kmeans[it2],centroides,n);
			distancias_tmp.push_back(distancia);
		}
		distancias_centroides.push_back(distancias_tmp);
		distancias_tmp.clear();
	}

	double minimo = 1000000;
	int pos_minima = 0;
	//en este doble for se va a buscar por cada vector de datos, la distancia mas corta hacia un centroide. Posteriormente éste se guarda en grupos_centroides.
#pragma omp parallel for num_threads(numeroHilos) private(pos_minima,minimo)
	for (int i = 0; i < distancias_centroides.size(); i++) {
		minimo = 1000000;
		for (int j = 0; j < k; j++) {									//cada columna es un centroide
			if ((distancias_centroides[i][j] <= minimo)) {
				minimo = distancias_centroides[i][j];
				pos_minima = j;											//guarda el centroide de distancia mas corta
			}	
		}
		if ( pos_minima != centroides_kmeans[pos_minima]) {								//evita repeticion debido a que esos centroides ya están agregados
#pragma omp critical
			grupos_centroides[pos_minima].push_back(i);							//guarda el indice del centroide. Es decir, grupos_centroides esta lleno de indices hacia el vector centroides que tienen las coordenadas.
		}
	}
	return grupos_centroides;
}

vector <double> reagrupe(vector<int> &valores) {																//saca las medias y genera nuevos centroides
	vector<double> centroide_media;								//en este vector se van a guardar los nuevos valores de x,y,z... 
	int j = 0;
	double suma = 0;
	for (int j = 0; j < n;j++) {
		suma = 0;
#pragma omp parallel for num_threads(numeroHilos) reduction(+:suma)
		for (int i = 0; i < valores.size(); i++) {
			if (i == 0) {
				suma += centroides[valores[0]][j];
			}
			else {
				suma += datos[valores[i]][j];
			}
		}
#pragma omp barrier
#pragma omp single
		centroide_media.push_back((suma / valores.size()));
	}
	return centroide_media;
}

//metodo temporal para ver resultados
void imprime_centroides() {
	cout << "Los centroides son: " << endl;
	for (int t = 0; t != centroides.size(); t++) {
		for (int it = 0; it < n; it++) {
			cout<< centroides[t][it] << ", " << endl;
		}
		cout << endl;
	}
	cout << endl;
}

void importe() {
	/* lee el archivo enteros */
	ifstream ee("enteros_erroneo.csv", ios::in);
	if (!ee)
		cout << "no encuentra el archivo de datos" << endl;
	vector< vector< double > > v;
	try {
		v = carga_valida_datos< double >(ee, stoi_wrapper);
	}
	catch (exception ee) {
		cout << "valor invalido o fuera de limite" << endl;
	}
	for (auto f : v) {
		vector<double> tmp;
		for (double x : f) {
			tmp.push_back(x);
			//	cout << x << ", ";
		}
		if (!tmp.empty()) {
			datos.push_back(tmp);
		}
		//cout << endl;
	}
}

void enseneVec() {
	for (int it = 0; it < m; it++) {
		for (auto it2 = 0; it2 < n; it2++) {
			cout << datos[it][it2] << ", ";
		}
		cout << endl << endl;
	}
}


void muestreDistancias() {
	int i = 1;
	cout << "\t" << "\t";
	while (i != centroides.size() + 1) {
		cout << "C#" << i << "\t";
		i++;
	}
	cout << endl;
	for (i = 0; i < datos.size(); i++) {
		cout << "Vector# " << i << "\t";
		for (int j2 = 0; j2 < centroides.size(); j2++) {
			cout << distancias_centroides[i][j2] << "\t";
		}
		cout << endl << endl;
	}
}

void genereCosto() {
	double valor = 0;
	double minimo = 100000;
#pragma omp barrier
#pragma omp for reduction(+:c_) private(valor, minimo)
	for (int i = 0; i < m; i++)
	{
		for (int cent = 0; cent < centroides.size(); cent++) {
			minimo = 100000;
			for (int j = 0; j < n; j++)
			{
				valor = pow((datos[i][j] - centroides[cent][j]), 2);			// Sumatoria de (x-x) al cuadrado (distancia)
				valor = sqrt(valor);
				if (valor < minimo) {
					minimo = valor;
				}
			}
#pragma omp atomic
			c_ += minimo;
		}
	}
#pragma omp single
	costo = c_;
}

template < typename T, class F >
vector< vector< T > > carga_valida_datos(ifstream& archivo, F t) throw (invalid_argument, out_of_range)
{
	vector< vector< T > > valores;
	vector< T > linea_valores;
	string linea;
	while (getline(archivo, linea)) {
		linea_valores.clear();
		stringstream ss(linea);
		string numero_S;
		T numero_T;
		while (getline(ss, numero_S, ',')) {
			try {
				numero_T = t(numero_S);
			}
			catch (exception e) {
				throw e;
			}
			linea_valores.push_back(numero_T);
		}
		valores.push_back(linea_valores);
	}
	return valores;
}