// ProyectoK-Means.cpp : Este archivo contiene la función "main". La ejecución del programa comienza y termina ahí.

#include <omp.h>
#include <conio.h>
#include <fstream>
#include <string>
#include <utility>
#include <random>
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

/*
*@Función: se encarga de leer el archivo .txt con la lista de palabras y guardarlas en un array de tipo Palabra
*@Param: recibe el archivo .txt
*/
void cargarVectores(char* nombreArchivo, int n, int m, vector<vector<double>>& matriz) {
	ifstream archivo(nombreArchivo);
	if (archivo.is_open()) {
		int i = 0, j = 0;
		string numero = "";
		char aux;
		while (!archivo.eof()) {
			aux = archivo.get();
			if (aux != ',' && aux != '\n') {
				numero += aux;
			}
			else if (aux == ',' || aux == '\n') {
				matriz[i][j] = stod(numero);
				numero = "";
				++j;
			}
			if (j == n) {
				j = 0;
				++i;
			}
		}
		archivo.close();
	}
}

double generarDouble(int min, int max) {
	double rand = 0.0;
	default_random_engine generator;
	uniform_real_distribution<double> distribution(min, max);
	rand = distribution(generator);
	return rand;
}

int generarInt(int min, int max) {
	int rand = 0;
	default_random_engine generator;
	uniform_int_distribution<int> distribution(min, max);
	rand = distribution(generator);
	return rand;
}

double calcularDistanciaEuclidiana(vector<double>& vec1, vector<double>& vec2) {
	double distancia = 0.0;
	int tam = vec1.size();
	for (int i = 0; i < tam; ++i) {
		distancia += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
	}
	return distancia;
}
double calcDisMin(vector<double>& puntoEnX, vector<vector<double>>& vectorCentroides) {
	double min = 0.0;
	double aux = 0.0;
	min = calcularDistanciaEuclidiana(puntoEnX, vectorCentroides[0]);
	if (vectorCentroides.size() > 1) {
		int tam = vectorCentroides.size() - 1;
		for (int i = 1; i <= tam; ++i) {
			aux = calcularDistanciaEuclidiana(puntoEnX, vectorCentroides[i]);
			if (aux < min) {
				min = aux;
			}
		}
	}
	return min;
}

int calcMinPos(vector<double>& puntoEnX, vector<vector<double>>& vectorCentroides) { //Posicion donde esta la distancia minima
	int posMin = 0;
	double aux = 0.0;
	double min = calcularDistanciaEuclidiana(puntoEnX, vectorCentroides[0]);
	if (vectorCentroides.size() > 1) {
		int tam = vectorCentroides.size() - 1;
		for (int i = 1; i <= tam; ++i) {
			aux = calcularDistanciaEuclidiana(puntoEnX, vectorCentroides[i]);
			if (aux < min) {
				min = aux;
				posMin = i;
			}
		}
	}
	return posMin;
}

double calcularProb(vector<double>& puntoEnX, vector<vector<double>>& vectorCentroides, double psi) {
	double prob = 0.0;
	double min = calcDisMin(puntoEnX, vectorCentroides);
	prob = (min / psi);
	return prob;
}

//calcular probabilidades de cada punto
void probabilidadesPuntos(vector<vector<double>>& vectorDatos, vector<vector<double>>& vectorCentroides, vector<double>& probabilidades, double psi, double l) {
	double prob = 0.0;
	int tamano = vectorDatos.size();
	for (int i = 0; i < tamano; ++i) {
		prob = l * calcularProb(vectorDatos[i], vectorCentroides, psi);
		probabilidades.push_back(prob);
	}
}

//****************************************************************************************
//calcular en probabilidad
//recalcular costo de puntos con cada centroide, debe ser el minimo. d(x,C) distancia de x al centroide mas cercano en C

int compararVectores(vector<double> vec1, vector<double> vec2) {
	int iguales = 1;
	int tam = vec1.size();
	for (int i = 0; i < tam && iguales; ++i) {
		if (vec1[i] != vec2[i]) {
			iguales = 0;
		}
	}
	return iguales;
}
int validarCentroide(vector<vector<double>>& centroides, vector<double>& unX) { //ver si el centroide no existe
	int valido = 1;
	int tam = centroides.size();
	for (int i = 0; i < tam ; ++i) {
		if (compararVectores(centroides[i], unX)) { // Entra si son iguales
			valido = 0;
		}
	}
	return valido;
}

void reclusterKmeans2(vector<vector<double>>& vectorCentroides, vector<vector<double>>& centroidesFinales, int k) {
	int posR = 0;
	posR = generarInt(0, vectorCentroides.size() - 1);
	centroidesFinales.push_back(vectorCentroides[posR]);
	int contador = 1; //empieza en 1 porque ya se inserto un centroide
	double fCosto = 0.0;

	int tamCentroides = vectorCentroides.size();
	for (int i = 0; i < tamCentroides; ++i) {
		fCosto += calcDisMin(vectorCentroides[i], centroidesFinales);
	}
	vector<double> probabilidades;
	double prob = 0.0;
	while (contador < k) {
		probabilidadesPuntos(vectorCentroides, centroidesFinales,probabilidades, fCosto, 1);
		posR = generarInt(0, vectorCentroides.size() - 1);
		prob = generarDouble(0, 1);
		if (prob < probabilidades[posR]) {
			if (validarCentroide(centroidesFinales, vectorCentroides[posR])) {
				centroidesFinales.push_back(vectorCentroides[posR]);
				++contador;
				fCosto = 0;
				for (int i = 0; i < tamCentroides; ++i) {
					fCosto += calcDisMin(vectorCentroides[i], centroidesFinales);
				}
			}
		}
	}
}
void kmeansParallelInit(vector<vector<double>>& vectorDatos, vector<vector<double>>& vectorCentroides, vector<vector<vector<double>>>& valoresCentroides, int hilos, int m, int k) {
	vector<double> probabilidades; //las probabilidades de cada x en X de ser escogido como centroide
	vector<vector<double>> vCentroidesPrime;
	double l = k / 2;
	double fCosto = 0.0;
	int posicionR = generarInt(0, m - 1);
	vectorCentroides.push_back(vectorDatos[posicionR]);

	//Para calcular el costo inicial con el centroide inicial aleatorio
	for (int i = 0; i < m; ++i) {
		fCosto += calcDisMin(vectorDatos[i], vectorCentroides);
	}

	int contador = 0;
	double escoger = 0.0; // para escoger el numero random entre 0-1 y elegir un centroide si es menor a una prob
	for (int i = 0; i < 5; ++i) {
		probabilidadesPuntos(vectorDatos, vectorCentroides, probabilidades, fCosto, l);
# pragma omp parallel num_threads(hilos) shared(contador) private(escoger, posicionR)
		{
			while (contador < l) {
				escoger = generarDouble(0, 1);
				posicionR = generarInt(0, m - 1);
				if (escoger < probabilidades[posicionR]) {
					if (validarCentroide(vectorCentroides, vectorDatos[posicionR])) {
# pragma omp critical
						{
							vectorCentroides.push_back(vectorDatos[posicionR]);
							++contador;
						}
					}
				}
			}
			fCosto = 0;
# pragma omp single //Critical?
			{
				for (int i = 0; i < m; ++i) {
					fCosto += calcDisMin(vectorDatos[i], vectorCentroides);
				}
			}
			contador = 0;
		}
	}
	//En teoria aqui ya se tienen los centroides iniciales, ahora hay que calcular pesos por centroide y luego recluster
	/*
	int posMin = 0;
# pragma omp parallel for num_threads(hilos)
	for (int i = 0; i < m; ++i) {
# pragma omp critical
		{
			posMin = calcMinPos(vectorDatos[i], vectorCentroides);
			vector<double> lol = vectorDatos[i];
			valoresCentroides[posMin].push_back(lol);
		}
	}*/
	vector<vector<double>> centroidesFinales;
	reclusterKmeans2(vectorCentroides,centroidesFinales, k);
	vectorCentroides.assign(centroidesFinales.begin(), centroidesFinales.end());
	valoresCentroides.resize(vectorCentroides.size());
}
void promediarGrupo(vector<vector<double>>& vectorGrupo, vector<double>& vectorPromedio) {
	int tamGrupo = vectorGrupo.size();
	int tamDato = vectorGrupo[0].size();
	vectorPromedio.assign((int)vectorPromedio.size(), 0);
	for (int i = 0; i < tamGrupo; ++i) 
	{
		for (int j = 0; i < tamDato; ++i) {
			vectorPromedio[j] += vectorGrupo[i][j];
		}
	}
	for (int i = 0; i < tamDato; ++i) {
		vectorPromedio[i] /= tamDato;
	}
}

double algoritmoLloyd(vector<vector<double>>& vectorDatos, vector<vector<double>>& vectorCentroides, vector<vector<vector<double>>>& puntosAsociadosC, int m, int k, double eps) {
	int posMin = 0;
	double fCosto = 0.0;
	for (int i = 0; i < m; ++i) {
		fCosto += calcDisMin(vectorDatos[i], vectorCentroides);
		posMin = calcMinPos(vectorDatos[i], vectorCentroides);
		puntosAsociadosC[posMin].push_back(vectorDatos[i]);
	}
	double fCostoPrime = fCosto;
	do {
		puntosAsociadosC.clear();
		puntosAsociadosC.resize(k);
		fCosto = fCostoPrime;
		for (int i = 0; i < puntosAsociadosC[i].size(); ++i) {
			promediarGrupo(puntosAsociadosC[i], vectorCentroides[i]);
		}
		fCostoPrime = 0;
		for (int i = 0; i < m; ++i) {
			fCostoPrime += calcDisMin(vectorDatos[i], vectorCentroides);
			posMin = calcMinPos(vectorDatos[i], vectorCentroides);
			puntosAsociadosC[posMin].push_back(vectorDatos[i]);
		}
	} while ((fCosto - fCostoPrime) > eps);
	return fCostoPrime;
}

void escribirResultados(vector<vector<vector<double>>>& vGrupos, vector<vector<double>>& centroides, double fCosto, double tPared) {
	int cantGrupos = vGrupos.size();
	int n = centroides.size();
	ofstream salida;
	salida.open("salida.csv");
	for (int i = 0; i < cantGrupos; ++i) 
	{
		salida << "Centroide:" << endl;
		for (int j = 0; j < n; ++j) {
			salida << centroides[i][j];
			if (j != n - 1) {
				salida << ", ";
			}
		}
		salida << endl<< "Datos:"<<endl;
		for (int j = 0; j < vGrupos[i].size(); ++j) 
		{
			for (int k = 0; k < n; ++k) {
				salida << vGrupos[i][j][k];
				if (k != n - 1) {
					salida << ", ";
				}
			}
			salida << endl;
		}
		cout << "\n\n\n";
	}
	salida << "Tiempo pared: "<< tPared <<"   Funcion de Costo: " << fCosto;
	salida.close();
}
void imprimirVectores(int n, int m, vector<vector<double>>& datos) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << datos[i][j] << "\t";
		}
		cout << endl;
	}
}
/*
arg[1] = archivo de entrada
arg[2] = n (dimension datos)
arg[3] = m (cantidad vectores)
arg[4] = k (cantidad centroides)
arg[5] = hilos
*/
int main(int argc, char** argv) {
	if (argc == 5) {
		int n = atoi(argv[2]);
		int m = atoi(argv[3]);
		int k = atoi(argv[4]);
		int hilos = atoi(argv[5]);
		int procesadores = omp_get_num_procs();
		hilos *= procesadores;
		double eps = 1.0;
		double fCosto = 0.0;
		vector<vector<double>> datos(m, vector<double>(n));
		vector<vector<double>> centroide;
		vector<vector<vector<double>>> valoresAsociadosC (k);
		cargarVectores(argv[1], n, m, datos);
		imprimirVectores(n, m, datos);
		auto ini = omp_get_wtime();
		kmeansParallelInit(datos, centroide, valoresAsociadosC, hilos, m, k);
		fCosto = algoritmoLloyd(datos, centroide, valoresAsociadosC,  m, k, eps);
		auto fin = omp_get_wtime();
		escribirResultados(valoresAsociadosC,centroide, fCosto, (double)(fin-ini));
	}
	return 0;
}