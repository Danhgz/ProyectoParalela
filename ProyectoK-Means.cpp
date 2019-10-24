#include <omp.h>
#include <conio.h>
#include <time.h>
#include <fstream>
#include <string>
#include <utility>
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
	double random = ((double)rand()) / (double)RAND_MAX;
	double diff = (double)max - (double)min;
	double r = random * diff;
	return (double)min + r;
}
int generarInt(int min, int max) {
	int randm = 0;
	randm = rand() % (max - min)+min;
	return randm;
}

double calcularDistanciaEuclidiana(vector<double>& vec1, vector<double>& vec2) {
	double distancia = 0.0;
	int tam = vec1.size();
	for (int i = 0; i < tam; ++i) {
		distancia += abs(vec1[i] - vec2[i]);
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
void probabilidadesPuntos(vector<vector<double>>& vectorDatos, vector<vector<double>>& vectorCentroides, vector<double>& probabilidades, double psi, int hilos, double l) {
	double prob = 0.0;
	int tamano = vectorDatos.size();
	#pragma omp parallel for num_threads(hilos) shared(l)
	for (int i = 0; i < tamano; ++i) {
		prob = l * calcularProb(vectorDatos[i], vectorCentroides, psi);
		probabilidades[i] = prob;
	}
}

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

void reclusterKmeans2(vector<vector<double>>& vectorCentroides, vector<vector<double>>& centroidesFinales, int hilos, int k) {
	int posR = 0;
	posR = generarInt(0, vectorCentroides.size() - 1);
	centroidesFinales.push_back(vectorCentroides[posR]);
	int contador = 1; //empieza en 1 porque ya se inserto un centroide
	double fCosto = 0.0;

	int tamCentroides = vectorCentroides.size();
	#pragma omp parallel for num_threads(hilos) reduction(+:fCosto)
	for (int i = 0; i < tamCentroides; ++i) {
		fCosto += calcDisMin(vectorCentroides[i], centroidesFinales);
	}
	vector<double> probabilidades;
	double prob = 0.0;
	while (contador < k) {
		probabilidades.clear();
		probabilidades.resize(tamCentroides);
		probabilidadesPuntos(vectorCentroides, centroidesFinales,probabilidades, fCosto, hilos, 1);
		posR = generarInt(0, vectorCentroides.size() - 1);
		prob = generarDouble(0, 1);
		if (prob < probabilidades[posR]) {
			if (validarCentroide(centroidesFinales, vectorCentroides[posR])) {
				centroidesFinales.push_back(vectorCentroides[posR]);
				++contador;
				fCosto = 0;
				#pragma omp parallel for num_threads(hilos) reduction(+:fCosto)
				for (int i = 0; i < tamCentroides; ++i) {
					fCosto += calcDisMin(vectorCentroides[i], centroidesFinales);
				}
			}
		}
	}
}

void kmeansParallelInit(vector<vector<double>>& vectorDatos, vector<vector<double>>& vectorCentroides, vector<vector<vector<double>>>& valoresCentroides, int hilos, int m, int k) {
	vector<double> probabilidades; //las probabilidades de cada x en X de ser escogido como centroide
	double l = (double)k / 2;
	double fCosto = 0.0;
	int posicionR = generarInt(0, m - 1);
	vectorCentroides.push_back(vectorDatos[posicionR]);

	//Para calcular el costo inicial con el centroide inicial aleatorio
	#pragma omp parallel for num_threads(hilos) reduction(+:fCosto)
	for (int i = 0; i < m; ++i) {
		fCosto += calcDisMin(vectorDatos[i], vectorCentroides);
	}
	int contador = 0;
	double escoger = 0.0; // para escoger el numero random entre 0-1 y elegir un centroide si es menor a una prob
	for (int i = 0; i < 5; ++i) 
	{
		probabilidades.clear();
		probabilidades.resize(m);
		probabilidadesPuntos(vectorDatos, vectorCentroides, probabilidades, fCosto, hilos, l);
		# pragma omp parallel num_threads(hilos) shared(contador, vectorCentroides) private(escoger, posicionR)
		{
			while (contador < l) {
				escoger = generarDouble(0, 1);
				posicionR = generarInt(0, m - 1);
				if (escoger < probabilidades[posicionR]) {
					if (validarCentroide(vectorCentroides, vectorDatos[posicionR])) {
						#pragma omp critical
						{
							vectorCentroides.push_back(vectorDatos[posicionR]);
							++contador;
						}
					}
				}
			}
		}
		fCosto = 0;
		#pragma omp parallel for num_threads(hilos) reduction(+:fCosto)
		for (int i = 0; i < m; ++i) {
			fCosto += calcDisMin(vectorDatos[i], vectorCentroides);
		}
		contador = 0;
	}
	vector<vector<double>> centroidesFinales;
	reclusterKmeans2(vectorCentroides,centroidesFinales, hilos, k);
	vectorCentroides.assign(centroidesFinales.begin(), centroidesFinales.end());
	valoresCentroides.resize(vectorCentroides.size());
}
int promediarGrupo(vector<vector<double>>& vectorGrupo, vector<double>& vectorPromedio) {
	int tamGrupo = vectorGrupo.size();
	int tamDato = vectorGrupo[0].size();
	vector<double> anterior(vectorPromedio);
	int cambio = 1;
	vectorPromedio.assign((int)vectorPromedio.size(), 0);
	for (int i = 0; i < tamGrupo; ++i) 
	{
		for (int j = 0; j < tamDato; ++j) {
			vectorPromedio[j] += vectorGrupo[i][j];
		}
	}
	for (int i = 0; i < tamDato; ++i) {
		vectorPromedio[i] /= tamGrupo;
	}
	if (compararVectores(anterior, vectorPromedio)) {
		cambio = 0;
	}
	return cambio;
}

double algoritmoLloyd(vector<vector<double>>& vectorDatos, vector<vector<double>>& vectorCentroides, vector<vector<vector<double>>>& puntosAsociadosC, int hilos, int m, int k, double eps) {
	int posMin = 0;
	double fCosto = 0.0;
	int centroidesAlterados;
	#pragma omp parallel for num_threads(hilos) shared(puntosAsociadosC) private(posMin) reduction(+:fCosto)
	for (int i = 0; i < m; ++i) {
		fCosto += calcDisMin(vectorDatos[i], vectorCentroides);
		posMin = calcMinPos(vectorDatos[i], vectorCentroides);
		#pragma omp critical
		{
			puntosAsociadosC[posMin].push_back(vectorDatos[i]);
		}
	}
	double fCostoPrime = fCosto;
	do {
		fCosto = fCostoPrime;
		centroidesAlterados = 0;
		#pragma omp parallel for num_threads(hilos) reduction(+:centroidesAlterados)
		for (int i = 0; i < (int)puntosAsociadosC.size(); ++i) {
			centroidesAlterados += promediarGrupo(puntosAsociadosC[i], vectorCentroides[i]);
		}
		puntosAsociadosC.clear();
		puntosAsociadosC.resize(k);
		fCostoPrime = 0;
		#pragma omp parallel for num_threads(hilos) shared(puntosAsociadosC) private(posMin) reduction(+:fCostoPrime)
		for (int i = 0; i < m; ++i) {
			fCostoPrime += calcDisMin(vectorDatos[i], vectorCentroides);
			posMin = calcMinPos(vectorDatos[i], vectorCentroides);
			#pragma omp critical
			{
				puntosAsociadosC[posMin].push_back(vectorDatos[i]);
			}
		}
	} while ( ((fCosto - fCostoPrime) > eps) && (centroidesAlterados != 0) );
	return fCostoPrime;
}

void escribirResultados(vector<vector<vector<double>>>& vGrupos, vector<vector<double>>& centroides, double fCosto, double tPared) {
	int cantGrupos = vGrupos.size();
	int n = centroides[0].size();
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
		for (int j = 0; j < (int)vGrupos[i].size(); ++j) 
		{
			for (int k = 0; k < n; ++k) {
				salida << vGrupos[i][j][k];
				if (k != n - 1) {
					salida << ", ";
				}
			}
			salida << endl;
		}
		salida << "Elementos: " << (int)vGrupos[i].size() << endl << endl;
	}
	salida << "Tiempo pared: "<< tPared <<"   Funcion de Costo: " << fCosto;
	salida.close();
}

/*
arg[1] = archivo de entrada
arg[2] = n (dimension datos)
arg[3] = m (cantidad vectores)
arg[4] = k (cantidad centroides)
arg[5] = hilos
arg[6] = Epsilon de parada
*/
int main(int argc, char** argv) {
	if (argc == 7) {
		srand((unsigned int)time(0));
		int n = atoi(argv[2]);
		int m = atoi(argv[3]);
		int k = atoi(argv[4]);
		int hilos = atoi(argv[5]);
		double eps = atof(argv[6]);		
		hilos *= omp_get_num_procs();
		double fCosto = 0.0;
		vector<vector<double>> datos(m, vector<double>(n));
		vector<vector<double>> centroide;
		vector<vector<vector<double>>> valoresAsociadosC (k);
		cargarVectores(argv[1], n, m, datos);
		auto ini = omp_get_wtime();
		kmeansParallelInit(datos, centroide, valoresAsociadosC, hilos, m, k);
		fCosto = algoritmoLloyd(datos, centroide, valoresAsociadosC, hilos, m, k, eps);
		auto fin = omp_get_wtime();
		escribirResultados(valoresAsociadosC,centroide, fCosto, (double)(fin-ini));
	}
	return 0;
}