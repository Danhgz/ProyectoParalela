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
		int i= 0,j = 0;
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

double calcularDistanciaEuclidiana(vector<double> vec1, vector<double> vec2) {
	double distancia = 0;
	for (int i = 0; i < vec1.size(); ++i) { 
		distancia += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
	}
	return distancia;
}

void seleccionInicialCentroides(vector<vector<double>>& vectorDatos, vector<vector<double>>& vectorCentroide, int n, int m, int k) {
	default_random_engine generator;
	uniform_int_distribution<int> distribution(0, m-1);
	vectorCentroide.push_back(vectorDatos[distribution(generator)]);
	vector<double> vectorCentroidePrime;
	double l = k / 2;
	double psi = 0;
	for (int i = 0; i < m; ++i) {
		psi += calcularDistanciaEuclidiana(vectorDatos[i],vectorCentroide[0]);
	}
	for (int i = 0; i < 5; i++) { //O(log(psi)?????
		
	}
}

void algoritmoLloyd() {

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
*/
int main(int argc, char** argv){
	if (argc == 5) {
		int n = atoi(argv[2]);
		int m = atoi(argv[3]);
		int k = atoi(argv[4]);
		vector<vector<double>> datos(m, vector<double>(n));
		vector<vector<double>> centroide;
		cargarVectores(argv[1], n, m, datos);
		imprimirVectores(n, m, datos);
		seleccionInicialCentroides(datos,centroide, n, m, k);
		algoritmoLloyd();
	}
	return 0;
}