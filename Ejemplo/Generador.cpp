#include "Generador.h"

double Generador::genereDistanciaCent(int pos, int c, vector < vector <double> > d, int n) {
	double valor = 0;
	for (int j = 0; j < n; j++)
	{
		valor += sqrt( pow((d[pos][j] - d[c][j]), 2) );			// Sumatoria de (x-x) al cuadrado (distancia)
	}
	return valor;
}

//void Generador::genereCosto(vector < vector <double> > centroides, vector < vector <double> > data, double &cost, int m_, int n_,double c_) {
//		double valor = 0;
//		double minimo = 100000;
//#pragma omp barrier
//#pragma omp for reduction(+:c_) private(valor, minimo)
//			for (int i = 0; i < m_; i++)
//			{
//				for (int cent = 0; cent < centroides.size(); cent++) {
//					minimo = 100000;
//					for (int j = 0; j < n_; j++)
//					{
//						valor = pow((data[i][j] - centroides[cent][j]), 2);			// Sumatoria de (x-x) al cuadrado (distancia)
//						valor = sqrt(valor);
//						if (valor < minimo) {
//							minimo = valor;
//						}
//					}
//#pragma omp atomic
//					c_ += minimo;
//				}
//			}
//#pragma omp single
//		cost = c_;
//}
//
double Generador::recalculeCosto(int pos, int pos_centroide,vector < vector <double> > vec, vector < vector <double> > centroides, int m_, int n_) {
	double c = 0;
	double valor = 0;
	double minimo = 100000;
	for (int j = 0; j < n_; j++)
	{
		valor = pow((vec[pos][j] - centroides[pos_centroide][j]), 2);			// Sumatoria de (x-x) al cuadrado (distancia)
		valor = sqrt(valor);
		if (valor < minimo) {
			minimo = valor;
		}
	}
	c = minimo;                  
	return c;
}

void Generador::imprime_dato(vector<int> tmp, vector < vector <double> > data, vector < vector <double> > centroides, int n_) {
	int cont = 0;
	for (auto it = tmp.begin(); it != tmp.end(); it++) {
		for (int j = 0; j < n_; j++) {
			if ( it == tmp.begin() ) {
				cout << centroides[*it][j] << ", ";
			}
			else {
				cout << data[*it][j] << ", ";
				cont++;
			}
		}
		cout << endl;
	}
	cont = cont / n_;
	cout << cont << endl;
}

void Generador::imprime_grupos(vector<vector <int>> g_centroid, vector < vector <double> > data, vector < vector <double> > centroides,int n_) {
	for (auto it = g_centroid.begin(); it != g_centroid.end(); it++) {
		imprime_dato(*it, data, centroides, n_);
		cout << endl;
	}
}

Generador::Generador()
{
}


Generador::~Generador()
{
}
