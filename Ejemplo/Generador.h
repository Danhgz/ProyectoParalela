#pragma once
#include <vector>
#include <omp.h>
#include <iostream>
using namespace std;

class Generador
{
public:
	int hilos;

	double genereDistanciaCent(int, int, vector < vector <double> >, int);
	//void genereCosto(vector < vector <double> >, vector < vector <double> >, double&, int, int,double);
	double recalculeCosto(int,int, vector < vector <double> >, vector < vector <double> >, int, int);
	void imprime_dato(vector<int>, vector < vector <double> >, vector < vector <double> >, int);
	void imprime_grupos(vector<vector <int>>, vector < vector <double> >, vector < vector <double> >,int);
	Generador();
	~Generador();
};

