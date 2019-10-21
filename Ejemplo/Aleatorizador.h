//
//  Aleatorizador.h
//  OMP_PryArribada
//
//  Created by administrador on 2/5/19.
//  Copyright © 2019 administrador. All rights reserved.
//

#ifndef Aleatorizador_h
#define Aleatorizador_h

#include <random>
#include <chrono>
#include <assert.h>     /* assert */
using namespace std;

class Aleatorizador {
public:
    typedef std::chrono::high_resolution_clock t_clock;
    static default_random_engine generador; // generador de números aleatorios
    static uniform_real_distribution<double> random_uniform_real; // generador basado en distribucion uniforme
    static double random_logistic(double location, double scale); // generador basado en distribucion logistic
    static void inicializar_generador_random(); // inicializador del generador
};

default_random_engine Aleatorizador::generador;
uniform_real_distribution<double> Aleatorizador::random_uniform_real(0.0, 1.0);

double Aleatorizador::random_logistic(double location, double scale)
{
    assert(scale > 0.);
    assert(location >= 0.);
    return location - scale * log(1. / random_uniform_real(generador) - 1.);
}

void Aleatorizador::inicializar_generador_random()
{
    t_clock::time_point beginning = t_clock::now();
    t_clock::duration lapso = t_clock::now() - beginning;
    unsigned semilla = lapso.count();
    generador.seed(semilla);
}
#endif /* Aleatorizador_h */
