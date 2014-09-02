/** Funciones de manejo de genes - H file
	usan/modifican genes únicamente.
*/
#ifndef GEN_H_INCLUDED
#define GEN_H_INCLUDED
#ifndef PARAMS_H_INCLUDED
#include "params.h"
#endif


void calcularLimThOneNode(unsigned indexPob, unsigned nodeNum, TConfig* conf);
/**		calcula los límites maxTh y minTh para un nodo de un genoma , también altera Th si se sale de los límites..
*/
void calcularLimThOneGenome(unsigned indexPob, TConfig* conf);
/**		calcula los límites maxTh y minTh para un enoma de la pob, también altera Th si se sale de los límites..
*/
void calcularLimTh(TConfig* conf);
/**		calcula los límites maxTh y minTh para cada nodo de cada genoma de la pob, también altera Th si se sale de los límites..
// se debe llamar DESPUES DE evaluarpob, es decir, una sola vez al final del ciclo ppal. y una vez después de primeragen.
	//Para cada indexpob, pera cada conex, acumule max y minTh para el nodo de destino de la conex.
	//si tsigma es -1,1 se debe colocar el max=abs(sum de negativos) + sumpositivos y el min=-max
	//tsigma=0 = sigma(0,1),		y = 1 / (1 + exp (- D * x))
	//		1 = sigma aprox (0,1)	y = 0.5 + x * (1 - abs(x) / 2), y=0 si x<=-1, y=1 si x>=1
	//		2 = elliot (0,1 )	y = (x / 2) / (1 + |x|) + 0.5
	//		3 = gauss(-1,1),	y = exp(- x * x)
	//		4 = elliot,(-1,1)	y = x / (1 + |x|)
	//		5 = tanh(-1,1),		y = 2 / (1 + exp(-2 * x)) - 1
*/
void perturbarPeso(unsigned index, float porcent, float porcentConex, unsigned usarQuadPorcent, unsigned quadBase, unsigned tipoProb, float minProb, TConfig* conf); //OPTIMIZADA: //TODO otros tipos de pert, quadPorcent y hacer Funcion perturbarTh si es necesaria
/**  Perturba uniformemente en un valor aleatorio máximo de un porcentaje del rango de pesos, un porcentaje del total de conexiones conexiones con
// probabilidad uniforme, lineal o cuadrática dependiente la antiguedad de la conexión. También se cuenta con la opción de decrementar el porcentaje
// de conexiones a perturbqar si se supera cierto número de conexiones en total para evitar grandes cambios de funcionamiento en genomas grandes.
// Parámetros: 	index = indexpob del genoma al que se le aplicará la perturbación.
//				porcent = (0,1) porcentaje del rango inicial de pesos máximo con el que se perturbarán (dependiendo de tipoProb) las conexiones seleccionadas
//				porcentConex = (0,1) porcentaje máximo de conexiones que se perturbarán uniformemente, el resto se coloca aleatorio en rango inicial  (def=0.9)
//				usarQuadPorcent = 	si es 1 usa el param quadBase y porcentConex para calcular el porcentaje de conexiones a perturbar, puede
//									servir paqra genomas grandes, donde se deben modificar menos conexiones a medida que el genoma incrementa el número de conexiones.
//				quadBase = número de conexiones del genoma a partir del cual se usa la regla de reducción cuadrática del porcentaje de conex perturbadas
//				tipoProb = 	0 = uniforme : selecciona aleatoriamente en distribución de probabilidad uniforme las conexiones a modificar (def=0)
//							1 = lineal : selecciona aleatoriamente en una distribución lineal que incrementa desde minProb en las conex más antiguas hasta 100% en las conex más nuevas.
//							2 = cuadrática : selecciona aleatoriamente en una distrib cuadratica que se incrementa desde minProb para la conex más antigua hasta 100% en la más nueva.
//							3 o más = exponencial : selecciona aleatoriamente en una distribuación exponencial que se incrementa desde minProb hasta 100% para la conex más nueva.
//				minProb = mínimo de probabilidad para la conexión más antigua que se e3escoge para las distribuciones lineal, cuadrática o exponencial.

*/
unsigned randomizarPesos(unsigned indexpob, float minPeso, float maxPeso, TConfig* conf);
/** 	randomiza los pesos de cada conexión del genoma especificado por indexpob
//Retorna 0 si hay error, 1 si ok.
//Parámetros: 	float minPeso = límite inferior del rango de números aleatorios
//				float maxPeso = límite superior del rango de números aleatorios
*/
#endif //GEN_H_INCLUDED
