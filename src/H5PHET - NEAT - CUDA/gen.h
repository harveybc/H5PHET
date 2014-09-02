/** Funciones de manejo de genes - H file
	usan/modifican genes �nicamente.
*/
#ifndef GEN_H_INCLUDED
#define GEN_H_INCLUDED
#ifndef PARAMS_H_INCLUDED
#include "params.h"
#endif


void calcularLimThOneNode(unsigned indexPob, unsigned nodeNum, TConfig* conf);
/**		calcula los l�mites maxTh y minTh para un nodo de un genoma , tambi�n altera Th si se sale de los l�mites..
*/
void calcularLimThOneGenome(unsigned indexPob, TConfig* conf);
/**		calcula los l�mites maxTh y minTh para un enoma de la pob, tambi�n altera Th si se sale de los l�mites..
*/
void calcularLimTh(TConfig* conf);
/**		calcula los l�mites maxTh y minTh para cada nodo de cada genoma de la pob, tambi�n altera Th si se sale de los l�mites..
// se debe llamar DESPUES DE evaluarpob, es decir, una sola vez al final del ciclo ppal. y una vez despu�s de primeragen.
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
/**  Perturba uniformemente en un valor aleatorio m�ximo de un porcentaje del rango de pesos, un porcentaje del total de conexiones conexiones con
// probabilidad uniforme, lineal o cuadr�tica dependiente la antiguedad de la conexi�n. Tambi�n se cuenta con la opci�n de decrementar el porcentaje
// de conexiones a perturbqar si se supera cierto n�mero de conexiones en total para evitar grandes cambios de funcionamiento en genomas grandes.
// Par�metros: 	index = indexpob del genoma al que se le aplicar� la perturbaci�n.
//				porcent = (0,1) porcentaje del rango inicial de pesos m�ximo con el que se perturbar�n (dependiendo de tipoProb) las conexiones seleccionadas
//				porcentConex = (0,1) porcentaje m�ximo de conexiones que se perturbar�n uniformemente, el resto se coloca aleatorio en rango inicial  (def=0.9)
//				usarQuadPorcent = 	si es 1 usa el param quadBase y porcentConex para calcular el porcentaje de conexiones a perturbar, puede
//									servir paqra genomas grandes, donde se deben modificar menos conexiones a medida que el genoma incrementa el n�mero de conexiones.
//				quadBase = n�mero de conexiones del genoma a partir del cual se usa la regla de reducci�n cuadr�tica del porcentaje de conex perturbadas
//				tipoProb = 	0 = uniforme : selecciona aleatoriamente en distribuci�n de probabilidad uniforme las conexiones a modificar (def=0)
//							1 = lineal : selecciona aleatoriamente en una distribuci�n lineal que incrementa desde minProb en las conex m�s antiguas hasta 100% en las conex m�s nuevas.
//							2 = cuadr�tica : selecciona aleatoriamente en una distrib cuadratica que se incrementa desde minProb para la conex m�s antigua hasta 100% en la m�s nueva.
//							3 o m�s = exponencial : selecciona aleatoriamente en una distribuaci�n exponencial que se incrementa desde minProb hasta 100% para la conex m�s nueva.
//				minProb = m�nimo de probabilidad para la conexi�n m�s antigua que se e3escoge para las distribuciones lineal, cuadr�tica o exponencial.

*/
unsigned randomizarPesos(unsigned indexpob, float minPeso, float maxPeso, TConfig* conf);
/** 	randomiza los pesos de cada conexi�n del genoma especificado por indexpob
//Retorna 0 si hay error, 1 si ok.
//Par�metros: 	float minPeso = l�mite inferior del rango de n�meros aleatorios
//				float maxPeso = l�mite superior del rango de n�meros aleatorios
*/
#endif //GEN_H_INCLUDED
