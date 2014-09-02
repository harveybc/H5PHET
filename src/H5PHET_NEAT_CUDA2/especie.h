/** Funciones para menejo de especies - H file
	usan/modifican especies �nicamente (pueden usars operaciones de genomas o genes).
*/
#ifndef ESPECIE_H_INCLUDED
#define ESPECIE_H_INCLUDED
#ifndef PARAMS_H_INCLUDED
#include "params.h"
#endif
void imprimirListasInnov(TConfig* conf);
/**     Imprime la listaInnovCon y listaInnovNodo con todos sus miembros y dimesiones.
*/
unsigned verificarListasInnov(TConfig* conf);
/** Retorna 0 si hay error en listaInnovCon o listaInnovNodo
*/
float calcularDist(unsigned indexpob1,unsigned indexpob2,float c1, float c2, float c3,unsigned eG_t, TConfig* conf);
/** 	Retorna la distancia entre dos genomas. uusando la formula D=((((C1)*E)/n)+(((C2)*D)/n)+C3*W)
//			Donde: 	C1,C2,C3 = constantes de proporcionalidad (mirar ejemplos para valores iniciales.)
//					n = N�mero de genes en el genoma m�s grande.
//					E = N�mero de genes de exceso
//					D = N�mero de genes disjounsigned  int
//					W = Promedio de diferencia de pesos entre genes correspondientes. Puede ser seteado a 1 pasra valores que no sean excesivamente grandes.
//Par�metros: 	indexpob1 , indexpob2, = index de los genomas en la conf->pob.
//				c1 = constante de proporcionalidad para Excess genes
//				c2 = constante de proporcionalidad para Disjounsigned genes
//				c3 = constante de proporcionalidad para el promedio de diferencias de pesos.
//				eG_t = threshold para considerar n�mero de genes "excesivamente" grandes (y hacer n=1)
*/
unsigned especieMinDist(unsigned indexpob,float c1, float c2, float c3,unsigned eG_t, TConfig* conf);
/** 	Retorna la especie del genoma al que se tiene la m�nima distancia entre los conf->representantes
//Par�metros: indexpob = indice del genoma para el cual se busca la mindist.
//				c1 = constante de proporcionalidad para Excess genes
//				c2 = constante de proporcionalidad para Disjounsigned genes
//				c3 = constante de proporcionalidad para el promedio de diferencias de pesos.
//				eG_t = threshold para considerar n�mero de genes "excesivamente" grandes (y hacer n=1)
*/
float distEspecieCercana(unsigned indexpob,float c1, float c2, float c3,unsigned eG_t, TConfig* conf);  //OPTIMIZADA
/**     retorna la distancia a la especie m�s cercana diferente a la del indexpob
//Par�metros: indexpob = indice del genoma para el cual se busca la mindist.
//				c1 = constante de proporcionalidad para Excess genes
//				c2 = constante de proporcionalidad para Disjounsigned genes
//				c3 = constante de proporcionalidad para el promedio de diferencias de pesos.
//				eG_t = threshold para considerar n�mero de genes "excesivamente" grandes (y hacer n=1)
*/
unsigned asignarEspecie(unsigned indexpob, float espThreshold,float c1, float c2, float c3,unsigned eG_t, TConfig* conf);
/** 	Usando especieMinDist obtiene la especie m�s ccompatible, y la compara su distancia con el threeshold, si es menor asigna la especie, si es mayor, crea
//una nueva y la asigna al genoma en cuesti�n
//Retorna el n�mero de la especie asignada.
//Par�metros:	indexpob	= indice del genoma al que se asignar� la especie
//				threeshold	= l�mite(inferior) de distancia para pertenecier a una especie (//TODO: Variaci�n din�mica de este par�m)
*/
unsigned evaluarEspecie(unsigned inicializar,unsigned especie,unsigned maxBufferSize,char *fileNameGTDv1, TConfig* conf); // OPTIMIZADA, NMR
/** Eval�a todos los genomas de una especie en los archivos de entrada y salida y calcula sus valores post Sigma y sus fitness.
*/
unsigned nuevaInnovNodo(unsigned nodoIn, unsigned nodoOut, TConfig* conf);
/** 	retorna el n�mero innovaci�n(n�mero de nodo) buscandolo en la lista de innovaciones de Nodo, si no lo encuentra, lo crea.
//Par�metros: nodoIn, nodoOut = nodos de conexi�n unsigned  interrumpida para agregar el nuevo nodo.
//RETORNA -1 si hay error, numero de innovaci�n si ok.
//OPTIMIZACI�N: Puede univicarse nuevaInnovNodo y nuevaInnovCon en una sola pasando el puntero a la lista correspondiente como par�metro.
*/
unsigned nuevaInnovCon(unsigned nodoIn, unsigned nodoOut, TConfig* conf);
/** 	retorna el n�mero de innovaci�n para una conexi�n buscandolo en la lista de innovaciones de conex. si no lo encuentra lo crea.
//retrona -1 si hay error
//Par�metros: nodoIn,nodoOut = nodos (innovNums ) de origen y destino de la conex.
//OPTIMIZACI�N: Puede univicarse nuevaInnovNodo y nuevaInnovCon en una sola pasando el puntero a la lista correspondiente como par�metro.
*/

#endif // ESPECIE_H_INCLUDED
