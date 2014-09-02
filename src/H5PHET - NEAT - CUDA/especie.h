/** Funciones para menejo de especies - H file
	usan/modifican especies únicamente (pueden usars operaciones de genomas o genes).
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
//					n = Número de genes en el genoma más grande.
//					E = Número de genes de exceso
//					D = Número de genes disjounsigned  int
//					W = Promedio de diferencia de pesos entre genes correspondientes. Puede ser seteado a 1 pasra valores que no sean excesivamente grandes.
//Parámetros: 	indexpob1 , indexpob2, = index de los genomas en la conf->pob.
//				c1 = constante de proporcionalidad para Excess genes
//				c2 = constante de proporcionalidad para Disjounsigned genes
//				c3 = constante de proporcionalidad para el promedio de diferencias de pesos.
//				eG_t = threshold para considerar número de genes "excesivamente" grandes (y hacer n=1)
*/
unsigned especieMinDist(unsigned indexpob,float c1, float c2, float c3,unsigned eG_t, TConfig* conf);
/** 	Retorna la especie del genoma al que se tiene la mínima distancia entre los conf->representantes
//Parámetros: indexpob = indice del genoma para el cual se busca la mindist.
//				c1 = constante de proporcionalidad para Excess genes
//				c2 = constante de proporcionalidad para Disjounsigned genes
//				c3 = constante de proporcionalidad para el promedio de diferencias de pesos.
//				eG_t = threshold para considerar número de genes "excesivamente" grandes (y hacer n=1)
*/
float distEspecieCercana(unsigned indexpob,float c1, float c2, float c3,unsigned eG_t, TConfig* conf);  //OPTIMIZADA
/**     retorna la distancia a la especie más cercana diferente a la del indexpob
//Parámetros: indexpob = indice del genoma para el cual se busca la mindist.
//				c1 = constante de proporcionalidad para Excess genes
//				c2 = constante de proporcionalidad para Disjounsigned genes
//				c3 = constante de proporcionalidad para el promedio de diferencias de pesos.
//				eG_t = threshold para considerar número de genes "excesivamente" grandes (y hacer n=1)
*/
unsigned asignarEspecie(unsigned indexpob, float espThreshold,float c1, float c2, float c3,unsigned eG_t, TConfig* conf);
/** 	Usando especieMinDist obtiene la especie más ccompatible, y la compara su distancia con el threeshold, si es menor asigna la especie, si es mayor, crea
//una nueva y la asigna al genoma en cuestión
//Retorna el número de la especie asignada.
//Parámetros:	indexpob	= indice del genoma al que se asignará la especie
//				threeshold	= límite(inferior) de distancia para pertenecier a una especie (//TODO: Variación dinámica de este parám)
*/
unsigned evaluarEspecie(unsigned inicializar,unsigned especie,unsigned maxBufferSize,char *fileNameGTDv1, TConfig* conf); // OPTIMIZADA, NMR
/** Evalúa todos los genomas de una especie en los archivos de entrada y salida y calcula sus valores post Sigma y sus fitness.
*/
unsigned nuevaInnovNodo(unsigned nodoIn, unsigned nodoOut, TConfig* conf);
/** 	retorna el número innovación(número de nodo) buscandolo en la lista de innovaciones de Nodo, si no lo encuentra, lo crea.
//Parámetros: nodoIn, nodoOut = nodos de conexión unsigned  interrumpida para agregar el nuevo nodo.
//RETORNA -1 si hay error, numero de innovación si ok.
//OPTIMIZACIÓN: Puede univicarse nuevaInnovNodo y nuevaInnovCon en una sola pasando el puntero a la lista correspondiente como parámetro.
*/
unsigned nuevaInnovCon(unsigned nodoIn, unsigned nodoOut, TConfig* conf);
/** 	retorna el número de innovación para una conexión buscandolo en la lista de innovaciones de conex. si no lo encuentra lo crea.
//retrona -1 si hay error
//Parámetros: nodoIn,nodoOut = nodos (innovNums ) de origen y destino de la conex.
//OPTIMIZACIÓN: Puede univicarse nuevaInnovNodo y nuevaInnovCon en una sola pasando el puntero a la lista correspondiente como parámetro.
*/

#endif // ESPECIE_H_INCLUDED
