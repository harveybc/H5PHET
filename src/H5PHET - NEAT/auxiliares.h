/** Funciones gen�ticas auxiliares (busqueda, diagn�stico, etc...) - H file
	funciones que no modifican el genoma ni la poblaci�n (estad�sticas, busquedas).
*/
#ifndef AUXILIARES_H_INCLUDED
#define AUXILIARES_H_INCLUDED
#ifndef PARAMS_H_INCLUDED
#include "params.h"
#endif

float correlac(float* Vc,TConfig* conf);
/** Calcula el coeficiente de correlaci�n de Pearson entre el vector Vc y el vector de valores de entrenamiento
    TODO: FALTA hacer correlac para que retorne promedio o multip de varios vectores de salida(pextender para varias salidas).
*/

void genSSNhdr1(int indexpob, TConfig * conf);
/** genera el header SNN para el genoma indexpob
*/

float randL(TConfig* conf);
/** retorna el siguiente valor de la lista de numeros aleatorios y apunta el contador al siguiente.
*/

int genOrdenEvalF1g(int indexpob, TConfig* conf);
/** genera la lista de conexiones en orden de evaluaci�n

*/
int ordenarListaConexF(hdrSNNv1 headerSNN, tConexDataF* listaConexData, int* ordenEval, int tamOrdenEval);
/** usando el arreglo ordenEval, organiza la lista de conexiones para que queden primero las
// que est�n de �ltimas en el arreglo ordenEval
*/
int genOrdenEvalF(int numGenomas, hdrSNNv1* headerSNN, tConexDataF** listaConexData, int** ordenEval, int* tamOrdenEval, TConfig * conf);
/** genera ordenEval[], orden de evaluaci�n de los nodos de cada genoma (inverso a recursivo empezando de salida)
*/
int recurOrganicer(int indexNodoOut, hdrSNNv1 headerSNN ,tConexDataF* listaConexData, int* ordenEval, char* valCalculado, int* cont, TConfig* conf);
/** genera la lista de nodos a evaluar al hacer un recorrido recursivo y lo guarda enordenEval
*/

void imprimirSeleccion(TConfig* conf);
/** imprime las listas usadas en seleccionCrossover para debigging
*/
float fError(float valor, float salida);
/** retorna el error entre el valor obtenido en una salida y la salida deseada (0,1) o (0,-1)
*/
unsigned contarGPEsp(unsigned indEspecie, TConfig* conf);
/**	retorna el n�mero de genomas en la poblaci�n que pertenecen a la especie indEspecie
*/
unsigned verificarMejor(TConfig* conf);
/** 	Verificar mejor
//retorna 0 si error , se debe usar despu�s de cada evaluaci�n de poblaci�n
*/
void swap(unsigned  int* a, unsigned  int* b);
/** 	unsigned  intercambia los valores apuntados por dos punteros.
*/
unsigned cargarGenoma(unsigned indexpob, char *filename, TConfig* conf);
/** 	lee un genoma desde un archivo y sobreescribe con �l un genoma de conf->pob , la memoria necesaria para el arreglo de nodos y conexiones se gestiona desde esta Funcion.
//Prerequisito:  debe haberse reservado memoria para el genoma en conf->pob[indexpob] (se hace con primeraGen())
//El formato de entrada es(sin separadores): Genoma, Genoma.nodo, genoma.conex las longitudes a escribir
//de cada estructura se basan en el tama�o de Genoma, GenNodoF, GenconexF y en los valores Genoma.totalNodos
//y Genoma.totalConexiones
//Par�metros:	indexpob = indice del arreglo de genomas conf->pob que se va a reemplazar por el le�do
//				filename = path y nombre de archivo del que se leer� el genoma
//Retorna 0 si hay error, 1 ok
*/
unsigned buscarMejorFitness(TConfig* conf);
/** 	busca el mejor fitness entre todos los elementos de la conf->poblaci�n y actualiza lista de representantes
//Retorna el indexpob del genoma con el mejor fitness entre toda la conf->poblaci�n
//Par�metros:	sPob = tama�o de la conf->poblaci�n.
*/
long unsigned calcularMemoriaUsada(unsigned sPob, TConfig* conf);
/** 	calcula la memoria usada por los arreglos del programa en bytes
//Par�metros:	sPob
//Retorna el n�mero de bytes usados por el programa en memoria.
*/
void calcularD(TConfig* conf);
/**		calcula los coeficientes usados en el tSigma escogido
*/

float fSigma(float fX, unsigned param, float fD, TConfig* conf);
/** 	retorna un float corespondiente a la funcion de activaci�n seleccionada con param para una entrada X
	//Param =	0 = sigma(0,1),		y = 1 / (1 + exp (- D * x))
	//		1 = sigma aprox (0,1)	y = 0.5 + x * (1 - abs(x) / 2), y=0 si x<=-1, y=1 si x>=1
	//		2 = tanh(-1,1),		y = 2 / (1 + exp(-2 * x)) - 1
	//		3 = gauss(-1,1),	y = exp(- x * x)
	//		4 = elliot,(-1,1)	y = x / (1 + |x|)
	//		5 = elliot (0,1 )	y = (x / 2) / (1 + |x|) + 0.5
*/
void imprimirGenoma(unsigned index, TConfig* conf);
/** 	imprime la principal informaci�n del genoma incluyendo nodos y conexiones
*/
unsigned verificarGenoma(unsigned index, TConfig* conf);
/** 	verifica que los m�ximos innov num y m�ximos numNodo y conex correspondan con los que se encuentran en el genoms
//retorna 0 si hay error, 1 si OK
*/
unsigned verificarPob(TConfig* conf);
/** 	verifica los innovnum de todos los genomas de la pob
//retorna 0 si hay error, 1 si OK
*/
unsigned buscarInnovNodo(unsigned indexpob, unsigned innovNum, TConfig* conf);
/** 	retorna el index del arreglo de nodos en un genoma en la conf->poblaci�n, retorna -1 si no lo encuentra
.*/
unsigned buscarInnovConex(unsigned indexpob, unsigned innovNum, TConfig* conf);
/** 	retorna el index del arreglo de conexiones en un genoma en la conf->poblaci�n, retorna -1 si no lo encuentra.
//Par�metros: 	indexpob 	= �ndice del genoma en el arreglo de conf->poblaci�n.
//				innovNum 	= n�mero de innovaci�n buscado.
*/
unsigned buscarInnovConexPorNodos(unsigned indexpob, unsigned innovIn, unsigned innovOut, TConfig* conf);//no usa funciones
/** 	retorna el index del arreglo de conexiones en un genoma en la conf->poblaci�n, retorna -1 si no lo encuentra.
//Par�metros: 	indexpob 	= �ndice del genoma en el arreglo de conf->poblaci�n.
//				innovIn 	= n�mero de innovaci�n de nodo de Entrada Buscado.
//				innovOut 	= n�mero de innovaci�n de nodo de Salida Buscado.
*/
int generarDesdeMasterTrainer(TConfig* conf );
/**     genera el archivos de entrenamiento de neat a partir del n�mero de archivos de entradas deseado
// si se necesita se pueden adem�s usar como entradas de neat las entradas diferenciales (2xnumEntradas).
// retorna 0 si error, 1 si ok
*/
#endif
