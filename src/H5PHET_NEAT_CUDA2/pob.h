/** Archivo de manejo de población, incluye ciclo principal - H file
	usan/modifican la poblacion únicamente (pueden usarse operaciones de especies, genomas o genes).
*/
#ifndef POB_H_INCLUDED
#define POB_H_INCLUDED
#ifndef PARAMS_H_INCLUDED
#include "params.h"
#endif

//void evalPobCUDA(int numEntradas, int numSalidas, int numBias, int numDatos, float* dataGTDf, tConexDataF** listaConexData, int* tamListaConexPost, float** valorC, float** valorTr, float* fitness);


int evalPob(int numGenomas, int numDatos, hdrGTDv1 headerGTD, hdrSNNv1* headerSNN, float** dataGTD, tConexDataF** listaConexData, float* fitness, int** ordenEval, int* tamListaConexPost, TConfig* conf);
/** evalua por mutiplicación de matrices la población proporcionada en un buffer usando
    otro buffer con los datos de entrenamiento.
    Params:
            numGenomas = número de genomas en el buffer de Conexiones
            numDatos = número de conjuntos de registros de entrenamiento que contiene el buffer **dataGTD
            headerGTD = header del archivo GTD de entrenamiento.
            headerSNN[numGenomas] =  arreglo de headers SNN de los genomas de la Pob
            dataGTD[numDatos][numEntradas+numSalidas] = buffer con datos GTDv1
            listaConexData[i<numGenomas][numConex[i]] = buffer con datos de conexiónes por genoma i.
            fitness[numGenomas] = valor de retorno, arreglo de fitness calculado para cada genoma.
*/

unsigned distribProc(TConfig* conf);
/** actualiza el mejor de cada especie desde NNP de reps para procesamiento distribuido
*/
unsigned cargarDesdeNNP(char* filename, unsigned especieDestino, unsigned cantidad, TConfig* conf);
/** carga una población inicial desde archivo NNP (Neural Network Population) al inicio del arreglo pob
*/
int inicializaciones(TConfig* conf);
/** // inicializa todas las variables de configuración
*/
unsigned cicloPrincipal(TConfig* conf);
/** 	Prerequisisto: Tener una conf->población inicial usando primeraGen()
//Realiza el ciclo principal de NEAT, realiza  evaluación, especiación,seleccion(y cruce) y mutación AN + AC (con probabilidades de entrada)
//hasta que se cumpla maxIteraciones o minFitness se alcance.
//Parámetros:
//			probMutAN = entre 0 y 1 prob de mutación AN
//			probMutAC = entre 0 y 1 prob de mutación AC
//			maxIteraciones = máximo número de veces que debe correrse el ciclo
//			minFitness = entre 0 y 1 mínimo fitness (promediado durante el arch de entrenamiento) necesario para detener el ciclo
//			maxMemoriaUsada = Max memoria que se puede utilizar en MBytes (si se supera, sale del ciclo)
//			fileNameGTDv1 = Path para el archivo de datos de entrenamiento en formato GTDv1.
//			filenameSNNv1= Path  para el mejor genoma encontrado al  omento de terminar el ciclo principal.
//			repTrain = número de repeticiones del archivo de entrenamiento(para cada genoma).
//			maxBufferSize = tamáño máximo del buffer de lectura de archivos de entrada., memoria usada = maxBufferSize*4 bytes
//Retorna -1 si hay error, indexpob de genoma con mayor fitness si se cumple una de las condiciones de parada.
*/
unsigned evaluarPob(unsigned inicializar,unsigned primero,unsigned maxBufferSize, char *fileNameGTDv1, TConfig* conf);
/** 	evalúa toda la población en todo el archivo de entradas y salidas, deja los valores resultantes en el campo valor de cada gen nodo de cada genoma
//y calcula el fitness basado en el que se va acumulando con cada evaluación de cada genoma.
//Los archivos de entrada y salida deben estar previamente abiertos para lectura binaria br
//retorna 0 si hay error, 1 si ok.
////TODO: parámetro repeticiones para pasar los archivos de entrada repetidas veces por las redes neuronales al realizar las evaluaciones.
*/
unsigned seleccionCrossover(TConfig* conf);
/** 	realiza selección y cruce en toda la población
//Inicializa el arreglo conf->numGenomasPorEspecie (conf->fInal de pag 394(426) de AI game prog) y lo modifica para reducir a la mitad el número de especies
//en conservación. Y reparte el resto entre las especies que nó están en conservación.
//Requiere haber evaluado fitness y haber realizado especiacion.
//Retorna 0 si hay error, 1 si ok
//Parámetros:	porcentRedConserv = (entre 0 y 1) porcentaje de genomas que se quitan al número calclulado para especies en conserv. recom=0.3
			//	porcentElim = (entre 0 y 1) porcentaje de genomas que se eliminan en cada generación, el restante porcentaje se reproduce por cruce.
			//	unsigned  intentosPareja= número de unsigned  intentos para buscar pareja aleatoriamente entre los padres antes de "matrimonio forzoso :) "
			//	super = float entre 0 y 1, Probabilidad de heredar los exess y disjounsigned  ints del menos apto (aparte de los que se heredan normalmente del más apto)
			//	promendiarProb = float entre 0 y 1 = probabilidad de que en caso de matching, se promedien los pesos en lugar de
				//seleccionarlos aleatoriamente entre los padres.
*/
float especiacion(TConfig* conf);
/** 	realiza la asignación de especies (sobreescribiendo las existntes) para todos los genomas de la población
//Se debe realizar DESPUES de la evaluación.
//Altera el valor de threshold en un máximo porcentaje especificado para alcanzar el número de especies requerido.
//si el número de especies requerido ya se alcanzó,  no crea nuevas especies sino que asigna a caga genoma la especie más cercana.
//actualiza en cada ejecución la lista de conf->representantes y el número de generaciones sin mejora de fitness por especie.
//también actualiza la lista de especies en conservación.
//Parámetros:	sPob = número de genomas en la conf->población
//				numEspDeseadas = número de especies deseadas
//				numConserv = número de especies en conservación.
//				threshold = threshold inicial
//				porcentVarTh = entre 0 y 1 = porcentaje de variación del threshold si no se ha alcanzado el máximo número de especies.
//				c1= constante de proporcionalidad en distencia para número de genes excess entre los padres
//				c2= constante de proporcionalidad en distancia para número de genes disjounsigned entre los padres
//				c3= constante de proporcionalidad en distancia para promedio de diferencias de pesos en matching genes de los padres
//				eG_t=(creo que no es necesario REVISAR) número de genes necesarios para considerar el genoma suficientemente grande y hacer n=1
//retorna 0 si hay error, 1 si ok
*/
unsigned primeraGen(unsigned tamPob, unsigned nEntradas, unsigned nSalidas, unsigned nBias, float minPeso, float maxPeso, unsigned maxMutacionesAN,unsigned maxMutacionesAC,unsigned maxIntentosMutAC, short unsigned useMutarAC, short unsigned useMutarAN,unsigned useRandomization, TConfig* conf);
/** 	genera una población inicial de genomas de nIn entradas, nOut salidas y nBias bias a  partir de la mutación
// (AN + AC) de un genoma inicial (en indice 0) totalmente conectado y que es el representante de la especie 0
//Además muta aleatoriamente los pesos de llas conexiones de todos los genomas.
//retorna 0 si hubo algún error en el posicionamiento en memoria de las estructuras de los genomas.
//necesita: funcion (//TODO) mutarAN(), mutar(AC), genomaInicial(), nuevoNodo(), nuevaConex(), especieMinDist(),asignarEspecie()
////TODO: se deben hacer dos arreglos para los genomas de los conf->representantes de cada especie de la generación actual y para
//los conf->representantes de las especies de la generación anterior.
//En cada generación se debe comparar cada genoma con los conf->representantes de la generación anterior para determinar la especie a la que pertenecen
//y se toma el genoma con menor error como representante de la especie en el arreglo actual.
//Se debe llevar un record para cada representante de cada especie del número de generaciones que lleva sin mejorar, para poder eliminaar de esta
//manera especies que se quedaron estancadas(EXCEPTO LA MEJOR(o n mejores?)).
//El número de hijos que puede producir cada especie depende del fitness de sus individuos comparado con el promedio de fitness total como se
//muestra en la pag 394 del libro de AI game programming. (great tool)
//Para la primera generación: Averiguar
//Algoritmo para especiar toda la conf->población después de mutación y cruce, en pag 54 de disertación de PhD
*/
unsigned inicializarPob(unsigned tamPob,unsigned nEntradas,unsigned nSalidas,unsigned nBias, TConfig* conf);
/** 	inicializa el vector de genomas (conf->poblacción)
//también obtiene memoria para arreglo de nodos y conexiones de tamaño nEntradas*nSalidas*nBias
//Obtiene memoria para cada nueva estructura Genoma.
//Retorna 0 si hay error
*/
int competencia(TConfig* conf);
/**     realiza competencia, que copmara con una copia guardada en indexPob+numEspecies+sizePob
// usa variable conf->porcentCompetencia como prob de ganancia del backup si es mejor que el actual,
// si es peor, simplemente se reemplaza.
// si es negativa, no hace competencia.
// solo compara si los correspondientes indexpob son de la misma especie( para evitar acaparamiento)
// se debe ejecutar después de especiacion().
// retorna 0 si hay error, 1 si ok
*/
#endif //POB_H_INCLUDED
