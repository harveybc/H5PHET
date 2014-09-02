/** Archivo de parámetros globales - H file

*/
#ifndef PARAMS_H_INCLUDED
#define PARAMS_H_INCLUDED

#include <stdlib.h> //para usar calloc, realloc, free, rand, etc...
#include <stdio.h> //Para usar printf, etc...
#include <time.h> //para usar clock() para hacer funcion delay.
#include <math.h> //Para funciones matemáticas (exp,abs, otras.)
#include <string.h> //para usar memset(), memmove, memcpy()...
#include <limits.h> //Para verificar los límites de unsigned y rand
//#include <omp.h>
#include <cuda.h>

//Estructuras:

typedef struct
{
    float valor;
    void* next;
} tRandList;

// estructura del encabezado de GTDv1 (Generic Training Data).
// formato GDTv1: [encabezado][d1_entrada_0]..[d1_entrada_n][d1_salida_0]..[d1_salida_n].._entrada_0]..[dn_entrada_0]..[dn_entrada_n][dn_salida_0]..[dn_salida_n]
typedef struct
{
 	char fileID[3]; // Debe ser siempre GTD para identificar el tipo de archivo.
 	unsigned char version;
 	unsigned usarSigned;
 	unsigned tamRegistros;
 	unsigned numEntradas;
 	unsigned numSalidas;
} hdrGTDv1;

// estructura del encabezado de SNNv1 (Simple Neural Network).
// formato SNNv1: [encabezado][unsigned conexIn[numConex]][unsigned conexOut[numconex]][char enabled[numconex]][int tamListaConex][listaConexData[tamListaConex]][float/double peso[numConex]]
typedef struct
{
 	char fileID[3]; // Debe ser siempre SNN para identificar el tipo de archivo.
 	char version;
 	unsigned usarSigned;
 	unsigned tamRegistros;
 	unsigned numEntradas;
 	unsigned numSalidas;
 	unsigned numBias;
 	unsigned numHiddens;
 	unsigned numConex;
 	double sigmaFactor;
 	double actThreshold;
 	double lastFitness; // usado para programación distribuida.
} hdrSNNv1;
//Estructura tConexDataF, usada para los datos de las conex de SNN
typedef struct
{
    unsigned conexIn;
    unsigned conexOut;
    char enabled;
    float peso;
} tConexDataF;

typedef struct
{
    unsigned conexIn;
    unsigned conexOut;
    char enabled;
    double peso;
} tConexDataD;


// estructura del encabezado de NNPno, v1 (Neural Network Population).

// es una colección de SNNs.
// formato NNPv1: [hdrNNPv1][SNNv1_0][SNNv1_1]...[SNNv1_n]
// donde SNNv1_x = red neuronal simple en formato SNNv1.
typedef struct
{
 	char fileID[3]; // Debe ser siempre NNP para identificar el tipo de archivo.
 	unsigned char version;
 	unsigned numGenomas;
} hdrNNPv1;

typedef struct   // tipo GenConexF gen de conexión.
{
    unsigned innovNum; //es el número de la conexión(puede quitarse si se usa el mimso indice del arreglo de nodos).
    unsigned nodoIn; //innovnum del Nodo de origen de la conexión.
    unsigned nodoOut; //innovnum del Nodo de destino de la conexión.
    float peso; //Peso de la conexión.
    char recurrente; //0=No recurrente, 1=Recurrente.
    char enabled; //0=Desactivado, 1=Activado.
    unsigned indexIn; //index del nodo de entrada en el arreglo de nodos. Necesarios porque al hacer nuevoNodo se hace realloc y los punteros quedan apuntando a un valor libre.
    unsigned indexOut; //index del nodo de salida en el arreglo de nodos. Se utilizan en evaluarGenoma, se modifican en nuevonodo y nuevaconex
    //TODO: Cambiar todos los short unsigned por short short unsigned  int. para ahorrar memoria.
} GenConexF;

typedef struct   // tipo GenNodoF gen de nodo. Para SMP se necesita además una lista de los nodos padres y un contador de padres.
{
    float thNodo; //Threshold de activación de la neurona.
    unsigned innovNum;//Obtenido al crear nuevo nodo buscando lista de ;
    short unsigned nodeFunction; //0=entrada,1=oculto,2=salida,3=bias.
    float maxTh; // Mayor valor posible de th para esta neurona, es la sumatoria de los pesos positivos que entran
    float minTh; // Menor valor posible de th para esta neurona, es la sumatoria de los pesos negativos	que entran
    float valor; //Valor de salida de cada nodo(para computar la ann sin matriz de pesos).
    //TODO: hacer realloc de arreglos de indexes cuando se agregue nodo para que quede en
    GenConexF** conexHijo; //arreglo de punteros a conexiones a cada hijo solo se usa antes de la eval.
    unsigned contHijos; //contador de hijos. Se usan solo en evaluarGenoma.
    short unsigned estadoC; // flag que dice si el valor No ha sido calculado (0), está calculandose (1) o ya está calculado (2)
    //TODO: Falta arreglo de punteros a nodos hijos, arreglo de puntero a conexión a hijo y contador de nodos hijos y en genoma nueva func recursiva calcularValorNodo()
} GenNodoF;

typedef struct   // tipo Genoma es un genoma que contiene lista de nodos y conexiones,etc...
{
    GenNodoF *nodo; //Arreglo dinámico de nodos.
    GenConexF *conex; //Arreglo dinámico de conexiones.
    unsigned especie; //Especie a la que pertenece el genoma en la conf->población. Se establece por funcion calcularEspecie() después de la segunda generación.
    unsigned totalNodos;//Mantiene el total de nodos en el gen (Maxindex+1)
    unsigned totalConexiones;//Mantiene el total de conexiones (Maxindex+1)
    unsigned maxInnovNumConex; //Mantiene el máimo número de innovación de genes de conexión.
    unsigned maxInnovNumNodo; //Mantiene el máximo número de  innovación para genes de nodo.
    int numHijos; //Es el número de hijos que debe tener cada genoma por cruce. (se inicializa en -1, se modifica con)
    float fitness;//Usado para almacenar el error del genoma respecto a la salida deseada.
} Genoma;

typedef struct   // tipo TNodoOut de nodos de salida con su respectivo innovNum usado en TListaInnov.
{
    unsigned innovNum;//Número de innovación correspondiente a la salida relacionada con el index (entrada).
    unsigned nodoOut; //
} TNodoOut;

typedef struct   // tipo TListaInnov para listas de innovaciones, la lista es un puntero a esta struct. el index es el número del nodo de entrada.
{
    unsigned numOut; // Para cada index(in), es el número disponible de salidas (para limitar el index de este arreglo).
    TNodoOut *nodoOut;// Es un arreglo de estructuras listOut que contienen el innovNum y el nodo de salida de la conexión.
} TListaInnov;
//TODO: pasar a cada Funcion exáctamente los parámetros que necesita en lugar de todo el conf.

typedef struct   // tipo TConfig de parámetros de NEAT
{
    int usarDistrib; // flag de uso de procesamiento distribuido (0=no, 1= todas las iter, 2=saltando una, etc...)
 	char* distriCmd; // comando para ejecutar después de distribProc
 	char* prefix; // prefijo para los filenames
    unsigned tmpIndexPob; // usado para guardar la posición del indexpob usado como genoma temporal (ultimo del arreglo pob).
    int cargarDistrib;
    FILE *logFile; //usado para handler de archivo de logs.
    float probInterSp; // (default=0.001) Probabilidad de que en una iteración se produzca reproduccón inter-especies.
    unsigned cargarNNP; // (default=0) 1=Carga genomas iniciales desde archivo NNP (NeuralNetworkPopulation)
    float porcentCompetencia; // (default=0.8) -1=no usada. especifica la prob de que el mejor genoma gane entre generaciones (para cada genoma)
    unsigned maxIntentosDistInicial; // (default=10) 0= nu usa maxdistancia máximo número de intentos para máxima distancia de genomas iniciales
    unsigned maxIntentosDist; //máximo (default=5)número de intentos para máxima distancia con cada renew de genoma.
    unsigned maxGeneracSinMejoraPG; // (default=8)máximo número de generaciones sin mejora antes de generarmaxdist
    float probMutMaxDist; //probabilidad de mutación por máximas distancias por cada genoma
	unsigned usarSigned; // 1 si se usarán valores con signo para entrenamiento, 0 si se usan valores sin signo
    unsigned numEntradas;//Especifica el múmero de neuronas de entrada
    unsigned numSalidas;//Especifica el número de neuronas de salida
    unsigned numOcultas;//Especifica el número de neuronas ocultas (se debe actualizar cada vez que se haga un nuevo gen)
    unsigned numBias;//Especifica el número de neuronas de entrada de bias (cte 1)
    float Fthreshold;//Especifica el threshold, debe ser 1 para sigma que retornan [0,1] y -1,1
    unsigned tSigma; //Por defecto usa funcion de activación elliot unitario
    float fSigmaD; //Por defecto usa un D=1 para la funcion sigma Normal float(0)
    float A,D,F; // coeficientes de sigma calculados a mano para que f(0.5)=0.5 f(0)=1 y f(-1)=0 en tsigma (0,1) y f(0.5)=0 y f(-1)=-1 para tsigma (-1,1)
    unsigned useFloat;//indica si se deben usar todos los valores de entrada/salida/calculos tipo float si useFloat=0 se usan tipo double(lento).
    unsigned sizePob; //Tamaño de la conf->población.
//	float *c_t ; //compatibility threshold (mínima distancia para nueva especie.), //TODO, corregir valor inicial basado en ejemplos citados al conf->fInal del archivo
    unsigned spEspecies; //Set Pounsigned para el número de especies que se debe tratar de mantener usando c_t
    float maxFitnessConservacion; ////TODO: hacer funcion para pasar de error a fitness teniendo como parámetro Fitness paa error=100, linealidad: 0=>1-e, 1=>1/x, 2=1/x^2, etc.
    unsigned maxEspeciesConservacion; //Es el máximo número de especies en conservación. si este número se supera debido a buén fitness de otras especies,
    //a las especies excedentes(de menor fitness) se les dá el número de especímenes calculado sin restricciones
    //TODO: falta calcular automáticamente los valores iniciales (rango para escoger aleatorio) de pesos de conexiones a nodos de salidas usando preanálisis de los valores
    // de salida, los thresholds de las neuronas de entrada se pueden calcular aleatoriamente en rangos determinados por los valores máximos y mínimos
    // de entrada que toma una muestra de entrenamiento. Todas las conexiones adicionales tendrán peso en rango de -1 a 1 o de 0 a 1 dependiendo
    //de tSigma.
    float pesoMinInicial; //peso mínimo para valores aleatorios asignados en la primera generación
    float pesoMaxInicial; //peso máximo para valores aleatorios asignados en la primera generación.
    Genoma *pob; //Arreglo conf->población, es un arreglo de genomas.
    unsigned *representantes; //Arreglo de indexpob de conf->representantes para cada especie.
    unsigned ultimaInnov; //Guarda el valor de la última innovación (acumulativa entre nodos y conexiones)
    FILE *fIn;//Puntero a archivo de entradas
    FILE *fOut;//Puntero a archivo de entradas
    char *fileNameIn; //nombre del archivo de entradas
    float porcentEnableds;
    float antMejor;
    unsigned minConexMutPeso;
    float maxDesvThEspecies;  //0 hasta 1 max desviación desde threshold de compatibilidad.
    float probMutPeso;
    float probMutTh;
    float porcentMutPeso;
    float porcentMutNTh;

    //TODO: separar en estructura control todas las variables globales usadas por el prog y en conf las de configuración para simplificar config.
    unsigned actualizarEnCambios; //0 solo actualiza una vez en ciclo ppal después de evaluarPob antes de realizar cruce , >0 evalúa en cada cambio de la población.
    unsigned contInnovNodo; // Es el número de elementos de la lista conf->listaInnovNodo (y también es el máximo innovnum-1)
    unsigned contInnovCon; // Es el número de elementos de la lista conf->listaInnovCon (y también es el máximo innovnum-1)
    TListaInnov *listaInnovCon;	//Mantiene una lista para cada generación de números únicos de innovación para conexiones junto
    // con su nodo de entrada y salida para identificarlos y para poder buscar el innovNum
    // indicado cuando se hace una nueva conexión.
    TListaInnov *listaInnovNodo;	//Mantiene una lista para cada generación de números únicos de innovación para nodos junto
    // con su nodo de entrada y salida de la conexión eliminada para crearlos,para poder buscar el
    // numero del nodo de origen y destino cuando se crea una nueva conexión.
    unsigned numEspecies; // Número de especies actualmente presentes en la conf->población (y tamaño del arreglo de conf->representantes)
    // OPTIMIZACIÓN: Hacer constantes los sizeof de tipos de datos conocidos parq no llamar la funcion sizeof cada vez que se los necesita.
    unsigned *contGeneracSinMejora; // Mantiene una lista de el número de generacones sin mejora en fitness. Se incrementa para cada especie al conf->fInal del
    // ciclo principal, se reinicia a 0 si un genoma se vuelve representante.
    // EXTINCIÓN: si se llega al maxGensSinCambioTh y si la especie No está en el arreglo conservaciónEsp (de tamaño conf->numConservacion),
    // se elimina el representante (de esta lista , de la lista de maxconf->pobPorEspecie y de conf->representantes) y
    // Se reemplazan los individuos con especieEliminada por un genoma inicial totalmente conectado aleatorio (incluyendo maxBias).
    // y se escoge un individuo generado como representante de la especie.
    // Si la especie está en conservación,se disminuye su máximo número de especímenes al 50% del calculado y (maxconf->pobPorEspecie)
    // eliminando aleatoriamente entre el 70% menos apto para dar lugar a más
    //especímenes de las especies nuevas.
    unsigned *conservacionEsp;//Mantiene una lista que contiene los números de especie que se encuentran en conservación se aumenta cuando un especímen pasa a ser
    //representante de una especie al evqaluar su fitness y si supera a maxFitnesConservación.
    unsigned numConservacion;//Es el número de elementos de conservaciónEsp
    unsigned *numGenomasPorEspecie; 	//Mantiene el número de genomas que debe contener cada especie, se renueva en cada generación usando formula de pag 55 de phddissert
    float probMutAN;
    float probMutAC;
    unsigned maxIteraciones;
    float minFitness;
    unsigned maxMemoriaUsada;
    char *fileNameGTDv1;
    char *fileNameSNNv1;
    char *fileNameNNPv1;
    char *fileNameNNPv1Load;
    char *fileNameLog;
    char* fileNameDistrib;
    unsigned repTrain; //repeticiones del archivo de entrenamiento
    unsigned maxBufferSize;
    float c1;
    float c2;
    float c3;
    unsigned eG_t;
    float super;
    float promediarPeso;
    unsigned maxIntentosMutarAC;
    float threshold;
    float porcentVarTh;
    float porcentElim;
    unsigned intentosPareja;
    unsigned **listaOrdenFitness; //crea una matriz de indexpob de genomas pertenecientes a cada especie y la ordena descendentemente por fitness
    //por velocidad, se reserva memoria para conf->sizePob elementos de cada especie).
    unsigned *actNumGenomasPorEspecie; //obtiene el número actual de genomas por gada especie.
    float *fitnessAvgPorEspecie;//para calcular el fitness average de cada especie.
    float minPorcentGenomasPorEspecie;
    unsigned mutacionesPorExterminio;//El número de mutaciones AN y AC que se aplican al genoma inicial de una especie que sustituye a una exterminada
    unsigned maxGeneracSinMejora;//máximo número de generaciones sin mejora que una especie que eno esté en conservación puede exixtir antes de ser reemplazada por un  genoma inicial.
    unsigned maxGenParaNoCruce; // máximo número de generaciones sin mejora en una especie para no realizar cruce.
    unsigned tipoPert; // Tipo de perturbación de pesos 0=uniforme, 1 lineal, 3 cuadrática mayor en la última conex
    unsigned iteracion; //número de la iteración actual, usada por la función actualizarPNodos
    float** dataGTDf; // buffer usado para almacenar los datos de entrenamiento en formato GTD
    tConexDataF** listaConexData; // buffer usado para almacenar las conex  de los genomas en formato SNN para evaluación
    hdrSNNv1* headerSNN; //buffer usado para almacenar el encabezado SNN de cada genoma de la Pob
    float* fitness; //arreglo que guarda los fitness de cada genoma de la Pob
    hdrGTDv1 headerGTD; //usado durante evaluación
    int numDatos; //usado durante evaluación
    int** ordenEval; //orden de evaluación de nodos por genoma (inverso al recursivo empezando en salidas)
    int* tamOrdenEval; //tamaño de ordenEval[i]
    int* tamListaConexPost; // tamaño de la lista de conex después de ordenamiento y simplific.
    int* contO;
    int maxNodos;// número máximo de nodos  por genoma.
    int maxConex; //número máximo de conex, necesario para reservar memoria para los arreglos fijos.
    int realSizePob; // número de genomas totales de la pob incluyendo backups de representantes y genomas auxiliares
    tRandList* randList;
    int tamRandList;
    tRandList* punteroRand;
    float* valoresC; // arreglo  de valores calculados para cálculo de correlación.
    float* cu_dataGTDf; // buffer usado para almacenar los datos de entrenamiento en formato GTD
    tConexDataF** cu_listaConexData; // buffer usado para almacenar las conex  de los genomas en formato SNN para evaluación
    int* cu_tamListaConexPost; // tamaño de la lista de conex después de ordenamiento y simplific.
    float** cu_valorC; // usado para calcular en device el fitness con correlación.
    float** cu_valorTr; // usado para calcular en device el fitness con correlación.
	float* cu_fitness; // usado como vector de retorno del kernel evalPobCUDA
	int** cu_conexIn; //usado para pasar las conexiones a device
	int** cu_conexOut; //usado para pasar las conexiones a device
	float** cu_conexPeso; //usado para pasar las conexiones a device

} TConfig;

//funciones
void help(float version);
/** 	Muestra por stdout la ayuda.
		retorna 0 si error, 1 si ok.
*/
int procParameters(int argc, char *argv[], float version, TConfig* conf);
/** 	Procesa los parámetros de cli.
		retorna 0 si error, 1 si ok.
*/
int inicializaciones(TConfig* conf);
/** 	Inicializa todas las variables de configuración
		retorna 0 si error, 1 si ok.
*/
#endif
