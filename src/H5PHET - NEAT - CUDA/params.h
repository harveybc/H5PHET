/** Archivo de par�metros globales - H file

*/
#ifndef PARAMS_H_INCLUDED
#define PARAMS_H_INCLUDED

#include <stdlib.h> //para usar calloc, realloc, free, rand, etc...
#include <stdio.h> //Para usar printf, etc...
#include <time.h> //para usar clock() para hacer funcion delay.
#include <math.h> //Para funciones matem�ticas (exp,abs, otras.)
#include <string.h> //para usar memset(), memmove, memcpy()...
#include <limits.h> //Para verificar los l�mites de unsigned y rand
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
 	double lastFitness; // usado para programaci�n distribuida.
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

// es una colecci�n de SNNs.
// formato NNPv1: [hdrNNPv1][SNNv1_0][SNNv1_1]...[SNNv1_n]
// donde SNNv1_x = red neuronal simple en formato SNNv1.
typedef struct
{
 	char fileID[3]; // Debe ser siempre NNP para identificar el tipo de archivo.
 	unsigned char version;
 	unsigned numGenomas;
} hdrNNPv1;

typedef struct   // tipo GenConexF gen de conexi�n.
{
    unsigned innovNum; //es el n�mero de la conexi�n(puede quitarse si se usa el mimso indice del arreglo de nodos).
    unsigned nodoIn; //innovnum del Nodo de origen de la conexi�n.
    unsigned nodoOut; //innovnum del Nodo de destino de la conexi�n.
    float peso; //Peso de la conexi�n.
    char recurrente; //0=No recurrente, 1=Recurrente.
    char enabled; //0=Desactivado, 1=Activado.
    unsigned indexIn; //index del nodo de entrada en el arreglo de nodos. Necesarios porque al hacer nuevoNodo se hace realloc y los punteros quedan apuntando a un valor libre.
    unsigned indexOut; //index del nodo de salida en el arreglo de nodos. Se utilizan en evaluarGenoma, se modifican en nuevonodo y nuevaconex
    //TODO: Cambiar todos los short unsigned por short short unsigned  int. para ahorrar memoria.
} GenConexF;

typedef struct   // tipo GenNodoF gen de nodo. Para SMP se necesita adem�s una lista de los nodos padres y un contador de padres.
{
    float thNodo; //Threshold de activaci�n de la neurona.
    unsigned innovNum;//Obtenido al crear nuevo nodo buscando lista de ;
    short unsigned nodeFunction; //0=entrada,1=oculto,2=salida,3=bias.
    float maxTh; // Mayor valor posible de th para esta neurona, es la sumatoria de los pesos positivos que entran
    float minTh; // Menor valor posible de th para esta neurona, es la sumatoria de los pesos negativos	que entran
    float valor; //Valor de salida de cada nodo(para computar la ann sin matriz de pesos).
    //TODO: hacer realloc de arreglos de indexes cuando se agregue nodo para que quede en
    GenConexF** conexHijo; //arreglo de punteros a conexiones a cada hijo solo se usa antes de la eval.
    unsigned contHijos; //contador de hijos. Se usan solo en evaluarGenoma.
    short unsigned estadoC; // flag que dice si el valor No ha sido calculado (0), est� calculandose (1) o ya est� calculado (2)
    //TODO: Falta arreglo de punteros a nodos hijos, arreglo de puntero a conexi�n a hijo y contador de nodos hijos y en genoma nueva func recursiva calcularValorNodo()
} GenNodoF;

typedef struct   // tipo Genoma es un genoma que contiene lista de nodos y conexiones,etc...
{
    GenNodoF *nodo; //Arreglo din�mico de nodos.
    GenConexF *conex; //Arreglo din�mico de conexiones.
    unsigned especie; //Especie a la que pertenece el genoma en la conf->poblaci�n. Se establece por funcion calcularEspecie() despu�s de la segunda generaci�n.
    unsigned totalNodos;//Mantiene el total de nodos en el gen (Maxindex+1)
    unsigned totalConexiones;//Mantiene el total de conexiones (Maxindex+1)
    unsigned maxInnovNumConex; //Mantiene el m�imo n�mero de innovaci�n de genes de conexi�n.
    unsigned maxInnovNumNodo; //Mantiene el m�ximo n�mero de  innovaci�n para genes de nodo.
    int numHijos; //Es el n�mero de hijos que debe tener cada genoma por cruce. (se inicializa en -1, se modifica con)
    float fitness;//Usado para almacenar el error del genoma respecto a la salida deseada.
} Genoma;

typedef struct   // tipo TNodoOut de nodos de salida con su respectivo innovNum usado en TListaInnov.
{
    unsigned innovNum;//N�mero de innovaci�n correspondiente a la salida relacionada con el index (entrada).
    unsigned nodoOut; //
} TNodoOut;

typedef struct   // tipo TListaInnov para listas de innovaciones, la lista es un puntero a esta struct. el index es el n�mero del nodo de entrada.
{
    unsigned numOut; // Para cada index(in), es el n�mero disponible de salidas (para limitar el index de este arreglo).
    TNodoOut *nodoOut;// Es un arreglo de estructuras listOut que contienen el innovNum y el nodo de salida de la conexi�n.
} TListaInnov;
//TODO: pasar a cada Funcion ex�ctamente los par�metros que necesita en lugar de todo el conf.

typedef struct   // tipo TConfig de par�metros de NEAT
{
    int usarDistrib; // flag de uso de procesamiento distribuido (0=no, 1= todas las iter, 2=saltando una, etc...)
 	char* distriCmd; // comando para ejecutar despu�s de distribProc
 	char* prefix; // prefijo para los filenames
    unsigned tmpIndexPob; // usado para guardar la posici�n del indexpob usado como genoma temporal (ultimo del arreglo pob).
    int cargarDistrib;
    FILE *logFile; //usado para handler de archivo de logs.
    float probInterSp; // (default=0.001) Probabilidad de que en una iteraci�n se produzca reproducc�n inter-especies.
    unsigned cargarNNP; // (default=0) 1=Carga genomas iniciales desde archivo NNP (NeuralNetworkPopulation)
    float porcentCompetencia; // (default=0.8) -1=no usada. especifica la prob de que el mejor genoma gane entre generaciones (para cada genoma)
    unsigned maxIntentosDistInicial; // (default=10) 0= nu usa maxdistancia m�ximo n�mero de intentos para m�xima distancia de genomas iniciales
    unsigned maxIntentosDist; //m�ximo (default=5)n�mero de intentos para m�xima distancia con cada renew de genoma.
    unsigned maxGeneracSinMejoraPG; // (default=8)m�ximo n�mero de generaciones sin mejora antes de generarmaxdist
    float probMutMaxDist; //probabilidad de mutaci�n por m�ximas distancias por cada genoma
	unsigned usarSigned; // 1 si se usar�n valores con signo para entrenamiento, 0 si se usan valores sin signo
    unsigned numEntradas;//Especifica el m�mero de neuronas de entrada
    unsigned numSalidas;//Especifica el n�mero de neuronas de salida
    unsigned numOcultas;//Especifica el n�mero de neuronas ocultas (se debe actualizar cada vez que se haga un nuevo gen)
    unsigned numBias;//Especifica el n�mero de neuronas de entrada de bias (cte 1)
    float Fthreshold;//Especifica el threshold, debe ser 1 para sigma que retornan [0,1] y -1,1
    unsigned tSigma; //Por defecto usa funcion de activaci�n elliot unitario
    float fSigmaD; //Por defecto usa un D=1 para la funcion sigma Normal float(0)
    float A,D,F; // coeficientes de sigma calculados a mano para que f(0.5)=0.5 f(0)=1 y f(-1)=0 en tsigma (0,1) y f(0.5)=0 y f(-1)=-1 para tsigma (-1,1)
    unsigned useFloat;//indica si se deben usar todos los valores de entrada/salida/calculos tipo float si useFloat=0 se usan tipo double(lento).
    unsigned sizePob; //Tama�o de la conf->poblaci�n.
//	float *c_t ; //compatibility threshold (m�nima distancia para nueva especie.), //TODO, corregir valor inicial basado en ejemplos citados al conf->fInal del archivo
    unsigned spEspecies; //Set Pounsigned para el n�mero de especies que se debe tratar de mantener usando c_t
    float maxFitnessConservacion; ////TODO: hacer funcion para pasar de error a fitness teniendo como par�metro Fitness paa error=100, linealidad: 0=>1-e, 1=>1/x, 2=1/x^2, etc.
    unsigned maxEspeciesConservacion; //Es el m�ximo n�mero de especies en conservaci�n. si este n�mero se supera debido a bu�n fitness de otras especies,
    //a las especies excedentes(de menor fitness) se les d� el n�mero de espec�menes calculado sin restricciones
    //TODO: falta calcular autom�ticamente los valores iniciales (rango para escoger aleatorio) de pesos de conexiones a nodos de salidas usando prean�lisis de los valores
    // de salida, los thresholds de las neuronas de entrada se pueden calcular aleatoriamente en rangos determinados por los valores m�ximos y m�nimos
    // de entrada que toma una muestra de entrenamiento. Todas las conexiones adicionales tendr�n peso en rango de -1 a 1 o de 0 a 1 dependiendo
    //de tSigma.
    float pesoMinInicial; //peso m�nimo para valores aleatorios asignados en la primera generaci�n
    float pesoMaxInicial; //peso m�ximo para valores aleatorios asignados en la primera generaci�n.
    Genoma *pob; //Arreglo conf->poblaci�n, es un arreglo de genomas.
    unsigned *representantes; //Arreglo de indexpob de conf->representantes para cada especie.
    unsigned ultimaInnov; //Guarda el valor de la �ltima innovaci�n (acumulativa entre nodos y conexiones)
    FILE *fIn;//Puntero a archivo de entradas
    FILE *fOut;//Puntero a archivo de entradas
    char *fileNameIn; //nombre del archivo de entradas
    float porcentEnableds;
    float antMejor;
    unsigned minConexMutPeso;
    float maxDesvThEspecies;  //0 hasta 1 max desviaci�n desde threshold de compatibilidad.
    float probMutPeso;
    float probMutTh;
    float porcentMutPeso;
    float porcentMutNTh;

    //TODO: separar en estructura control todas las variables globales usadas por el prog y en conf las de configuraci�n para simplificar config.
    unsigned actualizarEnCambios; //0 solo actualiza una vez en ciclo ppal despu�s de evaluarPob antes de realizar cruce , >0 eval�a en cada cambio de la poblaci�n.
    unsigned contInnovNodo; // Es el n�mero de elementos de la lista conf->listaInnovNodo (y tambi�n es el m�ximo innovnum-1)
    unsigned contInnovCon; // Es el n�mero de elementos de la lista conf->listaInnovCon (y tambi�n es el m�ximo innovnum-1)
    TListaInnov *listaInnovCon;	//Mantiene una lista para cada generaci�n de n�meros �nicos de innovaci�n para conexiones junto
    // con su nodo de entrada y salida para identificarlos y para poder buscar el innovNum
    // indicado cuando se hace una nueva conexi�n.
    TListaInnov *listaInnovNodo;	//Mantiene una lista para cada generaci�n de n�meros �nicos de innovaci�n para nodos junto
    // con su nodo de entrada y salida de la conexi�n eliminada para crearlos,para poder buscar el
    // numero del nodo de origen y destino cuando se crea una nueva conexi�n.
    unsigned numEspecies; // N�mero de especies actualmente presentes en la conf->poblaci�n (y tama�o del arreglo de conf->representantes)
    // OPTIMIZACI�N: Hacer constantes los sizeof de tipos de datos conocidos parq no llamar la funcion sizeof cada vez que se los necesita.
    unsigned *contGeneracSinMejora; // Mantiene una lista de el n�mero de generacones sin mejora en fitness. Se incrementa para cada especie al conf->fInal del
    // ciclo principal, se reinicia a 0 si un genoma se vuelve representante.
    // EXTINCI�N: si se llega al maxGensSinCambioTh y si la especie No est� en el arreglo conservaci�nEsp (de tama�o conf->numConservacion),
    // se elimina el representante (de esta lista , de la lista de maxconf->pobPorEspecie y de conf->representantes) y
    // Se reemplazan los individuos con especieEliminada por un genoma inicial totalmente conectado aleatorio (incluyendo maxBias).
    // y se escoge un individuo generado como representante de la especie.
    // Si la especie est� en conservaci�n,se disminuye su m�ximo n�mero de espec�menes al 50% del calculado y (maxconf->pobPorEspecie)
    // eliminando aleatoriamente entre el 70% menos apto para dar lugar a m�s
    //espec�menes de las especies nuevas.
    unsigned *conservacionEsp;//Mantiene una lista que contiene los n�meros de especie que se encuentran en conservaci�n se aumenta cuando un espec�men pasa a ser
    //representante de una especie al evqaluar su fitness y si supera a maxFitnesConservaci�n.
    unsigned numConservacion;//Es el n�mero de elementos de conservaci�nEsp
    unsigned *numGenomasPorEspecie; 	//Mantiene el n�mero de genomas que debe contener cada especie, se renueva en cada generaci�n usando formula de pag 55 de phddissert
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
    unsigned *actNumGenomasPorEspecie; //obtiene el n�mero actual de genomas por gada especie.
    float *fitnessAvgPorEspecie;//para calcular el fitness average de cada especie.
    float minPorcentGenomasPorEspecie;
    unsigned mutacionesPorExterminio;//El n�mero de mutaciones AN y AC que se aplican al genoma inicial de una especie que sustituye a una exterminada
    unsigned maxGeneracSinMejora;//m�ximo n�mero de generaciones sin mejora que una especie que eno est� en conservaci�n puede exixtir antes de ser reemplazada por un  genoma inicial.
    unsigned maxGenParaNoCruce; // m�ximo n�mero de generaciones sin mejora en una especie para no realizar cruce.
    unsigned tipoPert; // Tipo de perturbaci�n de pesos 0=uniforme, 1 lineal, 3 cuadr�tica mayor en la �ltima conex
    unsigned iteracion; //n�mero de la iteraci�n actual, usada por la funci�n actualizarPNodos
    float** dataGTDf; // buffer usado para almacenar los datos de entrenamiento en formato GTD
    tConexDataF** listaConexData; // buffer usado para almacenar las conex  de los genomas en formato SNN para evaluaci�n
    hdrSNNv1* headerSNN; //buffer usado para almacenar el encabezado SNN de cada genoma de la Pob
    float* fitness; //arreglo que guarda los fitness de cada genoma de la Pob
    hdrGTDv1 headerGTD; //usado durante evaluaci�n
    int numDatos; //usado durante evaluaci�n
    int** ordenEval; //orden de evaluaci�n de nodos por genoma (inverso al recursivo empezando en salidas)
    int* tamOrdenEval; //tama�o de ordenEval[i]
    int* tamListaConexPost; // tama�o de la lista de conex despu�s de ordenamiento y simplific.
    int* contO;
    int maxNodos;// n�mero m�ximo de nodos  por genoma.
    int maxConex; //n�mero m�ximo de conex, necesario para reservar memoria para los arreglos fijos.
    int realSizePob; // n�mero de genomas totales de la pob incluyendo backups de representantes y genomas auxiliares
    tRandList* randList;
    int tamRandList;
    tRandList* punteroRand;
    float* valoresC; // arreglo  de valores calculados para c�lculo de correlaci�n.
    float* cu_dataGTDf; // buffer usado para almacenar los datos de entrenamiento en formato GTD
    tConexDataF** cu_listaConexData; // buffer usado para almacenar las conex  de los genomas en formato SNN para evaluaci�n
    int* cu_tamListaConexPost; // tama�o de la lista de conex despu�s de ordenamiento y simplific.
    float** cu_valorC; // usado para calcular en device el fitness con correlaci�n.
    float** cu_valorTr; // usado para calcular en device el fitness con correlaci�n.
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
/** 	Procesa los par�metros de cli.
		retorna 0 si error, 1 si ok.
*/
int inicializaciones(TConfig* conf);
/** 	Inicializa todas las variables de configuraci�n
		retorna 0 si error, 1 si ok.
*/
#endif
