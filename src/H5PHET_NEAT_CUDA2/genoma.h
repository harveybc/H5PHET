/** Redes neuronales para NEAT - H file
	usan/modifican genomas únicamente (pueden usarse operaciones de genes de gen.h en un genoma).
*/
#ifndef GENOMA_H_INCLUDED
#define GENOMA_H_INCLUDED
#ifndef PARAMS_H_INCLUDED
#include "params.h"
#endif

float evalGenom(int indexpob,TConfig* conf);
/** Evalúa un genoma con los datos de dataGTDf, retorna el fitness.
*/

unsigned guardarGenomaSNN(unsigned indexpob, char *filename, TConfig* conf);
/** 	escribe el genoma indexpob de pob en un archivo.
//El formato de salida es(sin separadores): Genoma, Genoma.nodo, genoma.conex las longitudes a escribir
//de cada estructura se basan en el tamaño de Genoma, GenNodoF, GenconexF y en los valores Genoma.totalNodos
//y Genoma.totalConexiones
//Parámetros:	indexpob = indice del arreglo de genomas conf->pob que se va a guardar
//				filename = path y nombre de archivo en el que se guardará el genoma
//Retorna 0 si hay error, 1 si ok
*/

int snnDataLoader(const unsigned indexpob, const unsigned especie, FILE* fileIn, TConfig* conf);
/** carga un SNN en indexpob desde un archivo fileIn que debe estar abierto reconstruyendo el genoma.
*/

unsigned actualizarPNodos(unsigned index, TConfig* conf);
/**     actualiza los valores de los punteros a conexiones hijo para cada nodo de un genoma.
//      se debe llamar al inicio de evaluarGenoma.
//      retorna 0 si hay error, 1 si ok
*/
float calcularValorNodo(Genoma* pGenoma, int indexPob, GenNodoF* nodo, int indexNodo, TConfig* conf);//OPTIMIZADA, se puede dejar otra Funcion que calcule el valor solo con inidce en la struct GenNodoF para ahorrar memoria(que puede ser limitada) pero sacrificando velocidad.
/**     calcula recursivamente valores de salida (después de pasar por fsigma) de un nodo de un genoma indicado.
//      retorna el valor del nodo calculado, también lo asigna a valor y coloca en 1 valorCalculado.
//      parámetros:  indexPob = index del genoma en la población.
//      nodo = puntero al nodo para el cual se quiere calcular el valor.
*/
unsigned copiarGenoma(unsigned srcindexpob,unsigned dstindexpob, TConfig* conf);
/** 	copia el genoma de origen al genoma de destino (el de destino es borrado)
//SUPONE que ya se ha reservado memoria para el array conf->pob que incluye a laso dos elementos: funcion inicializarPob();
//Retorna 0 si hay error, 1 si Ok
	//copia el contenido de conf->pob[src] al contenido de  conf->pob[dest] (copia un Genoma)
	////TODO: verificar que todos los free se hagan antes de asignar los nuevos valores de los punteros
*/
unsigned crossover(unsigned indexpob1, unsigned indexpob2, unsigned indexpobOut ,float super, float promediarPob,float porcentEnableds, TConfig* conf);
/** 	realiza el cruce entre dos genomas dados sus indexpob y lo coloca en un elemento de conf->pob (puede ser uno de los padres)
//Toma calcula el fitness de los dos genomas, se heredan los matching genes randomly,
//Los disjounsigned y excess se heredan solo del fittest,
//Parámetros: 	indexpob1, indexpob2 = genomas a cruszar
//				indexpobOut = indexpob donde se debe colocar el genoma resultante del cruce.
//				super = float entre 0 y 1, Probabilidad de heredar los exess y disjounsigned  ints del menos apto (aparte de los que se heredan normalmente del más apto)
//				promendiarProb = float entre 0 y 1 = probabilidad de que en caso de matching, se promedien los pesos en lugar de
								//seleccionarlos aleatoriamente entre los padres.
//Retona 0 si hubo error, 1 si ok, coloca en la variable global tempGenoma el genoma generado por crossover
*/
unsigned evaluarGenoma(int indexPob, Genoma* pGenoma, unsigned primero, float *entradas, float *salidas, TConfig* conf);
/** 	para un genoma i obtiene los valores y los fitness NO ajustados(1-error) para cada NODO de un genoma (incluyendo las salidas)
//también acumula fitness en pob[index].fitness para ser procesado luego por evaluarPob
//Parámetros:
//				index = indice de genoma por parámetro y
//				primero = especifica si es la primera vez (=1) que se evalúa el genoma para inicializar si no es la primera vez=0
//				entradas = puntero a un arreglo de nEntradas valores float que serán las entradas a evaluar
//				salidas = puntero a un arreglo de nSalidas valores float que serán las salidas deseadas, respecto a las cuales se obtendrá el error y por tant el fitness = 1-error.
//				nEntradas = número de elementos en el arreglo entradas
//				nSalidas = número de elementos en el arreglo de salidas
//sus valores en 0 excepto los de las entradas.
//salida: retorna el genoma en el indice index evaluado para la entrada con la variable valor de la estructura GenNodoF evaluada.
//retorna 0 si hubo error.
//DEMORADA, probar velocidad haciendo primero funcion crearmatrix(genoma) y evaluarmatrix(), en lugar de evaluargenoma (puede ser otro parámetro de esta funcion)
////TODO:Colocarle nuevo parámetro para evaluar directamente(actual) o evaluar por generación de matrix y evaluación de matrix (como funciones?).
//contadores kjkjj*/
unsigned genomaPerfecto(unsigned index, TConfig* conf);
/** 	Coloca un genoma perfecto de xor en la posición deseada.
//parametros :index
*/
unsigned genomaInicial(unsigned index,unsigned nEntradas, unsigned nSalidas, unsigned nBias, unsigned primero, unsigned especie,TConfig* conf);
/** 	crea un nuevo genoma totalmente conectado asigna la especie 0 y lo ubica en el índice 0 de la conf->|población.
//parámetros: index=donde queda el genoma inicial, nEntradas, nSalidas, nBias.
//			primero si =1 se bora lista de innovaciones , representantes, conservacion, generaciones sin mejora, etc...
//retorna 0 si hubo error, 1 si la creación fué exitosa.
*/
unsigned mutarAC(unsigned indexpob, unsigned maxIntentos, TConfig* conf);
/** 	agrega una conexión al azar usando la funcion nuevaconexión
//adiciona la nueva conex.
//también realiza la asignación del innovNum de la conexión buscando en la lista(mediante in y out), sino existe, lo adiciona.
//Verifica que la conexión entre los nodos resultantes no exista.
//Se debe además verificar que en casos de conexiones recurrentes, el nodo no sea de entrada o bias.
//Se verifica si la conexión seleccionada ya existe, si esto ocurre randomiza de  nuevo la entrada y salida (en mutar AC)
//hasta un número máximo maxIntentosNuevacon para evitar que se entre en bucle inconf->fInito.
//si no encuentra una nueva conexión posible, retorna 0 pero no crea la conexión.
//y al haber deswcubierto solo conexiones redundantes, incrementa aleatoriamente el peso de la última conexión encontrada??
//Si se crea satisfactoriamemnte la conexión se le asigna un peso de 1
//La especie se asigna con la funcion calcularEspecie para la primera generación, luego se recalcula después de cada mutación+cruce.
//Parámetros: indexpob = indice de la conf->población del genoma a mutar.
//Retorna 0 si hay error, 1 si ok, 2 si no se pudo agregar por límite de conexiones.
*/
unsigned mutarAN(unsigned indexpob, TConfig* conf);
/** 	agrega un nodo al azar usando la funcion nuevoNodo, entran en la selección todas las conexiones existentes.
//reerva memoria para el nuevo tamaño del genoma con realloc y luego adiciona el nuevo nodo.
//también asigna el número de nodo relizando  la busqueda (mediante in y out) en la lista de innovaciones de nodo para ver si ya existe, sino, la adiciona a la  lista.
//La busqueda podría acelerarse con if anidado en lugar de if (in and out)
//también podría acelerarse si se hace un arreglo de estructuras  con index=in y estructuras (unsigned maxOut,unsigned *structout )  y en cada
//structout un out y un innovnum
//así solo se escogería la entrada y en los outs (sol hasta maxouts)se buscaría el out deseado y el innovnum.
//entonces se buscaría únicamente en maxouts en lugar de entre todas las innovaciones, reduciendo drásticamente el tiempo de búsqueda..
//Pero al crearse cada nodo debería actualizar mazout y hacere realloc para un nuevo structout.
//Y tqambién se usaría más memoria,antes para cada innov se usaban 12bytes(8 usando innovNum como index) en lista rápida
// se usan 20 bytes para la primera entrada y 8 para el resto. :) I can live with that.
//¿Es posible agregar un nodo bias?

//Retorna 0 si hubo error , 1 si ok.
//Parámetros:	indexpob	= index de conf->pob que se desea mutar.
*/
unsigned nuevaConex(unsigned indexpob, unsigned indexIn, unsigned indexOut, float peso, short unsigned recurrente, short unsigned enabled, TConfig* conf);
/** 	adiciona una nueva conexión a un genoma de la conf->población (param = indice del genoma en la conf->pob, todas
//las variables de struct GenConexF excepto innovNum que es una variable global para cada gen.).
//Parámetros: indexpob, IndexIn, IndexOut,recurrente, peso, enabled.
//Retorna 0 si hubo error, 1 si ok.
////TODO MODIFICAR LA FUNCION EVALUAR GENOMA PARA EVALUAR POR INNOVNUM DE NODOS y CONEXIONES en lugar de indexes.
	//obtener memoria con realloc en conf->pob[indexpob].conex para un tamaño (totalConexiones+1)sizeof(GenNodoF)
	*/
unsigned nuevoNodo(unsigned indexpob, unsigned indexElimInnovConex, TConfig* conf);
/** 	adiciona una nueva conexión a un genoma en la conf->población (param=indice de genoma en la conf->pob, todas las vars de la struct
//GenConexF excepto el innovNum que se adiciona automáticamente al último entre toda la conf->población (debe llevarse var global de ultInnov))
//Funcionamiento: nuevoNodo() escoge una conexión y la elimina,agrega un nodo y agragega un gen de conexión
//entre el (indice del)nodo de origen de la actual y el nuevo pero=1 y agrega otro gen de conexión entre el nuevo y el anterior de destino
//con peso=peso anterior, requiere numeros de innovacion del nodo y las conexiones resultantes y el indice de la conexión eliminada y de conf->pob.
//Parámetros indiceconf->pob,IndiceConexiónInnvovEliminada(en conf->pob[ip].nodo[]),funcion (1=oculta,3=bias),valorNuevoNodo,innovNumNuevoNodo,innovNumConexI_N,InnovNumN_O.
//Retorna 0 si hubo algún error, 1 si ok.
	//Coloca en disabled la conexión con número de innovación indexElimConex*/
#endif //GENOMA_H_INCLUDED

void genomaMasLejano(int indexPob,int intentos, TConfig* conf);
/**     busca el genoma más lejano a todos los de su especie entre la población  <i/tentos> veces y lo coloca en indexPob
// también debe tener distancia menor a 1/2 de distacia a especie más cercana.
*/
