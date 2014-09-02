/** H5PHET - NEAT
	Implementación en ANSI C 89 de la técnica: "Neuro Evolution of Augmenting Topologies" de Kenneth Stanley (2004).
	Usa archivo de entrada GTD (Generic Training Data) y genera archivo de salida SNN (Simple Neural Network).
	Parte de H5PHET.
	Por Harvey Bastidas.
*/
# include <sys/types.h>
//# include <sys/local.h>
#include "params.h"
#include "auxiliares.h"
#include "genoma.h"
#include "pob.h"


//TODO: pasar código por gnuindent para estandarizar la indentación.
////TODO: Funcion freeAll() libera memoria de todos los arreglos usados. :) TODO: problemas:
//		hay conexiones con el mismo innovnum repetidas en un genoma.
//		hay disminución de fitness entre evaluaciones
//		hay error de punteros cuando numbias=0
//		en la iteración 130 de un experimento se retornó fitness de 8.4
//		hay que adicionar como parámetro uso de mutar_tSigma y su respectiva probabilidad por genoma para mutar el tSigma de un nod, también correcciones en calcularValorNodo();
//		perturbación de pesos con variación máxima para el último mínima para nodos de entrada y bias y lineal o cuadrática para el unsigned  intérvalo
//		implementar mutación removerNodo y remover conexión (verificar si es posible y en que casos y condiciones debería usarse)
//      erroro de punteros después de e my raro(gen 100+)
//      error poco frecuente en we012 (Mirar donde y a que index de cone se hace malloc y ver porque falla al hacer free)
//          también verificar otras variables que se liberan
//      TODO: cálculo de distáncia mín para pertenecer a especie a partir de spEspecies , numDisjounsigned (o mejor porcent respecto más grande?), numExcess y difPEsos (al menos con esto se puede hacer) max y min como rángo de búsqueda
//      función verificarMejor retorna que se decrementó el fitness después de una extinción.
//		TODO: es probable que evaluargenoma no funcione bién, probarlo con genomaPerfecto.
//		TODO: función para evaluar un genoma en particular.diferente a evaluarGenoma porque usa evaluarpob para el streaming de entradas
//      TODO: después de corregir bug de we012 y errores de crossover hacer variables: maxThreadProc, usarCUDA y numThCuda (mutex?)
//              además colocar

/********* MAIN *********/
int main(int argc, char *argv[])
{

    TConfig conf; //parámetros del sistema

//Versión:

    float version=0.71; // H5PHET 0.71 Parte de EVA parte de TGV
                        // Entrada = Archivo GTD y NNP de inicio o para procesamiento distribuido
                        // Salida = Archivo SNN y NNP.

    unsigned j=1;
    unsigned k=1;
    unsigned i;
    unsigned theOne=0;
    //(genera un número entre 0 y 1e3 )
    i=(unsigned  int) time(NULL);
    j=j*(unsigned  int)(1+floor((i-10.0*floor((float)i/10.0))+0.5));
    j=j*(unsigned  int)(floor((i-1000.0*floor(i/1000.0))+0.5));
    j=j+(unsigned  int)floor((i-10.0*floor((float)i/10.0))+0.5)+1;
    j=j+(unsigned  int)floor((i-1000.0*floor(i/1000.0))+0.5);
    //Inicializa el pseudo-random seed con la hora actual.
    srand(j);
  /*  for (i=0; i<j; i++)
    {
        k=(unsigned )rand();
    }
	*/
	k=(unsigned )rand();
    // re-Inicializa el pseudo-random seed con el último rand obtenido.
    srand(k);
    // inicializa todas las variables de configuración a valores por defecto.
    inicializaciones(&conf);
	//verifica y asigna los parámetros de linea de comandos, si hay algún error, sale.
	if (!procParameters(argc,argv,version, &conf))
	{
	    printf("\nError M6 en main() llamando a procParameters\n");
        return(0);
	}
    // abre archivo de logs.
    // unsigned primeraGen(float tamPob, unsigned nEntradas, unsigned nSalidas, unsigned nBias, float minPeso, float maxPeso, unsigned maxMutacionesAN,unsigned maxMutacionesAC,unsigned maxIntentosMutAC, short unsigned useMutarAC, short unsigned useMutarAN,unsigned useRandomization){
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nH5PHET - NEAT V%2.2f\nCreando primera generación\n",version);
    if ((primeraGen(conf.sizePob,conf.numEntradas,conf.numSalidas,conf.numBias,conf.pesoMinInicial,conf.pesoMaxInicial,conf.mutacionesPorExterminio,conf.mutacionesPorExterminio,conf.maxIntentosMutarAC,0,0,1,&conf))==0) // TODO: Cuadrar en params.h los parámetros useMutarAN y useMutarAC de esta Funcion.
    {
        fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nError M1 en main() en función primeraGen(%u,%u,%u,%u,%1.1f,%1.1f,%u,%u,%u,%u,%u,%u)",conf.sizePob,conf.numEntradas,conf.numSalidas,conf.numBias,-3.0,3.0,3,1,100,0,0,1);
        fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nGuardando mejor genoma en c:\\mejorGenoma.txt");
        i=guardarGenomaSNN(buscarMejorFitness(&conf),conf.fileNameSNNv1,&conf); //TODO: esto es necesario en este punto?
        fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nSaliendo...");
        return(0);
    }
//ciclo principal
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>Iniciando ciclo principal\n");
    if ((cicloPrincipal(&conf))==UINT_MAX)
    {
        fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nError M2 en main() en funcion cicloPrincipal()");
        fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nGuardando mejor genoma en c:\\mejorGenoma.txt");
        i = guardarGenomaSNN(buscarMejorFitness(&conf),conf.fileNameSNNv1,&conf);
        imprimirGenoma(i,&conf);
        fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nSaliendo...");
        fclose(conf.fOut);
        fclose(conf.fIn);
        return(0);
    }

    theOne=buscarMejorFitness(&conf);

    /*	fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>unsigned  introduciendo GenomaPerfectocomo elemento 15");
    	if (genomaPerfecto(15,&conf)==0){ //necesaria esta evaluación antes de primera especiación?
    		fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nError 59 en funcion cicloPrincipal() llamando a evaluarPob()");
    		return(UINT_MAX);
    	}
    theOne=15;*/
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nLista de Representantes: ");
    for (i=0; i<conf.numEspecies; i++)
    {
        fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nEspecie %u(%u) = %7.7f",i,conf.representantes[i],conf.pob[conf.sizePob+i].fitness);
    }
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nGuardando mejor genoma %u con fitness = %11.11f en c:\\mejorGenoma.txt",theOne,conf.pob[theOne].fitness);
    i=guardarGenomaSNN(theOne,conf.fileNameSNNv1,&conf);

/*
    //Prueba la Funcion xor con el mejor genoma
    entr=(float*)malloc(sizeof(float)*8);
    if (entr==NULL)
        return(0);
    salid=(float*)malloc(sizeof(float)*4);
    if (salid==NULL);
    {
        free (entr);
        free (salid);//??para evitar warning de clocwork
        return(0);
    }
    entr[0]=-1;
    entr[1]=-1;
    entr[2]=-1;
    entr[3]=1;
    entr[4]=1;
    entr[5]=-1;
    entr[6]=1;
    entr[7]=1;
    salid[0]=1;
    salid[1]=1;
    salid[2]=1;
    salid[3]=1;
    actualizarPNodos(theOne,&conf);
    i=evaluarGenoma(theOne,0,entr,salid,&conf);
    salid[0]=conf.pob[theOne].nodo[3].valor;
    i=evaluarGenoma(theOne,0,entr+2,salid+1,&conf);
    salid[1]=conf.pob[theOne].nodo[3].valor;
    i=evaluarGenoma(theOne,0,entr+4,salid+2,&conf);
    salid[2]=conf.pob[theOne].nodo[3].valor;
    i=evaluarGenoma(theOne,0,entr+6,salid+3,&conf);
    salid[3]=conf.pob[theOne].nodo[3].valor;

    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nEvaluación de la tabla de XOR para el mejor genoma: \n");
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br> %3.3f XOR %3.3f = %3.3f \n",entr[0],entr[1],salid[0]);
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br> %3.3f XOR %3.3f = %3.3f \n",entr[2],entr[3],salid[1]);
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br> %3.3f XOR %3.3f = %3.3f \n",entr[4],entr[5],salid[2]);
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br> %3.3f XOR %3.3f = %3.3f \n",entr[6],entr[7],salid[3]);
*/

    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nMemoria Usada=%lu\n",calcularMemoriaUsada(conf.sizePob,&conf)/1024);
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nNúmero de especies iniciales: %u\n",conf.numEspecies);
    imprimirGenoma(theOne,&conf);
    fclose(conf.logFile);conf.logFile=fopen(conf.fileNameLog,"a+");fprintf(conf.logFile,"<br>\nH5PHET - NEAT finanizado correctamente.\n");
    ////TODO funcion freeAll(), libera la memoria requerida por cada arreglo de nodos y conexiones de cada genoma para toda la conf.población y libera conf.fInalmente *conf.pob
    return(0);
}

//Mashauritaki
/*
VER EVOLVING NEUR... Para parametros de XOR y //TODOntes de perturbación de pesos.
Inicializaciones de parámetros en pag 14 de Evolving Neural networks through augmenting topologies.pdf
La especición se hace basado en la dt el threshold de distancia mínimo para una nueva especie.




*/
