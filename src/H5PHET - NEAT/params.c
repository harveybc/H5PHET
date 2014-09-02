/** Archivo de parámetros globales - C file
 args:
 -t 96 -e8 -n8 -d1 -pD:\\h5phet\\ -xD:\\h5phet\\bringtheshit.bat
*/

#include "params.h"


int procParameters(int argc, char *argv[], float version, TConfig* conf)
{
// procesa los command-line parameters retorna 0 si hay error 1 si ok.
    // si no hay parámetros muestra error.
    char* tmp;
    hdrGTDv1 header;
/*    if (argc==1)
    {
        printf("\nH5PHET - NEAT v%2.2f\nUsando valores por defecto. \nUse: neat -h para obtener ayuda.\n",version);
    }*/
    // hace realloc de cada uno de lsos filenames tamaño strlen() prefix y 2   + 1
    conf->prefix =(char*)calloc(1,300);
    conf->distriCmd = (char*)calloc(1,300);
    tmp=(char*)calloc(1,300);
    conf->fileNameGTDv1=(char*)calloc(1,300);
    conf->fileNameSNNv1=(char*)calloc(1,300);
    conf->fileNameNNPv1=(char*)calloc(1,300);
    conf->fileNameNNPv1Load=(char*)calloc(1,300);
    conf->fileNameLog=(char*)calloc(1,300);
    conf->fileNameDistrib=(char*)calloc(1,300);
    // coloca valores por defecto
    strcpy(conf->fileNameGTDv1,"trdata.gtd\0"); // filename del archivo de entradas (binario)
    strcpy(conf->fileNameSNNv1,"mejor.snn\0"); // filename parael mejor genoma en formato SimpleNeuralNetwork (SNN)
    strcpy(conf->fileNameNNPv1,"reps.nnp\0"); // filename paralos mejores reps en formato NNP
    strcpy(conf->fileNameNNPv1Load,"reps_inic.nnp\0"); // filename para el NNP que se carga a las cargarNNp intentos.
    strcpy(conf->fileNameLog,"log.html\0"); // filename para el archivo de logs en html
    strcpy(conf->fileNameDistrib,"distrib\0"); // inicio del filename de los archivos para procesamiento distribuido.
    strcpy(conf->distriCmd,"D:\\Networking\\bringtheshit.bat\0"); // inicio del filename de los archivos para procesamiento distribuido.

    if ((!tmp)||(!conf->prefix)||(!conf->fileNameGTDv1)||(!conf->fileNameSNNv1)||(!conf->fileNameNNPv1)||(!conf->fileNameNNPv1Load)||(!conf->fileNameLog)||(!conf->fileNameDistrib))
    {
        printf("E0");
        return(0);
    }
    strcpy(conf->prefix,"d:\\h5phet\\\0");
    while ((argc > 1) && (argv[1][0] == '-'))
    {
        switch (argv[1][1])
        {
            // -p<path> prefijo para el path
            case 'p':
                strcpy(conf->prefix,&argv[1][2]);
            break;
            // -x<command> comando del sistema para ejecutar DESPUES de distribProc()
            case 'x':
                strcpy(conf->distriCmd,&argv[1][2]);
            break;
            // -s<unsigned> skip de iteraciones para ejecutar comendo con -x<comando>
            case 's':
                conf->usarDistrib = atoi(&argv[1][2])+1; //porque 0=noUsaDistrib, 1=sinSkip, 2=skip1, 3=skip2, etc...
            break;
            // -d<unsigned> número de archivos de proc distrib a cargar
            case 'd':
                conf->cargarDistrib = atoi(&argv[1][2]); //número de distrib a cargar.
                if (conf->usarDistrib==0) conf->usarDistrib=1;
            break;
            // -t<unsigned> tamaño dew población
            case 't':
                conf->sizePob=atoi(&argv[1][2]);
            break;
            // -e<unsigned> número de especies
            case 'e':
                conf->spEspecies=atoi(&argv[1][2]);
            break;
            // -n<unsigned> número de genomas a cargar de NNP, 0=disabled
            case 'n':
                conf->cargarNNP=atoi(&argv[1][2]);
            break;
            // -h función de ayuda
            case 'h':
                help(version);
                return(0);
            break;
        }
        // fin  switch
        ++argv;
        --argc;
    }// fin while
    // comprueba si hay error en algún param

    // adiciona el filename al prefix.
    // concatena prefix+filenames
    strcpy(tmp,conf->prefix);
    strcpy(tmp+strlen(conf->prefix),conf->fileNameGTDv1);
    strcpy(conf->fileNameGTDv1,tmp);
    strcpy(tmp,conf->prefix);
    strcpy(tmp+strlen(conf->prefix),conf->fileNameSNNv1);
    strcpy(conf->fileNameSNNv1,tmp);
    strcpy(tmp,conf->prefix);
    strcpy(tmp+strlen(conf->prefix),conf->fileNameNNPv1);
    strcpy(conf->fileNameNNPv1,tmp);
    strcpy(tmp,conf->prefix);
    strcpy(tmp+strlen(conf->prefix),conf->fileNameNNPv1Load);
    strcpy(conf->fileNameNNPv1Load,tmp);
     strcpy(tmp,conf->prefix);
    strcpy(tmp+strlen(conf->prefix),conf->fileNameLog);
    strcpy(conf->fileNameLog,tmp);
    strcpy(tmp,conf->prefix);
    strcpy(tmp+strlen(conf->prefix),conf->fileNameDistrib);
    strcpy(conf->fileNameDistrib,tmp);
    // Establece valosres de numEntradas y numSalidas
    // Abre archivo GTDv1
    if((conf->fIn=fopen(conf->fileNameGTDv1,"rb"))==NULL)
    {
        conf->logFile=fopen(conf->fileNameLog,"a+");

        fprintf(conf->logFile,"<br>\nError 57 en funcion evaluarPob() llamando a fopen(%s,\"br\")",conf->fileNameGTDv1);
        return(0);
    }
    // Lee encabezado GTDv1
    if (fread(&header, sizeof(hdrGTDv1), 1, conf->fIn)!=1)
    {
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>Error 58 en función evaluarPob() llamando a fread(%s)\n",conf->fileNameGTDv1);
        return(0);
    }
    conf->numEntradas = header.numEntradas;
    conf->numSalidas = header.numSalidas;
    conf->useFloat = header.tamRegistros==4? 1: 0;
    conf->usarSigned = header.usarSigned==0? 0: 1;
    // imprime encabezado GTDv1
    conf->logFile=fopen(conf->fileNameLog,"w");
    fprintf(conf->logFile,"<br>\n numEntradas=%i, numSalidas=%i, useFloat=%i, usarSigned=%i",conf->numEntradas,conf->numSalidas,conf->useFloat,conf->usarSigned);
    fclose(conf->fIn);
    free(tmp);
    return(1);
}

//Función help, imprime ayuda.
void help(float version)
{
    printf("\nANN V%3.3f\nRed neuronal artificial con entrenador por backpropagation.\nPor Harvey Bastidas.\n\nUso: ann <opciones>\nOpciones:\n",version);
    printf("\n-f usa valores float (32bits), si no se especifica, usa valores double (64bits)\n para archivos de entrada y calculos.");
    printf("\n-s<tipo_activacion> Selecciona el tipo de función de activación ");
    printf("\n <tipo_activacion> puede ser:\n     0= sigma, y=(0,1),            y = 1 / (1 + exp (- D * x))");
    printf("\n     1 = sigma aprox, y=(0,1),      y = 0.5 + x * (1 - abs(x) / 2), y=0 si x<=-1, y=1 si x>=1");
    printf("\n     2 = tanh, y=(-1,1),            y = 2 / (1 + exp(-2 * x)) - 1");
    printf("\n     3 = gauss, y=(-1,1),           y = exp(- x * x)");
    printf("\n     4 = elliot, y=(-1,1),          y = x / (1 + |x|)");
    printf("\n(def)5 = elliot unitario, y=(0,1 ), y = (x / 2) / (1 + |x|) + 0.5 \n");
    printf("-d<numero> valor de D solo para  tSigma=0, por defecto <numero>=1\n");
    printf("-b<filename> nombre de archivo bin32/64 (matriz cuadrada de pesos)\n");
    printf("-t<iteraciones> modo de entrenamiento , especifica el número \n de iteraciones,n  necesita un archivo de entradas (opción -e) y un archivo de \n genoma (opción -g) produce un archivo de pesos de salida (opción -m) \n");
    printf("-e<filename> archivo de entradas, debe contener vectores de \n 32/64bits (ver opción -f) del mismo tamaño que la matriz de pesos de entrada o\n el especificado en el genoma.\n");
    printf("-g<filename> archivo de genoma (formato por definir)\n");
    printf("-o<filename> archivo de salida (en formato por definir.)\n");
    printf("-m<filename> archivo de salida con matriz de pesos. (en formato por definir.)\n");
    printf("-i<numero> número de neuronas en la capa de entrada. \n");
    printf("-n<numero> número de neuronas en la capa de salida. \n");
    printf("-w<neuronas> fuerza creación de matriz para feedforward de <neuronas> neuronas.\n");
    printf("-c<capas> especifica número de capas ocultas default=1.\n");
    printf("-r<numero> umbral de activación (generalmente 0 para func de activación (-1,1) o 0.5 para (0,1) ). \n");
    printf("-h imprime ayuda.\n");
}

int inicializaciones(TConfig* conf)
{
// inicializa todas las variables de configuración, retorna 0 si error, 1 si ok
// TODO: colocar como const todos los parámetros no dinámicos para acelerar ejecución.
// TODO: renombrar params.h como datastr.h y hacer un nuevo config.h que se incluye en todos lso demás
//      en el cual se declaran como constantes todos los parámetros y los valores iniciales de variables.
//      para poder usarlos de forma GLOBAL dentro de las funciones.
//      VERIFICAR que lso nombres globales no correspondan con nombres de variables locales
    conf->tamRandList=1000000; //(10M?) tamaño de lista de números aleatorios pre-generados
    conf->maxConex = 100000; //(1M?) Usado para reservar memoria el principio para la lista de ordenEval
    conf->maxNodos = 10000; //(100k?) Usado para reservar memoria el principio para la lista de ordenEval
    conf->usarDistrib = 0; //(1) 0=No usa proc.distribuido. 1=Usa en todas las iteraciones, 2= salta una iterac, etc...
    conf->cargarDistrib = 0; //(2) número de NNPs numerados a leer para procesamiento distribuido
    conf->cargarNNP=0; //(0) número máximo de genomas a cargar desde NNP después de alcanzar spEspecies
    conf->probInterSp = 0.00001; // (0.001) probabilid por genoma de que se produzca un cruce entre especies diferentes
    conf->porcentCompetencia = -0.01   ; // probabilidad de que gane el mejor indexpob entre dos generaciones (solo de la misma especie).
    conf->maxIntentosDistInicial = 0; //número de intentos para obtener genoma incicial lo más heterogeneo posible.
    // (default=10) 0= nu usa maxdistancia máximo número de intentos para máxima distancia de genomas iniciales
    conf->maxIntentosDist= 0; //máximo (default=5)número de intentos para máxima distancia con cada renew de genoma.
    conf->maxGeneracSinMejoraPG= 800; // (default=8)máximo número de generaciones sin mejora antes de generarmaxdist
    conf->probMutMaxDist= -0.02; //probabilidad de mutación por máximas distancias por cada genoma

    conf->numOcultas = 0;			//Especifica el número de neuronas ocultas (se debe actualizar cada vez que se haga un nuevo gen)
    conf->numBias = 1;				//Especifica el número de neuronas de entrada de bias (cte 1)
    conf->Fthreshold = 1;			//Especifica el threshold, debe ser 1
    //TODO: verificar si es necesario que todos los fSigma(-1) = 0 (o -1 en salida con signo) y fSigma(0)=0
    //TODO: verificar si es mejor usar entradas y salidas con signo.
    conf->tSigma = 1003; 				// (default=4 gauss) Tipo de función de activación mayores a 1000 = rango de salida(-1,1) sino (0,1)
    conf->fSigmaD = 13; 			//Por defecto usa un D=9 para la funcion sigma Normal float(0)
    conf->sizePob = 96; 			//(36)Tamaño de la conf->población.F (mejor 20 y spEspecies=20)
    conf->spEspecies = 8; //(3) SQRT(sizepob)????			//Set Pounsigned para el número de especies que se debe tratar de mantener usando c_t
//Verificar ESto, para que siempre sean 3 individuos mínimo por especie.
    conf->minPorcentGenomasPorEspecie = 0.0002; //0.0002Es el mínimo de genomas por especie (COMO porcentaje de sizepob) que debe tener cada especie a pesar de tener mala fitness para evitar que haya especies con numGenomas=0
    conf->maxIteraciones = 500000;		//(500k) Máximo número de itraciones para terminar el ciclo principal
    conf->mutacionesPorExterminio = 0;	// (default 0) Número de mutaciones que se hacen al genoma randomizado después de exterminio.
    conf->maxGenParaNoCruce = 10; 	// (default 6) número máximo de generaciones sin mejora en especie para parar reproducción (como lo hace Stanley) solo mutan AC+AN y pesos.
    conf->maxGeneracSinMejora = 10000;		// (default 12) máximo número de generaciones sin mejora antes de exterminio de especie.
    conf->probMutPeso = 0.8; 			// (default=0.8) Probabilidad de mutación de pesos de conex de cada Genoma
    // TODO : verificar si se necesita este parámetro, además verificar uno por uno los parámetros a ver si se pueden calcular desde otro o si se necesitan
    conf->maxDesvThEspecies = 0.001; 	// (default=0.01 hasta 10 especies es ele 10% de desviación) Máxima desviación desde threshold de compat max 1 min 0
    conf->probMutTh = 0; //Probabilidad de mutación de threshold en cada genoma.
    conf->probMutAN = 0.01;	//(0.01) Probabilidad de mutación AN durante ciclo principal
    // razon de complejidad nodos/conex
    conf->probMutAC = 0.03; //(0.01)Probabilidad de mutación AC durante ciclo principal
    conf->threshold = 0.8; //Thrreshold de compatibilidad de especies usado para diferenciar entre una especie y otra
    conf->tipoPert = 0; // Tipo de perturbación de pesos 0=uniforme, 1 lineal, 2 cuadrática mayor en la última conex
    //TODO: verificar si funciona bién pertpesos quitando mutaciones y selección.
    conf->porcentMutPeso = 0.01; // (CAMBIADO default=0.01) Porcentaje aumentado a una mutaciónaleatoria normal
    conf->porcentMutNTh = 0.05; //porcentaje aumentado a una mutación de th
    conf->minFitness = 1; //Mínimo fitness necesario para terminar el ciclo principal
    conf->maxMemoriaUsada = 800; //Máximo tamaño en memoria para terminar el ciclo principal
    conf->porcentEnableds = 0.002;
    conf->antMejor = 0;
    conf->porcentElim = 0.5; //Porcentaje de población que se  reduce durante el proceso de selección para cada especie.
    conf->minConexMutPeso = 0; //
    conf->maxEspeciesConservacion = 0; //Es el máximo número de especies en conservación. si este número se supera debido a buén fitness de otras especies,
    conf->numConservacion = 0; //Es el número de elementos de conservaciónEsp
    conf->pesoMinInicial = -1.0; //peso mínimo para valores aleatorios asignados en la primera generación
    conf->pesoMaxInicial = 1.0; //peso máximo para valores aleatorios asignados en la primera generación.
    conf->pob = NULL; //Arreglo conf->población, es un arreglo de genomas.
    conf->representantes = NULL; //Arreglo de indexpob de conf->representantes para cada especie.
    conf->ultimaInnov = 0; // Guarda el valor de la última innovación (acumulativa entre nodos y conexiones)
    conf->fIn = NULL; // Puntero a archivo de entradas
    conf->fOut = NULL; // Puntero a archivo de entradas
    conf->contInnovNodo = 0; // Es el número de elementos de la lista conf->listaInnovNodo (y también es el máximo innovnum-1)
    //conf->contInnovCon; // Es el número de elementos de la lista conf->listaInnovCon (y también es el máximo innovnum-1)
    conf->actualizarEnCambios = 0; //0 solo actualiza una vez en ciclo ppal después de evaluarPob antes de realizar cruce , >0 evalúa en cada cambio de la población.
    conf->listaInnovCon = NULL;	//Mantiene una lista para cada generación de números únicos de innovación para conexiones junto
    // con su nodo de entrada y salida para identificarlos y para poder buscar el innovNum
    // indicado cuando se hace una nueva conexión.
    conf->listaInnovNodo = NULL; //Mantiene una lista para cada generación de números únicos de innovación para nodos junto
    // con su nodo de entrada y salida de la conexión eliminada para crearlos,para poder buscar el
    //numero del nodo de origen y destino cuando se crea una nueva conexión.
    conf->numEspecies = 0; //Número de especies actualmente presentes en la conf->población (y tamaño del arreglo de conf->representantes)
    //OPTIMIZACIÓN: Hacer constantes los sizeof de tipos de datos conocidos parq no llamar la funcion sizeof cada vez que se los necesita.
    conf->contGeneracSinMejora = NULL; //Mantiene una lista de el número de generacones sin mejora en fitness. Se incrementa para cada especie al conf->fInal del
    // ciclo principal, se reinicia a 0 si un genoma se vuelve representante.
    // EXTINCIÓN: si se llega al maxGensSinCambioTh y si la especie No está en el arreglo conservaciónEsp (de tamaño conf->numConservacion),
    // se elimina el representante (de esta lista , de la lista de maxconf->pobPorEspecie y de conf->representantes) y
    // Se reemplazan los individuos con especieEliminada por un genoma inicial totalmente conectado aleatorio (incluyendo maxBias).
    // y se escoge un individuo generado como representante de la especie.
    // Si la especie está en conservación,se disminuye su máximo número de especímenes al 50% del calculado y (maxconf->pobPorEspecie)
    // eliminando aleatoriamente entre el 70% menos apto para dar lugar a más
    //especímenes de las especies nuevas.
    conf->conservacionEsp = NULL; // Mantiene una lista que contiene los números de especie que se encuentran en conservación se aumenta cuando un especímen pasa a ser
    //representante de una especie al evqaluar su 1220435 fitness y si supera a maxFitnesConservación.
    conf->numGenomasPorEspecie = NULL; // Mantiene el número de genomas que debe contener cada especie, se renueva en cada generación usando formula de pag 55 de phddissert
/*
    conf->prefix = "d:\\h5phet\\\0";
    conf->fileNameGTDv1 = "trdata.gtd\0"; // filename del archivo de entradas (binario)
    conf->fileNameSNNv1 = "mejor.snn\0"; // filename parael mejor genoma en formato SimpleNeuralNetwork (SNN)
    conf->fileNameNNPv1 = "reps.nnp\0"; // filename paralos mejores reps en formato NNP
    conf->fileNameNNPv1Load = "reps_inic.nnp\0"; // filename para el NNP que se carga a las cargarNNp intentos.
    conf->fileNameLog = "log.html\0"; // filename para el archivo de logs en html
    conf->fileNameDistrib ="distrib\0"; // inicio del filename de los archivos para procesamiento distribuido.

    */
    // credirebaja
    //conf->fileNameGTDv1 = "/srv/www/htdocs/phppgadmin/h5phet/training/trdata.gtd\0"; // filename del archivo de entradas (binario)
    //conf->fileNameSNNv1 = "/srv/www/htdocs/phppgadmin/h5phet/test_mejor.snn\0"; // filename parael mejor genoma en formato SimpleNeuralNetwork (SNN)
    //conf->fileNameNNPv1 = "/srv/www/htdocs/phppgadmin/h5phet/test_reps.nnp\0"; // filename paralos mejores reps en formato NNP
    //conf->fileNameNNPv1Load = "/srv/www/htdocs/phppgadmin/h5phet/reps_inic.nnp\0"; // filename para el NNP que se carga a las cargarNNp intentos.
    //conf->fileNameLog = "/srv/www/htdocs/phppgadmin/h5phet/test_log.html\0"; // filename para el archivo de logs en html
    //conf->fileNameDistrib ="/srv/www/htdocs/phppgadmin/h5phet/distrib\0";

    // credirebaja
    //conf->fileNameGTDv1 = "/var/www/phppgadmin/h5phet/training/trdata.gtd\0"; // filename del archivo de entradas (binario)
    //conf->fileNameSNNv1 = "/var/www/phppgadmin/h5phet/test_mejor.snn\0"; // filename parael mejor genoma en formato SimpleNeuralNetwork (SNN)
    //conf->fileNameNNPv1 = "/var/www/phppgadmin/h5phet/test_reps.nnp\0"; // filename paralos mejores reps en formato NNP
    //conf->fileNameNNPv1Load = "/var/www/phppgadmin/h5phet/reps_inic.nnp\0"; // filename para el NNP que se carga a las cargarNNp intentos.
    //conf->fileNameLog = "/var/www/phppgadmin/h5phet/test_log.html\0"; // filename para el archivo de logs en html
    //conf->fileNameDistrib ="/var/www/phppgadmin/h5phet/distrib\0";


    // Formato H5PHET-SNN: int numIns,numBias,numOuts,numHiddens,int numConex, int conexIn[numConex],conexOut[numConex],double peso[numConex]
    //                     el valor de los elementos de los arreglos es el index de cada Nodo In u out de la conex,
    //                     el orden de los valores es nodosIn,nodosBias,nodosOut,nodosHidden.
    conf->repTrain = 0; //repeticiones del archivo de entrenamiento durante cada fase de evaluación
    conf->maxBufferSize = 70000000; //tamaño (en floats) del buffer de lectura
    conf->c1 = 1; //factor de proporcionalidad en la formula de distancia para disjounsigned  ints
    conf->c2 = 1; //factor de proporcionalidad en la formula de distancia para excess
    conf->c3 = 0.4; //orig 0.4 factor de proporcionalidad en la formula de distancia para promedio de diferencia de pesos
    conf->eG_t = 100000; //threshold para considerar número de genes "excesivamente" grandes (y hacer N=1 en la formula de distancia)
    conf->super = -0.1; //Probabilidad de usar excess y disjounsigned del menos apto además de los que se heredan normalmente del más apto durante cruce
    conf->promediarPeso = 0.01; //Probabilidad de promediar el peso en caso de matching Genes durante cruce
    conf->maxIntentosMutarAC = 3; //Máximo número de busqueda aleatoria de conexiones no existentes para realizar mutaciónAC
    conf->porcentVarTh = 0.05; //Taza de variación del threshold para alcanzar el número de especies requerido
    conf->intentosPareja = 4; //unsigned  intentos de búsqueda aleatoria de pareja antes de que se recurra a busqueda secuencial decremental en fitness
    // Lee el header del archivo GTDv1 para sobreescribir valores de numEntradas y numSalidas
    return(1);
}
