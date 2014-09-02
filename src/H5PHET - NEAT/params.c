/** Archivo de par�metros globales - C file
 args:
 -t 96 -e8 -n8 -d1 -pD:\\h5phet\\ -xD:\\h5phet\\bringtheshit.bat
*/

#include "params.h"


int procParameters(int argc, char *argv[], float version, TConfig* conf)
{
// procesa los command-line parameters retorna 0 si hay error 1 si ok.
    // si no hay par�metros muestra error.
    char* tmp;
    hdrGTDv1 header;
/*    if (argc==1)
    {
        printf("\nH5PHET - NEAT v%2.2f\nUsando valores por defecto. \nUse: neat -h para obtener ayuda.\n",version);
    }*/
    // hace realloc de cada uno de lsos filenames tama�o strlen() prefix y 2   + 1
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
            // -d<unsigned> n�mero de archivos de proc distrib a cargar
            case 'd':
                conf->cargarDistrib = atoi(&argv[1][2]); //n�mero de distrib a cargar.
                if (conf->usarDistrib==0) conf->usarDistrib=1;
            break;
            // -t<unsigned> tama�o dew poblaci�n
            case 't':
                conf->sizePob=atoi(&argv[1][2]);
            break;
            // -e<unsigned> n�mero de especies
            case 'e':
                conf->spEspecies=atoi(&argv[1][2]);
            break;
            // -n<unsigned> n�mero de genomas a cargar de NNP, 0=disabled
            case 'n':
                conf->cargarNNP=atoi(&argv[1][2]);
            break;
            // -h funci�n de ayuda
            case 'h':
                help(version);
                return(0);
            break;
        }
        // fin  switch
        ++argv;
        --argc;
    }// fin while
    // comprueba si hay error en alg�n param

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
        fprintf(conf->logFile,"<br>Error 58 en funci�n evaluarPob() llamando a fread(%s)\n",conf->fileNameGTDv1);
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

//Funci�n help, imprime ayuda.
void help(float version)
{
    printf("\nANN V%3.3f\nRed neuronal artificial con entrenador por backpropagation.\nPor Harvey Bastidas.\n\nUso: ann <opciones>\nOpciones:\n",version);
    printf("\n-f usa valores float (32bits), si no se especifica, usa valores double (64bits)\n para archivos de entrada y calculos.");
    printf("\n-s<tipo_activacion> Selecciona el tipo de funci�n de activaci�n ");
    printf("\n <tipo_activacion> puede ser:\n     0= sigma, y=(0,1),            y = 1 / (1 + exp (- D * x))");
    printf("\n     1 = sigma aprox, y=(0,1),      y = 0.5 + x * (1 - abs(x) / 2), y=0 si x<=-1, y=1 si x>=1");
    printf("\n     2 = tanh, y=(-1,1),            y = 2 / (1 + exp(-2 * x)) - 1");
    printf("\n     3 = gauss, y=(-1,1),           y = exp(- x * x)");
    printf("\n     4 = elliot, y=(-1,1),          y = x / (1 + |x|)");
    printf("\n(def)5 = elliot unitario, y=(0,1 ), y = (x / 2) / (1 + |x|) + 0.5 \n");
    printf("-d<numero> valor de D solo para  tSigma=0, por defecto <numero>=1\n");
    printf("-b<filename> nombre de archivo bin32/64 (matriz cuadrada de pesos)\n");
    printf("-t<iteraciones> modo de entrenamiento , especifica el n�mero \n de iteraciones,n  necesita un archivo de entradas (opci�n -e) y un archivo de \n genoma (opci�n -g) produce un archivo de pesos de salida (opci�n -m) \n");
    printf("-e<filename> archivo de entradas, debe contener vectores de \n 32/64bits (ver opci�n -f) del mismo tama�o que la matriz de pesos de entrada o\n el especificado en el genoma.\n");
    printf("-g<filename> archivo de genoma (formato por definir)\n");
    printf("-o<filename> archivo de salida (en formato por definir.)\n");
    printf("-m<filename> archivo de salida con matriz de pesos. (en formato por definir.)\n");
    printf("-i<numero> n�mero de neuronas en la capa de entrada. \n");
    printf("-n<numero> n�mero de neuronas en la capa de salida. \n");
    printf("-w<neuronas> fuerza creaci�n de matriz para feedforward de <neuronas> neuronas.\n");
    printf("-c<capas> especifica n�mero de capas ocultas default=1.\n");
    printf("-r<numero> umbral de activaci�n (generalmente 0 para func de activaci�n (-1,1) o 0.5 para (0,1) ). \n");
    printf("-h imprime ayuda.\n");
}

int inicializaciones(TConfig* conf)
{
// inicializa todas las variables de configuraci�n, retorna 0 si error, 1 si ok
// TODO: colocar como const todos los par�metros no din�micos para acelerar ejecuci�n.
// TODO: renombrar params.h como datastr.h y hacer un nuevo config.h que se incluye en todos lso dem�s
//      en el cual se declaran como constantes todos los par�metros y los valores iniciales de variables.
//      para poder usarlos de forma GLOBAL dentro de las funciones.
//      VERIFICAR que lso nombres globales no correspondan con nombres de variables locales
    conf->tamRandList=1000000; //(10M?) tama�o de lista de n�meros aleatorios pre-generados
    conf->maxConex = 100000; //(1M?) Usado para reservar memoria el principio para la lista de ordenEval
    conf->maxNodos = 10000; //(100k?) Usado para reservar memoria el principio para la lista de ordenEval
    conf->usarDistrib = 0; //(1) 0=No usa proc.distribuido. 1=Usa en todas las iteraciones, 2= salta una iterac, etc...
    conf->cargarDistrib = 0; //(2) n�mero de NNPs numerados a leer para procesamiento distribuido
    conf->cargarNNP=0; //(0) n�mero m�ximo de genomas a cargar desde NNP despu�s de alcanzar spEspecies
    conf->probInterSp = 0.00001; // (0.001) probabilid por genoma de que se produzca un cruce entre especies diferentes
    conf->porcentCompetencia = -0.01   ; // probabilidad de que gane el mejor indexpob entre dos generaciones (solo de la misma especie).
    conf->maxIntentosDistInicial = 0; //n�mero de intentos para obtener genoma incicial lo m�s heterogeneo posible.
    // (default=10) 0= nu usa maxdistancia m�ximo n�mero de intentos para m�xima distancia de genomas iniciales
    conf->maxIntentosDist= 0; //m�ximo (default=5)n�mero de intentos para m�xima distancia con cada renew de genoma.
    conf->maxGeneracSinMejoraPG= 800; // (default=8)m�ximo n�mero de generaciones sin mejora antes de generarmaxdist
    conf->probMutMaxDist= -0.02; //probabilidad de mutaci�n por m�ximas distancias por cada genoma

    conf->numOcultas = 0;			//Especifica el n�mero de neuronas ocultas (se debe actualizar cada vez que se haga un nuevo gen)
    conf->numBias = 1;				//Especifica el n�mero de neuronas de entrada de bias (cte 1)
    conf->Fthreshold = 1;			//Especifica el threshold, debe ser 1
    //TODO: verificar si es necesario que todos los fSigma(-1) = 0 (o -1 en salida con signo) y fSigma(0)=0
    //TODO: verificar si es mejor usar entradas y salidas con signo.
    conf->tSigma = 1003; 				// (default=4 gauss) Tipo de funci�n de activaci�n mayores a 1000 = rango de salida(-1,1) sino (0,1)
    conf->fSigmaD = 13; 			//Por defecto usa un D=9 para la funcion sigma Normal float(0)
    conf->sizePob = 96; 			//(36)Tama�o de la conf->poblaci�n.F (mejor 20 y spEspecies=20)
    conf->spEspecies = 8; //(3) SQRT(sizepob)????			//Set Pounsigned para el n�mero de especies que se debe tratar de mantener usando c_t
//Verificar ESto, para que siempre sean 3 individuos m�nimo por especie.
    conf->minPorcentGenomasPorEspecie = 0.0002; //0.0002Es el m�nimo de genomas por especie (COMO porcentaje de sizepob) que debe tener cada especie a pesar de tener mala fitness para evitar que haya especies con numGenomas=0
    conf->maxIteraciones = 500000;		//(500k) M�ximo n�mero de itraciones para terminar el ciclo principal
    conf->mutacionesPorExterminio = 0;	// (default 0) N�mero de mutaciones que se hacen al genoma randomizado despu�s de exterminio.
    conf->maxGenParaNoCruce = 10; 	// (default 6) n�mero m�ximo de generaciones sin mejora en especie para parar reproducci�n (como lo hace Stanley) solo mutan AC+AN y pesos.
    conf->maxGeneracSinMejora = 10000;		// (default 12) m�ximo n�mero de generaciones sin mejora antes de exterminio de especie.
    conf->probMutPeso = 0.8; 			// (default=0.8) Probabilidad de mutaci�n de pesos de conex de cada Genoma
    // TODO : verificar si se necesita este par�metro, adem�s verificar uno por uno los par�metros a ver si se pueden calcular desde otro o si se necesitan
    conf->maxDesvThEspecies = 0.001; 	// (default=0.01 hasta 10 especies es ele 10% de desviaci�n) M�xima desviaci�n desde threshold de compat max 1 min 0
    conf->probMutTh = 0; //Probabilidad de mutaci�n de threshold en cada genoma.
    conf->probMutAN = 0.01;	//(0.01) Probabilidad de mutaci�n AN durante ciclo principal
    // razon de complejidad nodos/conex
    conf->probMutAC = 0.03; //(0.01)Probabilidad de mutaci�n AC durante ciclo principal
    conf->threshold = 0.8; //Thrreshold de compatibilidad de especies usado para diferenciar entre una especie y otra
    conf->tipoPert = 0; // Tipo de perturbaci�n de pesos 0=uniforme, 1 lineal, 2 cuadr�tica mayor en la �ltima conex
    //TODO: verificar si funciona bi�n pertpesos quitando mutaciones y selecci�n.
    conf->porcentMutPeso = 0.01; // (CAMBIADO default=0.01) Porcentaje aumentado a una mutaci�naleatoria normal
    conf->porcentMutNTh = 0.05; //porcentaje aumentado a una mutaci�n de th
    conf->minFitness = 1; //M�nimo fitness necesario para terminar el ciclo principal
    conf->maxMemoriaUsada = 800; //M�ximo tama�o en memoria para terminar el ciclo principal
    conf->porcentEnableds = 0.002;
    conf->antMejor = 0;
    conf->porcentElim = 0.5; //Porcentaje de poblaci�n que se  reduce durante el proceso de selecci�n para cada especie.
    conf->minConexMutPeso = 0; //
    conf->maxEspeciesConservacion = 0; //Es el m�ximo n�mero de especies en conservaci�n. si este n�mero se supera debido a bu�n fitness de otras especies,
    conf->numConservacion = 0; //Es el n�mero de elementos de conservaci�nEsp
    conf->pesoMinInicial = -1.0; //peso m�nimo para valores aleatorios asignados en la primera generaci�n
    conf->pesoMaxInicial = 1.0; //peso m�ximo para valores aleatorios asignados en la primera generaci�n.
    conf->pob = NULL; //Arreglo conf->poblaci�n, es un arreglo de genomas.
    conf->representantes = NULL; //Arreglo de indexpob de conf->representantes para cada especie.
    conf->ultimaInnov = 0; // Guarda el valor de la �ltima innovaci�n (acumulativa entre nodos y conexiones)
    conf->fIn = NULL; // Puntero a archivo de entradas
    conf->fOut = NULL; // Puntero a archivo de entradas
    conf->contInnovNodo = 0; // Es el n�mero de elementos de la lista conf->listaInnovNodo (y tambi�n es el m�ximo innovnum-1)
    //conf->contInnovCon; // Es el n�mero de elementos de la lista conf->listaInnovCon (y tambi�n es el m�ximo innovnum-1)
    conf->actualizarEnCambios = 0; //0 solo actualiza una vez en ciclo ppal despu�s de evaluarPob antes de realizar cruce , >0 eval�a en cada cambio de la poblaci�n.
    conf->listaInnovCon = NULL;	//Mantiene una lista para cada generaci�n de n�meros �nicos de innovaci�n para conexiones junto
    // con su nodo de entrada y salida para identificarlos y para poder buscar el innovNum
    // indicado cuando se hace una nueva conexi�n.
    conf->listaInnovNodo = NULL; //Mantiene una lista para cada generaci�n de n�meros �nicos de innovaci�n para nodos junto
    // con su nodo de entrada y salida de la conexi�n eliminada para crearlos,para poder buscar el
    //numero del nodo de origen y destino cuando se crea una nueva conexi�n.
    conf->numEspecies = 0; //N�mero de especies actualmente presentes en la conf->poblaci�n (y tama�o del arreglo de conf->representantes)
    //OPTIMIZACI�N: Hacer constantes los sizeof de tipos de datos conocidos parq no llamar la funcion sizeof cada vez que se los necesita.
    conf->contGeneracSinMejora = NULL; //Mantiene una lista de el n�mero de generacones sin mejora en fitness. Se incrementa para cada especie al conf->fInal del
    // ciclo principal, se reinicia a 0 si un genoma se vuelve representante.
    // EXTINCI�N: si se llega al maxGensSinCambioTh y si la especie No est� en el arreglo conservaci�nEsp (de tama�o conf->numConservacion),
    // se elimina el representante (de esta lista , de la lista de maxconf->pobPorEspecie y de conf->representantes) y
    // Se reemplazan los individuos con especieEliminada por un genoma inicial totalmente conectado aleatorio (incluyendo maxBias).
    // y se escoge un individuo generado como representante de la especie.
    // Si la especie est� en conservaci�n,se disminuye su m�ximo n�mero de espec�menes al 50% del calculado y (maxconf->pobPorEspecie)
    // eliminando aleatoriamente entre el 70% menos apto para dar lugar a m�s
    //espec�menes de las especies nuevas.
    conf->conservacionEsp = NULL; // Mantiene una lista que contiene los n�meros de especie que se encuentran en conservaci�n se aumenta cuando un espec�men pasa a ser
    //representante de una especie al evqaluar su 1220435 fitness y si supera a maxFitnesConservaci�n.
    conf->numGenomasPorEspecie = NULL; // Mantiene el n�mero de genomas que debe contener cada especie, se renueva en cada generaci�n usando formula de pag 55 de phddissert
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
    conf->repTrain = 0; //repeticiones del archivo de entrenamiento durante cada fase de evaluaci�n
    conf->maxBufferSize = 70000000; //tama�o (en floats) del buffer de lectura
    conf->c1 = 1; //factor de proporcionalidad en la formula de distancia para disjounsigned  ints
    conf->c2 = 1; //factor de proporcionalidad en la formula de distancia para excess
    conf->c3 = 0.4; //orig 0.4 factor de proporcionalidad en la formula de distancia para promedio de diferencia de pesos
    conf->eG_t = 100000; //threshold para considerar n�mero de genes "excesivamente" grandes (y hacer N=1 en la formula de distancia)
    conf->super = -0.1; //Probabilidad de usar excess y disjounsigned del menos apto adem�s de los que se heredan normalmente del m�s apto durante cruce
    conf->promediarPeso = 0.01; //Probabilidad de promediar el peso en caso de matching Genes durante cruce
    conf->maxIntentosMutarAC = 3; //M�ximo n�mero de busqueda aleatoria de conexiones no existentes para realizar mutaci�nAC
    conf->porcentVarTh = 0.05; //Taza de variaci�n del threshold para alcanzar el n�mero de especies requerido
    conf->intentosPareja = 4; //unsigned  intentos de b�squeda aleatoria de pareja antes de que se recurra a busqueda secuencial decremental en fitness
    // Lee el header del archivo GTDv1 para sobreescribir valores de numEntradas y numSalidas
    return(1);
}
