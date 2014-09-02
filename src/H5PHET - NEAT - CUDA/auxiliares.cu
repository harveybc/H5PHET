/** Funciones genéticas auxiliares (busqueda, diagnóstico, etc...) - C file
	funciones que no modifican el genoma ni la población (estadísticas, búsquedas).
*/

#ifndef PARAMS_H_INCLUDED
#include "params.h"
#define PARAMS_H_INCLUDED
#endif
#include "auxiliares.h"



float correlac(float* Vc,TConfig* conf)
// retorna el coeficiente de correlación de dos arreglos de datos x e y de tamaño numDatos.
// la matriz de correlaciones correlacM debe inicializarse cn valores >10 la primera vez.
// parámetros:
//      x,y = vectores a comparar
//      xm,ym = media de x y media de y, se deben proporcionar por optimización de pre-calculo de medias
//      numDatos = tamaño de los vectores
{
    int i=0;
    float xm=0;
    float ym=0;
    float sum0=0;
    float sum1=0;
    float sum2=0;
    // calcula las medias
    for (i=0;i<conf->numDatos;i++)
    {
        xm+=conf->dataGTDf[i][conf->headerGTD.numEntradas];
        ym+=Vc[i];
    }
    xm/=conf->numDatos;
    ym/=conf->numDatos;
    // calcula sum0(0,nD,(Xi-Xm)*(Yi-Ym)),sum1(0,nD,sqr(Xi-Xm)) y sum2(0,nD,sqr(Yi-Ym)
    for (i=0;i<conf->numDatos;i++)
    {
        sum0+=((conf->dataGTDf[i][conf->headerGTD.numEntradas]-xm)*(Vc[i]-ym));
        sum1+=((conf->dataGTDf[i][conf->headerGTD.numEntradas]-xm)*(conf->dataGTDf[i][conf->headerGTD.numEntradas]-xm));
        sum2+=((Vc[i]-ym)*(Vc[i]-ym));
    }
    // retorna la correlación
    return(sum0/(sqrt(sum1)*sqrt(sum2)));
}


void genSSNhdr1(int indexpob, TConfig * conf)
{
    conf->headerSNN[indexpob].fileID[0] = 'S';
    conf->headerSNN[indexpob].fileID[1] = 'N';
    conf->headerSNN[indexpob].fileID[2] = 'N';
    conf->headerSNN[indexpob].version = 1;
    conf->headerSNN[indexpob].usarSigned = conf->tSigma>=1000? 1 : 0;
    conf->headerSNN[indexpob].tamRegistros = conf->useFloat==0? 8 : 4;
    conf->headerSNN[indexpob].numEntradas = conf->numEntradas;
    conf->headerSNN[indexpob].numSalidas = conf->numSalidas;
    conf->headerSNN[indexpob].numBias = conf->numBias;
    conf->headerSNN[indexpob].numHiddens = conf->pob[indexpob].totalNodos-(conf->numEntradas+conf->numSalidas+conf->numBias);
    conf->headerSNN[indexpob].numConex = conf->pob[indexpob].totalConexiones;
    conf->headerSNN[indexpob].sigmaFactor = (double)conf->A;
    conf->headerSNN[indexpob].actThreshold = (double) conf->Fthreshold;
    conf->headerSNN[indexpob].lastFitness = (double) conf->pob[indexpob].fitness; // usado para programac distribuida
}

float randL(TConfig* conf)
{
    float valor;
    valor=(*conf->punteroRand).valor;
    conf->punteroRand=(tRandList*)(*conf->punteroRand).next;
    return(valor);
}

int genOrdenEvalF1g(int indexpob, TConfig* conf)
{
    int j,temp;
    int maxN=0;
    char* valCalculado;
    // inicializa los campos del header de SNN
    conf->headerSNN[indexpob].fileID[0] = 'S';
    conf->headerSNN[indexpob].fileID[1] = 'N';
    conf->headerSNN[indexpob].fileID[2] = 'N';
    conf->headerSNN[indexpob].version = 1;
    conf->headerSNN[indexpob].usarSigned = conf->tSigma>=1000? 1 : 0;
    conf->headerSNN[indexpob].tamRegistros = conf->useFloat==0? 8 : 4;
    conf->headerSNN[indexpob].numEntradas = conf->numEntradas;
    conf->headerSNN[indexpob].numSalidas = conf->numSalidas;
    conf->headerSNN[indexpob].numBias = conf->numBias;
    conf->headerSNN[indexpob].numHiddens = conf->pob[indexpob].totalNodos-(conf->numEntradas+conf->numSalidas+conf->numBias);
    conf->headerSNN[indexpob].numConex = conf->pob[indexpob].totalConexiones;
    conf->headerSNN[indexpob].sigmaFactor = (double)conf->A;
    conf->headerSNN[indexpob].actThreshold = (double) conf->Fthreshold;
    conf->headerSNN[indexpob].lastFitness = (double) conf->pob[indexpob].fitness; // usado para programac distribuida
    //calcula el número de nodos
    maxN=conf->headerSNN[indexpob].numEntradas+conf->headerSNN[indexpob].numBias+conf->headerSNN[indexpob].numSalidas+conf->headerSNN[indexpob].numHiddens;
    // libera memoria de listaConexData[indexpob]
    free(conf->listaConexData[indexpob]);
    // reserva memoria para listaConexData[indexpob],tam=numConex[indexpob]*sizeof(tConexDataF)
    conf->listaConexData[indexpob]=(tConexDataF*)malloc((conf->pob[indexpob].totalConexiones)*sizeof(tConexDataF));
//printf("%d,",indexpob);
    if (!conf->listaConexData[indexpob])
    {
        printf("\nError 66.7 en genOrdenEvalF1g llamando a malloc()");
        return(0);
    }
//TODO: FALTA VER PORQUE NO SE PUDO LIBERAR?    // libera memoria de ordenEval
// TODO: HACER FUNCION EVALGENOM() para evaluar solo1
// TODO: guardar después de evaalgenom para verificar y setear el leastfitness.
    valCalculado=(char*)malloc(maxN*sizeof(char));
    if (!valCalculado)
    {
        printf("\nError 61.1 en genOrdenEvalF1g llamando a malloc()");
        return(0);
    }
    // inicializa en 0 valcalculado para todos excepto bias
    for (j=0; j<maxN; j++)
    {
        valCalculado[j]=0;
    }
    // inicializa valcalculado en 1 para bias
    for (j=0; j<conf->headerSNN[indexpob].numBias; j++)
    {
        valCalculado[conf->headerSNN[indexpob].numEntradas+j]=1;
    }
    for (j=0; j<conf->pob[indexpob].totalConexiones; j++)
    {
        // genera listaConexData[indexpob][j]=tConexDataF[j]
        conf->listaConexData[indexpob][j].conexIn=conf->pob[indexpob].conex[j].indexIn;
        conf->listaConexData[indexpob][j].conexOut=conf->pob[indexpob].conex[j].indexOut;
        conf->listaConexData[indexpob][j].enabled=conf->pob[indexpob].conex[j].enabled;
        conf->listaConexData[indexpob][j].peso=conf->pob[indexpob].conex[j].peso;
    }
    // inicializa el contador global de posición
    conf->tamOrdenEval[indexpob]=0;
    temp=conf->headerSNN[indexpob].numEntradas+conf->headerSNN[indexpob].numBias;
    // para cada salida
    for(j=0; j<conf->headerSNN[indexpob].numSalidas; j++)
    {
        // llama recurOrganicer(indexNodoOut, tConexDataF* listaConexData, int* ordenEval)
        recurOrganicer(temp+j, conf->headerSNN[indexpob] ,conf->listaConexData[indexpob],  conf->ordenEval[indexpob], valCalculado,&(conf->tamOrdenEval[indexpob]));
    }
    // coloca tamListaConexPost[indexpob] en el valor obtenido del ordenamiento.
    conf->tamListaConexPost[indexpob]=ordenarListaConexF(conf->headerSNN[indexpob], conf->listaConexData[indexpob] ,conf->ordenEval[indexpob],conf->tamOrdenEval[indexpob]);
    //libera valcalculado
    free(valCalculado);
    return(1);
}

int ordenarListaConexF(hdrSNNv1 headerSNN, tConexDataF* listaConexData, int* ordenEval, int tamOrdenEval)
// usando el arreglo ordenEval, organiza la lista de conexiones para que queden primero las
// que están de últimas en el arreglo ordenEval
{
    int i,j;
    int tmp=0;
    tConexDataF tempConex;
    // para cada valor de ordenEval leyendolo desde el final.
    for (i=(tamOrdenEval-1); i>=0; i--)
    {
        // para todas las conexiones j
        for (j=0; j<headerSNN.numConex; j++)
        {
            // si conexOut = ordenEval[i]
            if (listaConexData[j].conexOut==ordenEval[i])
                // si está enabled
                if (listaConexData[j].enabled)
                {
                    // hace swap
                    tempConex=listaConexData[j];
                    listaConexData[j]=listaConexData[tmp];
                    listaConexData[tmp]=tempConex;
                    // incrementa contador
                    tmp++;
                }
        }
    }
    return(tmp);
}

int recurOrganicer(int indexNodoOut, hdrSNNv1 headerSNN ,tConexDataF* listaConexData, int* ordenEval, char* valCalculado, int* cont)
// usa arreglo conf->valCalculado[] para marcar todos como no calculados al principio.
{
    int i;
    // si está marcado como calculado, retorna 1
    if (valCalculado[indexNodoOut]) return(1);
    // marca indexNodoOut como Calculado
    valCalculado[indexNodoOut]=1;
    // adiciona indexNodoOut a ordenEval[*cont]
    ordenEval[*cont]=indexNodoOut;
    // incrementa *cont
    *cont=*cont+1;
    // busca entre las conex las que tengan como conexOut a indexNodoOut
    for (i=0; i<headerSNN.numConex; i++)
    {
        // si tiene conexOut=indexNodoOut
        if (listaConexData[i].conexOut==indexNodoOut)
            // si es enabled y no es bias (SOLO PUEDE HABER UN NODO BIAS)
            if ((listaConexData[i].enabled)&&(indexNodoOut!=headerSNN.numEntradas))
            {
                // si el nodo conexIn  está calculado,
                if (!valCalculado[listaConexData[i].conexIn])
                    // llama a recurOrganicer
                    recurOrganicer(listaConexData[i].conexIn, headerSNN ,listaConexData, ordenEval, valCalculado, cont);
            }
    }
    //retorna 1
    return(1);
}

int genOrdenEvalF(int numGenomas, hdrSNNv1* headerSNN, tConexDataF** listaConexData, int** ordenEval, int* tamOrdenEval, TConfig* conf)
// genera ordenEval[], orden de evaluación de los nodos de cada genoma (inverso a recursivo empezando de salida)
{
    int i;
    //determina el número máximo de nodos
    for (i=0; i<numGenomas; i++)
    {
        genOrdenEvalF1g(i,conf);
    }
    // libera el vector temporal valCalculado
    return(1);
}

void imprimirSeleccion(TConfig* conf)
{
// Imprime los areglos: numGenomasPorEspecia, actNumGenomasPorEspecie, y la matriz listaordenFitness
// SOLO se puede usar dentro de seleccionCrossover debido a que al final se liberan los punteros.
    int i=0;
    int j=0;
    // imprime numGenomasPorEspecie
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>\nnumGenomasPorEspecie:    ");
    for (i=0; i<conf->numEspecies; i++)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>%u=%u,",i,conf->numGenomasPorEspecie[i]);
    }
    // imprime actNumGenomasPorEspecie
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>\nactNumGenomasPorEspecie: ");
    for (i=0; i<conf->numEspecies; i++)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>%u=%u,",i,conf->actNumGenomasPorEspecie[i]);
    }
    // imprime listaOrdenFitness
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>\nlistaOrdenFitness:");
    for (j=0; j<conf->numEspecies; j++)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nE%2i = ",j);
        for (i=0; i<conf->actNumGenomasPorEspecie[j]; i++)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>%u,",conf->listaOrdenFitness[j][i]);
        }
    }
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>\n");
}

float fError(float valor, float salida)
{
    // retorna -1 si hay error.
    // función de error para un genoma
    //TODO: falta implementar cuando rango de salidas es (-1,1)
    //TODO: PROBANDO Cambiado de versión 0.66
    float tmpA=2.435;
    float a=fabs(valor-salida)/2;
    //fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>Valor=%3.3f  salida=%3.3f\n",valor,salida);
    //calculo del error a partir de la diferencia de salidas (Diferente para cada problema)
    //a = a<=0.1? a: a<=0.2?sqrt(a):sqrt(sqrt(a));
    // A = -2.30258509299/((100.0-MinPorcentGanancia)/100);

    // a=1-exp(-10*a*a);

    // Error para forex con error de tres lineas -0.5,0,0.5:
    if ((valor>0.5)&&(salida<0.5)) tmpA*=7;
    if ((valor<0.5)&&(salida>0.5)) tmpA*=7;
    if ((valor>0)&&(salida<0)) tmpA*=7;
    if ((valor<0)&&(salida>0)) tmpA*=7;
    if ((valor>-0.5)&&(salida<-0.5)) tmpA*=7;
    if ((valor<-0.5)&&(salida>-0.5)) tmpA*=7;
    a=1-exp(-tmpA*a*a);



    //Comprueba que el error esté entre 0 y 1 YA que en selección Crossover se parte de esto, valores superiores causan errores.

    //return( a>=0 ? a<=2 ? a : -1 : -1 );
    return(a);
}

unsigned buscarMejorFitness(TConfig* conf)  // OPTIMIZADA
{
//Busca el mejor fitness entre todos los elementos de la conf->población
//Se debe llamar al final de laFuncion evaluarPob evaluarPob
//Resetea los representantes a valor -1 para indicar que no se han asignado, luego se recalculan
//Retorna el indexpob del genoma con el mejor fitness entre toda la conf->población -1 si hay error
//Parémetros:	sPob = tamaño de la conf->población.
    float temp=-100000;//fitness inicial =-100000
    unsigned i;
    unsigned j=0;
//TODO: probando qutando calculo de representantes y limitando la búsqueda a los representantes actuales.
    /*	for (i=0;i<conf->numEspecies;i++){	 //resetea todos los valores de los representantes a -1 para indicar que no se han asignado.
    		conf->representantes[i]=UINT_MAX;
    	}
    	for (i=0;i<conf->sizePob;i++){ //calcula los representantes
    		if (conf->representantes[conf->pob[i].especie]!=UINT_MAX){//si representantes de especie de i es diferente de -1 entonces: verifica si i es mejor que su representante.
    			if(conf->pob[conf->representantes[conf->pob[i].especie]].fitness<conf->pob[i].fitness){
    				conf->representantes[conf->pob[i].especie]=i;
    			}
    		}
    		else{ //si representantes de especie de i es -1 entonces: hace el representante de especie =i.
    			conf->representantes[conf->pob[i].especie]=i;
    		}
    	}
    */
    //busca el mejor entre todos y los backups de representantes.
    for (i=0; i<(conf->sizePob+conf->numEspecies); i++)
    {
        if(conf->pob[i].fitness>temp)
        {
            temp=conf->pob[i].fitness;
            j=i;
        }
    }
    return(j);
}

unsigned contarGPEsp(unsigned indEspecie, TConfig* conf)  // OPTIMIZADA
{
    // retorna el número de genomas en la población que pertenecen a la especie indEspecie
    unsigned i;
    unsigned temp=0;
    for (i=0; i<conf->sizePob; i++)
        if(conf->pob[i].especie==indEspecie)
            temp++;
    return temp;
}

unsigned verificarMejor(TConfig* conf)  // OPTIMIZADA, //TODO: se debe quitar cuando ya no se necesite
{
//Verificar mejor
//retorna 0 si error, se debe usar después de cada especiación de población
    float a;
    a=conf->pob[buscarMejorFitness(conf)].fitness;
    if (conf->antMejor>a)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError en Funcion verificarMejor()");
        return(0);
    }
    else
        conf->antMejor=a;
//TODO: probando quitar esta comprobación debido a que puede dar un falso positivo si muta el representante y todavía no se actualiza el guardado en especiacion
    /*	for (i=0;i<conf->sizePob;i++){
    		if (conf->pob[i].fitness>conf->pob[conf->representantes[conf->pob[i].especie]].fitness){
    			fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 243 en Funcion VerificarMejor() comparando el fitness de cada genoma con el del representante de su esp.");
    			return(0);
    		}

    	}
    */
    return(1);
}

void swap(unsigned  int* a, unsigned  int* b)  // OPTIMIZADA
{
//unsigned  intercambia los valores apuntados por dos punteros.
    unsigned temp=*a;
    *a=*b;
    *b=temp;
}

unsigned cargarGenoma(unsigned indexpob, char *filename, TConfig* conf)  // OPTIMIZADA
{
//Lee un genoma desde un archivo y sobreescribe con él un genoma de conf->pob , la memoria necesaria para el arreglo de nodos y conexiones se gestiona desde esta Funcion.
//Prerequisito:  debe haberse reservado memoria para el genoma en conf->pob[indexpob] (se hace con primeraGen())
//El formato de entrada es(sin separadores): Genoma, Genoma.nodo, genoma.conex las longitudes a escribir
//de cada estructura se basan en el tamaño de Genoma, GenNodoF, GenconexF y en los valores Genoma.totalNodos
//y Genoma.totalConexiones
//Parámetros:	indexpob = indice del arreglo de genomas conf->pob que se va a reemplazar por el leído
//				filename = path y nombre de archivo del que se leerá el genoma
//Retorna 0 si hay error, 1 ok
    FILE *fileIn;
    size_t leidos=0;
    GenNodoF* temp;
    //Abre el archivo para lectura
    if ((fileIn=fopen(filename,"rb"))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 49 en funcion guardarMejorFitness(%u,%s) llamando a fopen(%s,\"br\")\n",indexpob,filename,filename);
        return(0);
    }
    if (feof(fileIn)!=0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 243-1 en Funcion cargarGenoma() llamando a la Funcion feof() antes de leer archivo, archivo vacío.\n");
        fclose(fileIn);
        return(0);
    }
    leidos=fread(&(conf->pob[indexpob]),sizeof(Genoma),1,fileIn); //lee Genoma
    if (feof(fileIn)!=0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 243-2 en Funcion cargarGenoma() llamando a la Funcion feof() después de leer genoma, genoma incompleto, faltan nodos.\n");
        fclose(fileIn);
        return(0);
    }
    if ((temp=(GenNodoF *)malloc(sizeof(GenNodoF)*conf->pob[indexpob].totalNodos))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 50 en funcion guardarMejorFitness(%u,%s) llamando a malloc(%u)\n",indexpob,filename,(unsigned  int)(conf->pob[indexpob].totalNodos*sizeof(GenNodoF)));
        fclose(fileIn);
        return(0);
    }
    //if (conf->pob[indexpob].nodo!=NULL)
    free(conf->pob[indexpob].nodo);
    if (temp!=NULL)
        conf->pob[indexpob].nodo=temp;
    leidos=fread(conf->pob[indexpob].nodo,1,conf->pob[indexpob].totalNodos*sizeof(GenNodoF),fileIn); //lee nodos NO es necesario leer listas de punteros ya que se actualizan entes de evaluación.
    if (!feof(fileIn))
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 243-3 en Funcion cargarGenoma() llamando a la Funcion feof() después de leer genoma, genoma incompleto, faltan conexiones.\n");
        fclose(fileIn);
        return(0);
    }
    if (conf->pob[indexpob].conex!=NULL) free((void *)conf->pob[indexpob].conex);
    if ((conf->pob[indexpob].conex=(GenConexF *)malloc(sizeof(GenConexF)*conf->pob[indexpob].totalConexiones))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 51 en funcion guardarMejorFitness(%u,%s) llamando a malloc(%u)\n",indexpob,filename,(unsigned  int)(conf->pob[indexpob].totalConexiones*sizeof(GenConexF)));
        fclose(fileIn);
        return(0);
    }
    leidos=fread(conf->pob[indexpob].conex,1,sizeof(GenConexF)*conf->pob[indexpob].totalConexiones,fileIn); //lee conexiones
    if (leidos<conf->pob[indexpob].totalConexiones*sizeof(GenConexF))
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 243-4 en Funcion cargarGenoma() llamando a la Funcion feof() después de leer genoma, genoma incompleto, conexiones incompletas.\n");
        fclose(fileIn);
        return(0);
    }
    //Cierra el archivo
    if (fclose(fileIn)!=0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 50 en funcion cargarGenoma(%u,%s) llamando a fclose(fileout))\n",indexpob,filename);
        return(0);
    }
    return(1);
}

unsigned guardarRepresentantesNNP(char *filename, TConfig* conf)  // OPTIMIZADA
// Guarda los representantes en un archivo en formato NNPv1
// retorna 0 si hubo error, 1 si ok.
{
    hdrNNPv1 headerNNP;
    hdrSNNv1 header;
    int i,j,k,result;
    float tmpPesoF;
    double tmpPesoD;
    char tmpChar;
    FILE* fileNNP;
    unsigned escritos; // para verificar el número de elementos escritos en el archivo
    unsigned* reps; //para ordenar por orden descendente de fitness

    // si el número de especies es <1 retorna error;
    if (conf->numEspecies<1)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 51 en funcion guardarRepresentantesNPP(), numEspecies==0\n");
        return(0);
    }
    // abre el archivo para escritura
    if ((fileNNP=fopen(filename,"wb"))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 51.36 en funcion guardarRepresentantesNPP() llamando a fopen()\n");
        return(0);
    }
    // inicializa los campos del header de SNN
    headerNNP.fileID[0] = 'N';
    headerNNP.fileID[1] = 'N';
    headerNNP.fileID[2] = 'P';
    headerNNP.version = 1;
    headerNNP.numGenomas = conf->numEspecies;
    // guarda el encabezado en el archivo NNP
    escritos=fwrite(&headerNNP, sizeof(hdrNNPv1),1,fileNNP);
    if (escritos<1)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 51.37 en guardarRepresentantesNNP() llamando a fwrite()\n");
        fclose(fileNNP);
        return(0);
    }
    // reserva memoria para reps
    reps = (unsigned*) malloc(conf->numEspecies*sizeof(unsigned));
    if (reps==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 51.38 en guardarRepresentantesNNP() llamando a malloc()\n");
        return(0);
    }
    // llena reps con los index de los representantes
    for (k=0; k<conf->numEspecies; k++)
    {
        reps[k]=conf->representantes[k];
    }
    // ordena por fitness los representantes en reps
    for (k=1; k<conf->numEspecies; k++)
    {
        for (i=0; i<(conf->numEspecies-k); i++)
        {
            if (conf->pob[reps[i]].fitness<conf->pob[reps[i+1]].fitness)
            {
                //intercambia elementos i y i+1 del arreglo reps
                j=reps[i];
                reps[i]=reps[i+1];
                reps[i+1]=j;
            }
        }
    }

    // para cada representante
    for (k=0; k<conf->numEspecies; k++)
    {
        // genera encabezado SNNv1
        header.fileID[0] = 'S';
        header.fileID[1] = 'N';
        header.fileID[2] = 'N';
        header.version = 1;
        header.usarSigned = conf->tSigma>=1000? 1 : 0;
        header.tamRegistros = conf->useFloat==0? 8 : 4;
        header.numEntradas = conf->numEntradas;
        header.numSalidas = conf->numSalidas;
        header.numBias = conf->numBias;
        header.numHiddens = conf->pob[reps[k]].totalNodos-(conf->numEntradas+conf->numSalidas+conf->numBias);
        header.numConex = conf->pob[reps[k]].totalConexiones;
        header.sigmaFactor = (double)conf->A;
        header.actThreshold = (double) conf->Fthreshold;
        header.lastFitness = (double) conf->pob[reps[k]].fitness;
        // guarda encabezado SNN
        escritos=fwrite(&header, sizeof(hdrSNNv1),1,fileNNP);
        if (escritos<1)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 51.38 en guardarRepresentantesNNP() llamando a fwrite()\n");
            fclose(fileNNP);
            return(0);
        }
        // Escribe los arreglos en orden: int conexIn[numConex],conexOut[numConex],enabled[],double peso[numConex]
        // para conexIn
        for (i=0; i<conf->pob[reps[k]].totalConexiones; i++)
        {
            // escribe los 3 arrays con los datos de la conex
            escritos=fwrite(&(conf->pob[reps[k]].conex[i].indexIn), sizeof(unsigned),1,fileNNP);
        }
        // para conexOut
        for (i=0; i<conf->pob[reps[k]].totalConexiones; i++)
        {
            // escribe los 4 arrays con los datos de la conex
            escritos=fwrite(&(conf->pob[reps[k]].conex[i].indexOut), sizeof(unsigned),1,fileNNP);
        }
        // para enabled
        for (i=0; i<conf->pob[reps[k]].totalConexiones; i++)
        {
            tmpChar = conf->pob[reps[k]].conex[i].enabled;
            escritos = fwrite(&tmpChar, sizeof(char),1,fileNNP);
        }
        // tamaño de lista de orden de evaluación de conexiones
        escritos = fwrite(&(conf->tamListaConexPost[reps[k]]), sizeof(int),1,fileNNP);
        if (escritos==0)
        {
            printf("ERROR");
            exit(0);
        }
        // genera y escribe lista de evaluación de conexiones
        for (i=0; i<conf->tamListaConexPost[reps[k]]; i++)
        {
            for(j=0; j<conf->pob[reps[k]].totalConexiones; j++)
            {
                if ((conf->listaConexData[reps[k]][i].conexIn==conf->pob[reps[k]].conex[j].indexIn)&&(conf->listaConexData[reps[k]][i].conexOut==conf->pob[reps[k]].conex[j].indexOut))
                    result=j;
            }
            escritos = fwrite(&result, sizeof(int),1,fileNNP);
        }
        // escribe  listaConexData[indexpob][tamListaConexPost[indexpob]] que es la lista de orden de evaluación de conex.
        //escritos = fwrite(conf->listaConexData[conf->representantes[i]], conf->tamListaConexPost[conf->representantes[i]]*sizeof(int),1,fileNNP);
        // para peso
        for (i=0; i<conf->pob[reps[k]].totalConexiones; i++)
        {
            if (header.tamRegistros==4)//para float
            {
                //hace cast para sacar double a partir de float.
                tmpPesoF=(float)conf->pob[reps[k]].conex[i].peso;
                // escribe los 3 arrays con los datos de la conex
                escritos=fwrite(&tmpPesoF, sizeof(float),1,fileNNP);
            }
            else //para double
            {
                //hace cast para sacar double a partir de float.
                tmpPesoD=(double)conf->pob[reps[k]].conex[i].peso;
                // escribe los 3 arrays con los datos de la conex
                escritos=fwrite(&tmpPesoD, sizeof(double),1,fileNNP);
            }
        }
    }
    free(reps);
    //cierra el archivo NNP
    if (fclose(fileNNP)!=0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 50 en funcion guardarRepresentantesNNP() llamando a fclose()");
        return(0);
    }
    return(1);
}

long unsigned calcularMemoriaUsada(unsigned sPob, TConfig* conf)  //OPTIMIZADA - //TODO Adicionar consumo de arreglos de punteros a nodos, conexiones, enabled y peso en estructura GenNodoF
{
//Calcula la memoria usada por los arreglos del programa en bytes
//Parámetros:	sPob
//Retorna el número de bytes usados por el programa en memoria.
    unsigned  long int acum;
    unsigned i;
    acum=(unsigned  int)(sPob*sizeof(Genoma));//Calcula bytes usados por genomas
    for (i=0; i<sPob; i++) //calcula bytes usados por arreglos de nodos y conexiones
        acum=acum+(unsigned  int)((conf->pob[i].totalNodos*sizeof(GenNodoF))+(conf->pob[i].totalConexiones*sizeof(GenConexF)));
    acum=(unsigned  long int)(acum+(conf->numEspecies*sizeof(unsigned  int)));//lista de conf->representantes
    acum=(unsigned  long int)(acum+(conf->contInnovNodo*sizeof(TListaInnov)));//lista de innovaciones de nodos
    acum=(unsigned  long int)(acum+(conf->contInnovNodo*sizeof(TNodoOut)));//y sus nodos de salida
    acum=(unsigned  long int)(acum+(conf->contInnovNodo*sizeof(TListaInnov)));
    acum=(unsigned  long int)(acum+(conf->contInnovCon*sizeof(TNodoOut)));//lista de innovaciones de conexiones
    acum=(unsigned  long int)(acum+(conf->numEspecies*sizeof(unsigned  int)));//lista de conf->contGeneracSinMejora
    acum=(unsigned  long int)(acum+(conf->numConservacion*sizeof(unsigned  int)));//lista de conservacion

    return(acum);
}

void calcularD(TConfig* conf)  //TODO, falta calcular params paa aproximaciones
{
    // para el tsigma escogido, calcula los parámetros (coeficientes ) de la función fSigma seleccionada.
    if (conf->tSigma==0)  //con y(-1)=0, y(0)=1 y y(0.5)=0.5
    {
        conf->A = 10;
        conf->D = 1.04;
        conf->F = -0.2;
    }
    if (conf->tSigma==4)  //con y(-1)=0, y(0)=1 y y(0.5)=0.5
    {
        conf->A = 2.435;
        conf->D = 1.09626166044;
        conf->F = -0.09626166044;
    }
    if (conf->tSigma==1003)  //con y(-1)=0, y(0)=1 y y(0.5)=0.5
    {
        conf->A = 2.435;
        conf->D = 1.09626166044;
        conf->F = -0.09626166044;
    }

}

float fSigma(float fX, unsigned param, float fD, TConfig* conf)  //NO OPTIMA//TODO: usar INLINE si es posible
{
//fSigma, retorna un float corespondiente a la funcion de activación seleccionada con param para una entrada X
//TODO: verificar rangos de entrada y salida de cada aprox sigma
    float y=0;
    //Param =	0 = sigma(0,1),		y = 1 / (1 + exp (- D * x)) -> corregido
    //		1 = sigma aprox (0,1)	y = 0.5 + x * (1 - abs(x) / 2), y=0 si x<=-1, y=1 si x>=1
    //		2 = elliot (0,1)	y = (x / 2) / (1 + |x|) + 0.5
    //		3 = binario (0,1)	y = x>=0 ? 1:0
    //		4 = gauss(0,1),		y = exp(- x * x)
    //		1000 = tanh(-1,1),		y = 2 / (1 + exp(-2 * x)) - 1
    //		1000 = elliot,(-1,1)	y = x / (1 + |x|)
    //		1002 = binario (-1,1) 	y = x=0 ? 1: -1;
    //		1003 = gaussAn (-1,1)	y= 2*exp(- x*x))-1
    if (param==0)
        return ((conf->D / (1.0 + exp (- conf->A * (fX+0.5))))+conf->F);
    if (param==1)
    {
        if (fX<=-1)
        {
            return (0);
        }
        if (fX>=1.0)
        {
            return (1.0);
        }
        if ((fX!=-1.0)&&(fX<1.0))
            y = 0.5 + fX * (1.0 - (fabs(fX) / 2.0));
        return (y);
    }
    if (param==2)
        return ((fX / 2.0) / (1.0 + fabs(fX)) + 0.5);
    if (param==3)
    {
        return (fX>=0 ? fX>=1 ? 0: cos(3.14*fX): fX<=-1 ? 0: fX+1) ;
        //return (fX>=0 ? fX>=1 ? 0: -fX+1: fX<=-1 ? 0: fX+1) ;  TODO: probando función de activación sinusoidal, podría usarse con valores fasoriales de entradas
    }
    if (param==4)
        //return (conf->D*exp(-conf->A*(fX * fX))+conf->F);
        return(exp(-conf->A*(fX * fX)));
    //return (exp(-3*(fX * fX)));
    if (param==1000)
        return (2.0 / (1.0 + exp(fX*-2.0)) - 1.0);

    if (param==1001)
        return (fX / (1.0 + fabs(fX)));
    if (param==1002)
    {
        return (fX>=0 ? 1.0:-1.0);
    }
    if (param==1003)
        return 2*(exp(-conf->A*(fX * fX)))-1;
    return(0);
}

void imprimirGenoma(unsigned index, TConfig* conf)  //OPTIMIZADA - //TODO: generar bmp del genoma
{
//imprime la principal información del genoma incluyendo nodos y conexiones
    unsigned i;
    printf("\n ordenEval\n");
    /*   for (i=0;i<conf->contO[index];i++)
       {
           printf("\nOrdenEval[%2d] =  %2d",i,conf->ordenEval[index][i]);
       }*/
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>\nNúmero de nodos = %u, MaxInnovNumNodo = %u\nNúmero de conexiones = %u, MaxInnovNumConex = %u\nListado de nodos Nindex=function,valor,threshold\n",conf->pob[index].totalNodos,conf->pob[index].maxInnovNumNodo,conf->pob[index].totalConexiones,conf->pob[index].maxInnovNumConex);
    for (i=0; i<conf->pob[index].totalNodos; i++)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>N%u=%u,%1.1f,%1.1f,",i,conf->pob[index].nodo[i].nodeFunction,conf->pob[index].nodo[i].valor,conf->pob[index].nodo[i].thNodo);
    }
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>\nListado de conexiones Cindex = IndexIn,IndexOut,peso,enabled\n");
    for (i=0; i<conf->pob[index].totalConexiones; i++)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>C%u=%u,%u,%1.1f,%u,",i,conf->pob[index].conex[i].indexIn,conf->pob[index].conex[i].indexOut,conf->pob[index].conex[i].peso,conf->pob[index].conex[i].enabled);
    }
    /*fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nEstructura de punteros. Nindex=indexHijo0,indexHijo1....\n");
    for (i=0; i<conf->pob[index].totalNodos; i++)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br> N%u(%u)=",i,conf->pob[index].nodo[i].contHijos);
        for (j=0; j<conf->pob[index].nodo[i].contHijos; j++)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>%u,",conf->pob[index].nodo[i].conexHijo[j]->indexIn);
        }
    }*/
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>\n");
}

unsigned verificarGenoma(unsigned index, TConfig* conf)  // OPTIMIZADA - //TODO: quitarla cuando no haya errores
{
// Verificar Genoma, verifica que los máximos innov num y máximos numNodo y conex correspondan con los que se encuentran en el genoma
// retorna 0 si hay error, 1 si OK
    unsigned i;
    // verifica maxInnovnumNodo
    for (i=0; i<conf->pob[index].totalNodos; i++)
    {
        if (i>0)
            if (conf->pob[index].nodo[i].innovNum==0)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError en pob[%u] verificando innovNum de gen nodo[%u]!=0  es igual a 0 el maxInnovNumNodo",index,i);
                return(0);
            }
        if (conf->pob[index].nodo[i].innovNum>conf->pob[index].maxInnovNumNodo)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError en pob[%u] verificando no corresponde el maxInnovNumNodo[%u]",index,i);
            return(0);
        }
    }
    // verifica maxInnovNumConex
    for (i=0; i<conf->pob[index].totalConexiones; i++)
    {
        /*	if (i>0)
        		if (conf->pob[index].nodo[i].innovNum==0){
        			fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError en nodo %u verificando innovNum de gen nodo!=0  es igual a 0 el maxInnovNumNodo\n",index);
        			return(0);
        		} */
        if (conf->pob[index].conex[i].innovNum>conf->pob[index].maxInnovNumConex)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError en pob[%u] verificando no corresponde el maxInnovNumConex[%u]",index,i);
            return(0);
        }
    }

    return(1);
}

unsigned verificarPob(TConfig* conf)  //OPTIMIZADA - //TODO: quitarla cuando no haya errores
{
//VerificarPob: verifica los innovnum de todos los genomas de la pob
//retorna 0 si hay error, 1 si OK
    unsigned i;
    for (i=0; i<conf->sizePob; i++)
    {
        if (verificarGenoma(i,conf)==0)
            return(0);
    }
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>Verificado ");
    return(1);
}

unsigned buscarInnovNodo(unsigned indexpob, unsigned innovNum, TConfig* conf)  // OPTIMIZADA
{
//Retorna el index del arreglo de nodos en un genoma en la conf->población, retorna -1 si no lo encuentra.
//Parámetros: 	indexpob 	= índice del genoma en el arreglo de conf->población.
//				innovNum 	= número de innovación buscado.
    unsigned i;
    //para cada i entre 0 y maxnodos busca el que tenga nodo.innovnum == al buscado
    unsigned totalNodos = conf->pob[indexpob].totalNodos;
    for (i=0; i<totalNodos; i++)
        if(conf->pob[indexpob].nodo[i].innovNum==innovNum)
            return(i);
    return(UINT_MAX);
}

unsigned buscarInnovConex(unsigned indexpob, unsigned innovNum, TConfig* conf)  //OPTIMIZADA
{
//Retorna el index del arreglo de conexiones en un genoma en la conf->población, retorna -1 si no lo encuentra.
//Parámetros: 	indexpob 	= índice del genoma en el arreglo de conf->población.
//				innovNum 	= número de innovación buscado.
    unsigned i;
    unsigned totalConex = conf->pob[indexpob].totalConexiones; //para acelerar evaluación en for.
    //para cada i entre 0 y maxnodos busca el que tenga nodo.innovnum == al buscado
    for (i=0; i<totalConex; i++)
        if(conf->pob[indexpob].conex[i].innovNum==innovNum)
            return(i);
    return(UINT_MAX);
}

unsigned buscarInnovConexPorNodos(unsigned indexpob, unsigned innovIn, unsigned innovOut, TConfig* conf)  //OPTIMIZADA
{
//Retorna el index del arreglo de conexiones en un genoma en la conf->población, retorna -1 si no lo encuentra.
//Parámetros: 	indexpob 	= índice del genoma en el arreglo de conf->población.
//				innovIn 	= número de innovación de nodo de Entrada Buscado.
//				innovOut 	= número de innovación de nodo de Salida Buscado.
    unsigned i;
    unsigned totalConex=conf->pob[indexpob].totalConexiones; //para acelerar evaluación de for.
    //para cada i entre 0 y maxnodos busca el que tenga nodo.innovnum == al buscado
    for (i=0; i<totalConex; i++)
        if (conf->pob[indexpob].conex[i].nodoIn==innovIn)
            if (conf->pob[indexpob].conex[i].nodoOut==innovOut)
                return(i);
    return(UINT_MAX);
}
// sudor y lágrimas.
