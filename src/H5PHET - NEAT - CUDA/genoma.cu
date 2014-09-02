/** Redes neuronales para NEAT - H file
	usan/modifican genomas únicamente (pueden usarse operaciones de genes de gen.h en un genoma).
*/

#ifndef PARAMS_H_INCLUDED
#include "params.h"
#define PARAMS_H_INCLUDED
#endif
#ifndef AUXILIARES_H_INCLUDED
#include "auxiliares.h"
#define AUXILIARES_H_INCLUDED
#endif
#ifndef GEN_H_INCLUDED
#include "gen.h"
#define GEN_H_INCLUDED
#endif
#ifndef ESPECIE_H_INCLUDED
#include "especie.h"
#define ESPECIE_H_INCLUDED
#endif
#include "genoma.h"

float evalGenom(int indexpob,TConfig* conf)
{


    int i=indexpob;
    int j,k,totalNodos,indSalidas,indBias,indInS;
    // float matrixF[200][200]; // mil nodos máximo
    float valor[500]; // valor de salida de la neurona
    // float valorAnt[2000]; // valor de salida de la neurona
    float acum;
    // calcula el total de Nodos.
    totalNodos=conf->headerSNN[i].numEntradas+conf->headerSNN[i].numSalidas+conf->headerSNN[i].numBias+conf->headerSNN[i].numHiddens;
    // inicializa el fitness del genoma  en 0
    conf->fitness[i]=0;
    // calcula límites de for para optimización.
    indSalidas=(conf->headerGTD.numEntradas+conf->headerSNN[i].numBias+conf->headerSNN[i].numSalidas);
    indBias=(conf->headerGTD.numEntradas+conf->headerSNN[i].numBias);
    indInS=conf->headerGTD.numEntradas*sizeof(float);
    // reordena la lista de conexiones para que quede simialr a evaluación recursiva.
    // para los que estén al mismo nivel, debe ir primero el que esté más cerca a las entradas.
    // evalúa los datos de GTD en las conex y produce el vector de fitness
    // coloca en 0 valor[] y valorAnt[] superiores a las entradas y bias
    for (k=indBias;k<totalNodos;k++)
    {
        valor[k]=0;
    }
    //coloca bias en 1 ://TODO: falta para varias bias, es necesario? NO
    valor[conf->headerGTD.numEntradas]=1;
    for (j=0;j<conf->numDatos;j++)
    {
       // coloca los valores de entrada de valor[] en valorAnt[] y copia los nuevos a valor[]
// INICIO KERNEL CUDA_1, CAMBIAR  memcpy
// HACER SYNC ANTES DE HACER EL MEMCPY DESDE LOCAL A SHARED de dataGTDf[][] y los acum de Vc y Vt para posterior cálculo de error
       // memcpy(valorAnt,valor,indInS);
        memcpy(&(valor[0]),conf->dataGTDf[j],indInS);
        // coloca bias con valor 1 .
        valor[conf->headerGTD.numEntradas]=1;
        // calcula el fSigma de las entradas
        for (k=0;k<conf->headerGTD.numEntradas;k++)
        {
            valor[k]=2*((float)exp(-(float)conf->A*((valor[k]-(float)conf->Fthreshold)*(valor[k]-(float)conf->Fthreshold))))-1;;
        }
        // calcula el nuevo valor[] con el genoma i y dataGTD
        acum=0;
        for (k=0;k<conf->tamListaConexPost[i];k++)
        //for (k=0;k<headerSNN[i].numConex;k++)
        {
            acum+=(conf->listaConexData[i][k].peso*valor[conf->listaConexData[i][k].conexIn]);
            // si el conexOut no es el último
            if (k<(conf->tamListaConexPost[i]-1))
            {
                // si el ConexOut siguiente es diferente al actual
                if (conf->listaConexData[i][k].conexOut!=conf->listaConexData[i][k+1].conexOut)
                {
                    // calcula el sigma y actualiza valorAnt
                    valor[conf->listaConexData[i][k].conexOut]=2*((float)exp(-(float)conf->A*((acum-(float)conf->Fthreshold)*(acum-(float)conf->Fthreshold))))-1;
                    //reinicializa el acumulador
                    acum=0;
                }
            }
        }
        // saca el Fsigma del último nodo.

        valor[conf->listaConexData[i][conf->tamListaConexPost[i]-1].conexOut]=2*((float)exp(-(float)conf->A*(acum-(float)conf->Fthreshold)*(acum-(float)conf->Fthreshold)))-1;
        // copia valor a valorAnt
//        memcpy(valorAnt,valor,totalNodos*sizeof(float));
        // guarda los valores calculados en el arreglo valoresC[]
        for (k=indBias;k<indSalidas;k++)
        {
           conf->valoresC[j]=valor[k];
        }
        //FALTA: también calcular acumuladores de Vt y Vc para sacar medias al final desde host y luego llamar otro kernel que
        // calcula los errores. Lo mejor sería un kernel que calcule los vectores en un ciclo hasta numDatos/buffarsize
        // incluyendo transferencia desde buffer grande a shared y viceversa, luego
//FIN KERNEL CUDA_1 c, cambiar memcpy
    }
    // calcula el fitness como el coeficiente de correlación de Pearson de los 2 vectores (1 si =es,-1 si inversos, 0 si diferentes)
//INICIO KERNEL CUDA_2 (o parte final de CUDA_1)
    conf->fitness[i]=correlac(conf->valoresC,conf);
//FIN KERNEL CUDA_2 para calculo de fitness.

    return(conf->fitness[i]);
}

unsigned guardarGenomaSNN(unsigned indexpobA, char *filename, TConfig* conf)  // OPTIMIZADA
{
//Escribe un genoma deseado de conf->pob en un archivo
//El formato de salida es (sin separadores): Genoma, Genoma.nodo, genoma.conex las longitudes a escribir
//de cada estructura se basan en el tamaño de Genoma, GenNodoF, GenconexF y en los valores Genoma.totalNodos
//y Genoma.totalConexiones
//Parámetros:	indexpob = indice del arreglo de genomas conf->pob que se va a guardar
//				filename = path y nombre de archivo en el que se guardará el genoma
//Retorna 0 si hay error, 1 si ok
    int result;
    FILE *fileOut;
    size_t escritos=0;
    // Formato SNN: [header] unsigned conexIn[numConex], unsigned conexOut[numConex], unsigned ,double peso[numConex]
    int i,j,indexpob=indexpobA;
    double tmpPesoD;
    float tmpPesoF;
    char tmpChar;
//    float tmpFit=0;
    //si es un rep y si fitness <que backup, llama a las copiasdel backup
    if (indexpob>conf->sizePob)
    {
        if (conf->pob[indexpob].fitness>conf->pob[conf->representantes[indexpob-conf->sizePob]].fitness)
            copiarGenoma(indexpob,conf->representantes[indexpob-conf->sizePob],conf);
        indexpob=conf->representantes[indexpob-conf->sizePob];
    }
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
    // actualiza la lista de conexiones en orden de evaluación
//    tmpFit=evalGenom(indexpob,conf);
    printf("FitnessMejorGuardado=%7.7f\n",conf->pob[indexpob].fitness);
    // abre el archivo de salida para escritura
    if ((fileOut=fopen(filename,"wb"))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 49 en funcion guardarMejorFitness(%u,%s) llamando a fopen(%s,\"wb\")\n",indexpob,filename,filename);
        return(0);
    }
    // Escribe el encabezado SNN
    escritos=fwrite(&(conf->headerSNN[indexpob]), sizeof(hdrSNNv1),1,fileOut);
    // Escribe los arreglos en orden: int conexIn[numConex],conexOut[numConex],double peso[numConex]
    // para conexIn
    for (i=0; i<conf->pob[indexpob].totalConexiones; i++)
    {
        // escribe los 3 arrays con los datos de la conex
        escritos=fwrite(&(conf->pob[indexpob].conex[i].indexIn), sizeof(unsigned),1,fileOut);
    }
    // para conexOut
    for (i=0; i<conf->pob[indexpob].totalConexiones; i++)
    {
        // escribe los 3 arrays con los datos de la conex
        escritos=fwrite(&(conf->pob[indexpob].conex[i].indexOut), sizeof(unsigned),1,fileOut);
    }
    // para enableds
    for (i=0; i<conf->pob[indexpob].totalConexiones; i++)
    {
        tmpChar = conf->pob[indexpob].conex[i].enabled;
        escritos = fwrite(&tmpChar, sizeof(char),1,fileOut);
    }
    // tamaño de lista de orden de evaluación de conexiones
    if (conf->tamListaConexPost[indexpob]==0)
    {
        printf("ERROIR");
        exit(0);
    }
    escritos = fwrite(&(conf->tamListaConexPost[indexpob]), sizeof(int),1,fileOut);
    // genera y escribe lista de evaluación de conexiones
    for (i=0;i<conf->tamListaConexPost[indexpob];i++)
    {
        result=-1;
        for(j=0;j<conf->pob[indexpob].totalConexiones;j++)
        {
            if ((conf->listaConexData[indexpob][i].conexIn==conf->pob[indexpob].conex[j].indexIn)&&(conf->listaConexData[indexpob][i].conexOut==conf->pob[indexpob].conex[j].indexOut))
                result=j;
        }
        if (result==-1)
        {
            printf("\nError 56 en guardarGenomaSNN(), conex no encontrada");
            exit(0);
        }
        escritos = fwrite(&result, sizeof(int),1,fileOut);
    }
    //escritos = fwrite(conf->listaConexData[indexpob], conf->tamListaConexPost[indexpob]*sizeof(int),1,fileOut);
    // para peso
    for (i=0; i<conf->pob[indexpob].totalConexiones; i++)
    {
        if (conf->headerSNN[indexpob].tamRegistros==4)//para float
        {
            //hace cast para sacar double a partir de float.
            tmpPesoF=(float)conf->pob[indexpob].conex[i].peso;
            // escribe los 3 arrays con los datos de la conex
            escritos=fwrite(&tmpPesoF, sizeof(float),1,fileOut);
        }
        else //para double
        {
            //hace cast para sacar double a partir de float.
            tmpPesoD=(double)conf->pob[indexpob].conex[i].peso;
            // escribe los 3 arrays con los datos de la conex
            escritos=fwrite(&tmpPesoD, sizeof(double),1,fileOut);
        }
    }
    //cierra el archivo SNN.
    if (fclose(fileOut)!=0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 50 en funcion guardarMejorFitness(%u,%s) llamando a fclose(fileout))\n",indexpob,filename);
        return(0);
    }
    // libera memoria de los arreglos.
    return(1);
}
////TODO: el fitness se modifica en cada Funcion que modifique el genoma pero se debe usar conf.actualizarEnCambios para controlar esto.
//pero la Funcion debe llamar a evaluarpob y luego actualiza representantes.


int snnDataLoader(const unsigned indexpob, const unsigned especie, FILE* fileIn, TConfig* conf)
// Carga un snn en indexpob desde un archivo fileIn que debe estar abierto reconstruyendo el genoma.
{
    hdrSNNv1 headerSNN;
    unsigned* conexIn;
    unsigned* conexOut;
    char* enabled;
    float* pesosF=NULL;
    double* pesosD=NULL;
    int i, j , tmp, leidos,tmpNumNodos;

    // lee encabezado SNNv1
    leidos=fread(&headerSNN,sizeof(hdrSNNv1),1,fileIn);
    if (leidos<1)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 344.5 en snnDataLoader() llamando a fread()");
        return(1);
    }
    //verifica fileId y version de headerSNN
    if((headerSNN.fileID[0]!='S')||(headerSNN.fileID[1]!='N')||(headerSNN.fileID[2]!='N')||(headerSNN.version!=1)||(headerSNN.numEntradas>32000)||(headerSNN.numBias>32000)||(headerSNN.numConex>4000000)||(headerSNN.numHiddens>32000)||(headerSNN.numSalidas>32000)||(headerSNN.sigmaFactor>10000)||(headerSNN.actThreshold>1000)||(headerSNN.tamRegistros>256)||(headerSNN.usarSigned>100))
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 345 en snnDataLoader() error en encabezado SNNv1 %c%c%c v=%d,ne=%d,%d,%d,%d,%d,%d",headerSNN.fileID[0],headerSNN.fileID[1],headerSNN.fileID[2],headerSNN.version,headerSNN.numEntradas,headerSNN.numBias,headerSNN.numHiddens,headerSNN.numSalidas,headerSNN.numConex,headerSNN.tamRegistros);
        return(0);
    }
    // reserva memoria para los arreglos  a leer o leer secuencialmente? ver formato SNN
    tmp=headerSNN.numConex;
    if (tmp<4000000)
    {
        conexIn = (unsigned *) malloc(tmp*sizeof(unsigned));
        if (conexIn==NULL)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 346 en snnDataLoader() llamando a malloc()");
            return(0);
        }
        conexOut = (unsigned *) malloc(tmp*sizeof(unsigned));
        if (conexOut==NULL)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 347 en snnDataLoader() llamando a malloc()");
            return(0);
        }
        enabled = (char *) malloc(tmp*sizeof(sizeof(char)));
        if (enabled==NULL)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 347.1 en snnDataLoader() llamando a malloc()");
            return(0);
        }
        if (headerSNN.tamRegistros==4)
        {
            pesosF = (float*) malloc(tmp*sizeof(float));
            if (pesosF==NULL)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 348 en snnDataLoader() llamando a malloc()");
                return(0);
            }
        }
        else if (headerSNN.tamRegistros==8)
        {
            pesosD = (double*) malloc(tmp*sizeof(double));
            if (pesosD==NULL)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 349 en snnDataLoader() llamando a malloc()");
                return(0);
            }
        }
        else
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 350 en snnDataLoader() temRegistros invalido.");
            return(0);
        }
    }
    //sino imprime error y retorna
    else
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 350.2 en snnDataLoader() numConex>4M.");
        return(1);
    }
    // lee los 4 arreglos.
    leidos=fread(conexIn,sizeof(unsigned),headerSNN.numConex,fileIn);
    leidos+=fread(conexOut,sizeof(unsigned),headerSNN.numConex,fileIn);
    leidos+=fread(enabled,sizeof(char),headerSNN.numConex,fileIn);
    leidos+=fread(&(conf->tamListaConexPost[indexpob]),sizeof(int),1,fileIn);
    leidos+=fread(conf->listaConexData[indexpob],sizeof(int),conf->tamListaConexPost[indexpob],fileIn);
    if (headerSNN.tamRegistros==4)
    {
        leidos+=fread(pesosF,sizeof(float),headerSNN.numConex,fileIn);
    }
    else if (headerSNN.tamRegistros==8)
    {
        leidos+=fread(pesosD,sizeof(double),headerSNN.numConex,fileIn);
    }
    else
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 350.3 en snnDataLoader() tamregistros desconocido");
        free(conexIn); free(conexOut); if (pesosD) free(pesosD); if(pesosF) free(pesosF);
        return(1);
    }
    // verifica si se leyeron correctamente.
    if (leidos!=((4*headerSNN.numConex)+1+conf->tamListaConexPost[indexpob]))
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 350.5 en snnDataLoader() llamando a fread()");
        free(conexIn); free(conexOut); if (pesosD) free(pesosD); if(pesosF) free(pesosF);
        return(1);
    }
    // hasta aquí se puede retornar sinmodificar el genoma por tanto, se retorna 1.
    // crea nuevo genoma con numEntradas, numSalidas, numBias,etc..
    if (genomaInicial(indexpob, headerSNN.numEntradas, headerSNN.numSalidas, headerSNN.numBias, 0, especie, conf)==0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 351 en snnDataLoader() llamando a genomaInicial()");
        return(0);
    }
    // calcula el total de nodos del SNN
    tmpNumNodos=headerSNN.numEntradas+headerSNN.numSalidas+headerSNN.numBias+headerSNN.numHiddens;
    // verifica que valor leido no sea demasiado grande o corrupto
    if (tmpNumNodos>640000)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 352 en snnDataLoader() número de nodos de un SNN demasiado grande>64k.");
        return(0);
    }
    // reconstruye el genoma indexpob a partir de los arreglos leidos
    // para cada conexion i
    for (i=0;i<headerSNN.numConex;i++)
    {
        // si conexIn==pob[indexpob].numNodos;
        if (conexIn[i]==conf->pob[indexpob].totalNodos)
        {
            // busca en las conex actuales la que tiene ConexIn[i+1],conexOut[i]
            j=0;
            while ((j<i)&&((conexIn[i+1]!=conf->pob[indexpob].conex[j].indexIn)||(conexOut[i]!=conf->pob[indexpob].conex[j].indexOut)))
                j++;
            // si no la encontró, imprime error
            if (j==i)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 352.1 en snnDataLoader() conexión no encontrada..");
                return(0);
            }
            // nuevoNodo(ConexIn[i+1],conexOut[i])
            if (nuevoNodo(indexpob,j,conf)==0)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 352.1 en snnDataLoader() llamando a nuevoNodo..");
                return(0);
            }
            // actualizar peso y enabled de conex i
            conf->pob[indexpob].conex[i].enabled=enabled[i];
            if (headerSNN.tamRegistros==4)
            {
                conf->pob[indexpob].conex[i].peso = pesosF[i];
            }
            else if (headerSNN.tamRegistros==8)
            {
                conf->pob[indexpob].conex[i].peso = (float)pesosD[i];
            }
        }
        // si el nodo de entrada ya existe,
        else
        {
            // busca en conex actuales la que tenga conexIn[i] y conexOut[i]
            j=0;
            while ((j<i)&&((conexIn[i]!=conf->pob[indexpob].conex[j].indexIn)||(conexOut[i]!=conf->pob[indexpob].conex[j].indexOut)))
                j++;
            //si existeConex(conexin a conexout)
            if (j!=i)
            {
                //actualiza peso y enabled
                conf->pob[indexpob].conex[i].enabled=enabled[i];
                if (headerSNN.tamRegistros==4)
                {
                    conf->pob[indexpob].conex[i].peso = pesosF[i];
                }
                else if (headerSNN.tamRegistros==8)
                {
                    conf->pob[indexpob].conex[i].peso = (float)pesosD[i];
                }
            }
            // si no existe
            else
            {
                // verifica que el nodoIn y el Out Existan
                if ((conexIn[i]<conf->pob[indexpob].totalNodos)&&(conexOut[i]<conf->pob[indexpob].totalNodos))
                {
                    // crea la nueva conex(conexin a conexout)
                    nuevaConex(indexpob,conexIn[i],conexOut[i],((headerSNN.tamRegistros==4)? pesosF[i]: (float)pesosD[i]),0,enabled[i],conf);
                }
            }
        }
    }
    conf->pob[indexpob].fitness=(float)headerSNN.lastFitness;
    free(conexIn);
    free(conexOut);
    free(enabled);
    if (pesosF)
        free(pesosF);
    if (pesosD)
        free (pesosD);
    return(1);
}

 float calcularValorNodo(Genoma* pGenoma, int indexPob,GenNodoF* nodo, int indexNodo, TConfig* conf)  //OPTIMIZADA TODO: debería ir en gen.c
{
    // calcula recursivamente valores de salida (después de pasar por fsigma) de un nodo de un genoma indicado.
    // retorna el valor del nodo calculado, también lo asigna a valor y coloca en 1 valorCalculado.
    // parámetros:  indexPob = index del genoma en la población.
    //              nodo = puntero al nodo para el cual se quiere calcular el valor.
    unsigned i;
    unsigned contHijos = nodo->contHijos; // para que no se tenga que que hacer operación de consulta cada vez durante el for durante la condición.
    GenNodoF* pNodoHijo; //usado para acelerar el calculo
    float acum = 0;
    // verifica si el valor del nodo ya ha sido calculado.
    /*	if (nodo->estadoC == 2){ // si el valor del nodo ya fué calculado, solo lo retorna.
            return(nodo->valor);
    	}
    */
    if (nodo->estadoC == 1)  // si el valor del nodo ha empezado ha calcularse, pero aún no se tiene resultado retona el valor anterior.
    {
        // TODO: Para multiprocesamiento, se debe verificar si el nodo hijo corresponde a alguno de los padres
        // para procesamiento de un solo hilo, se retorna el valor anterior si es solicitado por
        return(nodo->valor);
    }
    // marca el valor del nodo como calculado, debe ir aquí por si recursivamente es solicitado este nodo, se retorne su valor anterior
    // también por esta razón se debe llamar calcularValorNodo exclusivamente en las salidas para que se calcule hacia abajo.
    nodo->estadoC = 1;
    // llena el valor de conf->ordenEval[indexpob][conf->contO[indexpob]] y incrementa conf->contO[indexpob]
//printf("oE[%d][%d]=%d,",indexPob,conf->contO[indexPob],indexNodo);
    conf->ordenEval[indexPob][conf->contO[indexPob]]=indexNodo;
    conf->contO[indexPob]++;
    // :)
    // FALTA: implementar ordenador de ordenFitness (sirve el que ya hice?)
    // FALTA: modificar guardarGenoma para agregar el arreglo ordenEval[indexpob][contO[indexpob]]
    // FALTA: // verificar si no es necesario para experimento cargarSNN.
    // FALTA: modificar SNNeval
    // para todos los hijos del nodo.
//PROBANDO, la sig linea no iba
    if (nodo->nodeFunction==0)//nodo->valor;
    {
//        if (contHijos==0)
 //           return(nodo->valor);
  //      else
            acum=nodo->valor;
    }
    for (i=0; i<contHijos; i++)
    {
        // si la conex con el hijo está enabled
        if (nodo->conexHijo[i]->enabled==1)  //TODO: Cuando esté funcionando, no colocar conexiones disableds en actualizarPNodos y quitar esta comprobación
        {
            // si el valor del nodo ya ha sido calculado
            pNodoHijo = &(pGenoma->nodo[nodo->conexHijo[i]->indexIn]);
            if (pNodoHijo->estadoC==1)
            {
                acum += nodo->conexHijo[i]->peso * pNodoHijo->valor;
            }
            else  //Si el valor del nodo hijo no ha sido calculado ERROR, falta INDEXNODO de nodo HIJO
            {
                acum += nodo->conexHijo[i]->peso*calcularValorNodo(pGenoma,indexPob,pNodoHijo,nodo->conexHijo[i]->indexIn,conf);
            }
        }
    }
    // retorna el fsigma.
    return(nodo->valor = fSigma(acum-nodo->thNodo,conf->tSigma,conf->fSigmaD,conf));
}

unsigned actualizarPNodos(unsigned index, TConfig* conf) //OPTIMIZADA
{
    //actualiza los valores de los punteros a conexiones hijo para cada nodo de un genoma.
    //se debe llamar al inicio de evaluarGenoma.
    //retorna 0 si hay error, 1 si ok
    unsigned i;
    unsigned j,k;
    Genoma* pGenoma = &(conf->pob[index]);
    GenNodoF* pNodo; //usado para acelerar los calculos
    //para cada gen nodo del genoma se reserva memoria para arreglo de conexHijo.
    for (i=0; i<pGenoma->totalNodos; i++)
    {
        pNodo = &(pGenoma->nodo[i]);
        if ( pNodo->nodeFunction != 3)
        {
            if ((pNodo->conexHijo = (GenConexF**)malloc(sizeof(GenConexF*)*pNodo->contHijos))==NULL)
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 1101 en Funcion actualizarPNodos() llamando a malloc() o contHijos=0");
                return(0);
            }
        }
        else
        {
            pNodo->contHijos=0;
            pNodo->conexHijo=NULL;
        }
        //para cada gen conex del genoma: busca las que tienen como indexOut a i y adiciona un puntero a esta conexión como hijo del nodo
        k=0;
        for (j=0; j<pGenoma->totalConexiones; j++)
        {
            if(pGenoma->conex[j].indexOut==i)
            {
                //adiciona al arreglo de conexHijo el ptr a la conex j
                /*				if (pGenoma->nodo[i].contHijos<=k){
                                    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 503.5 en actualizarPNodos() :pGenoma->nodo[i].contHijos<=k ");
                                    return(0);
                				}
                */
                pGenoma->nodo[i].conexHijo[k] = &(pGenoma->conex[j]);
                k++;
            }
        }
    }

    /*    //para cada gen conex del genoma:
        for (i=0;i<pGenoma->totalConexiones;i++){
            //adiciona al arreglo de conexHijo el ptr a la conex i
            indexNodoPadre=pGenoma->conex[i].indexOut;
            //para cada hijo
            for (j=0;j<conf->pob[index].nodo[indexNodoPadre].contHijos;j++){
                pGenoma->nodo[indexNodoPadre].conexHijo[j] = &(pGenoma->conex[i]);
            }
        }
    */
    return(1);
}

unsigned copiarGenoma(unsigned srcindexpob,unsigned dstindexpob,TConfig* conf)  //OPTIMIZADA
{
//Copia el genoma de origen al genoma de destino (el de destino es borrado)
//SUPONE que ya se ha reservado memoria para el array conf->pob que incluye a los dos elementos: funcion inicializarPob();
//Retorna 0 si hay error, 1 si Ok
    Genoma* pGenomaDst=&(conf->pob[dstindexpob]);
    Genoma* pGenomaSrc=&(conf->pob[srcindexpob]);
    // si fuente y destino son iguales, solo retorna 1
    if (srcindexpob==dstindexpob)
        return(1);
    //copia el contenido de conf->pob[src] al contenido de  conf->pob[dest] (copia un Genoma)
    if (pGenomaDst->nodo!=NULL) free((void *)pGenomaDst->nodo);
    if (pGenomaDst->conex!=NULL) free((void *)pGenomaDst->conex);
    //Copia los valores de la estructura genoma en el indexPob de destino.
    *pGenomaDst=*pGenomaSrc;
    //reserva memoria en el puntero conf->pob[dest].nodo con tamaño conf->pob[dest].totalNodos
    if ((pGenomaDst->nodo=(GenNodoF *) malloc(pGenomaSrc->totalNodos * ((unsigned  int)sizeof(GenNodoF))))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 31 en funcion copiarGenoma(%u,%u) llamando a malloc(%u).\n",srcindexpob,dstindexpob,pGenomaSrc->totalNodos * ((unsigned  int)sizeof(GenNodoF)));
        return(0);
    }
    //y copia el contenido del puntero conf->pob[src].nodo al contenido de conf->pob[dest].nodo
    pGenomaDst->nodo=(GenNodoF*)memcpy(pGenomaDst->nodo,pGenomaSrc->nodo, pGenomaSrc->totalNodos * ((unsigned  int)sizeof(GenNodoF)));
    if ((pGenomaDst->conex=(GenConexF *) malloc(pGenomaSrc->totalConexiones * ((unsigned  int)sizeof(GenConexF))))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 32 en funcion copiarGenoma(%u,%u) llamando a malloc(%u).\n",srcindexpob,dstindexpob,pGenomaSrc->totalConexiones * ((unsigned  int)sizeof(GenConexF)));
        return(0);
    }
    //y copia el contenido del puntero conf->pob[src].conex al contenido de conf->pob[dest].conex
    pGenomaDst->conex=(GenConexF*)memcpy(pGenomaDst->conex,pGenomaSrc->conex, pGenomaSrc->totalConexiones * ((unsigned  int)sizeof(GenConexF)));
    return(1);

}

/* se puede dejar otra Funcion que calcule el valor solo con indice en la struct GenNodoF para ahorrar memoria(que puede ser limitada) pero sacrificando velocidad.
Estado anterior de unsigned evaluarGenoma(unsigned index ,unsigned primero, float *entradas, float *salidas, unsigned nEntradas, unsigned nSalidas, TConfig* conf){ //OPTIMIZANDO
// Funcion evaluarGenoma() para un genoma i obtiene los valores y los fitness NO ajustados(1-error) para cada NODO de un genoma (incluyendo las salidas)
// también acumula fitness en pob[index].fitness para ser procesado luego por evaluarPob
// Parámetros:
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
//contadores
	unsigned i=0;
	unsigned j=0;
	float acum=0;//para almacenar entradasoralmente la sumatoria de cada nodo.
	unsigned tmp=0;
	unsigned inicOcultas=(conf->numEntradas+conf->numSalidas+conf->numBias);//posición inicial de nodos ocultos
	//si es primera vez (primero=1) hace cero todos los demás conf->pob[index].GenNodo[i].Valor.
	if (primero==1){
		fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>Reseteando valores de salida de neuronas");
		for (i=nEntradas;i<conf->pob[index].totalNodos;i++){
			if (conf->pob[index].nodo[i].nodeFunction!=3)
				conf->pob[index].nodo[i].valor=0;
			}
	}
	//Coloca los elementos leidos en los valores de los nodos correspondientes de entrada
	for (i=0;i<nEntradas;i++){
//		conf->pob[index].nodo[i].valor=entradas[i]; QUITADO POR TESTING DE VALOR DE SALIDA DE NEURONAS DE ENTRADA
// //TODO: para cada parámetro (pesos y thresholds) se debe guardar un registro de si el último cambio aumentó el fitness y el signo
//del incremento, si en el último incremento hubo mejora, continuar con el mismo signo
//sinó, usar el signo contrario para el próximo incremento aleatorio.
		conf->pob[index].nodo[i].valor=fSigma(entradas[i]-conf->pob[index].nodo[i].thNodo,conf->tSigma,conf->fSigmaD);
	}

//TESTING: calcula loa valores de salida de las neuronas de entrada


//Calcula los valores  de salida de las neuronas ocultas.(Las entradas siempre están al principio)
	//Para todos los que tengan GenNodo(i).nodeFunction=1 (ocultas) busca todos los que tengan GenNodo(j).nodoOut=i y los suma en acum (multiplicados por el peso de la conexión).
	for (i=inicOcultas;i<conf->pob[index].totalNodos;i++){
		acum=0;
		if(conf->pob[index].nodo[i].nodeFunction==1){//para nodos ocultos////TODO ESTA CONDICION SE PUEDE QUITAR
			for (j=0;j<conf->pob[index].totalConexiones;j++){//para cada conex
				if(conf->pob[index].conex[j].enabled==1){ //Si la conexión está enabled
					if(conf->pob[index].conex[j].nodoOut==conf->pob[index].nodo[i].innovNum){//si la conex j tiene como nodo de salida el nodo i.innovnum
						if((tmp=buscarInnovNodo(index,conf->pob[index].conex[j].nodoIn,conf))==UINT_MAX){//busca el index del nodo de entrada de la conex j
							fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 8 en Funcion evaluarGenoma(%u,%u,%u,%u) llamando a buscarInnovNodo(%u,%u)\n",index,primero,nEntradas,nSalidas,index,conf->pob[index].conex[j].nodoIn);
							return(0);
						}
						acum+=(conf->pob[index].conex[j].peso*conf->pob[index].nodo[tmp].valor);
					}
				}
			}
//Coloca en el valor de cada nodo i sigma de acum,threshold.
			conf->pob[index].nodo[i].valor=fSigma(acum-conf->pob[index].nodo[i].thNodo,conf->tSigma,conf->fSigmaD);
		}
	}
//Calcula los valores de las neuronas de salida (Las entradas siempre están al principio).
	//Para todos los que tengan GenNodo(i).function=0 (entradas) busca todos los que tengan GenNodo(j).nodoOut=i y los suma en acum (multiplicados por el peso de la conexión).
	for (i=(conf->numEntradas+conf->numBias);i<inicOcultas;i++){
		acum=0;
		if(conf->pob[index].nodo[i].nodeFunction==2){//para nodos de salida//ESTA CONDICION SE PUEDE QUITAR //TODO
			for (j=0;j<conf->pob[index].totalConexiones;j++){
				if(conf->pob[index].conex[j].enabled==1){ //Si la conexión está enabled
					if(conf->pob[index].conex[j].nodoOut==conf->pob[index].nodo[i].innovNum){
						if((tmp=buscarInnovNodo(index,conf->pob[index].conex[j].nodoIn,conf))==UINT_MAX){
							fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 9 en Funcion evaluarGenoma(%u,%u) llamando a buscarInnovNodo(%u,%u)\n",index,primero,index,conf->pob[index].conex[j].nodoIn);
							imprimirGenoma(index,conf);
							return(0);
						}
						acum+=(conf->pob[index].conex[j].peso*conf->pob[index].nodo[tmp].valor);
					}
				}
			}
//Coloca en el valor de cada nodo i sigma de acum,threshold.
			conf->pob[index].nodo[i].valor=fSigma(acum-conf->pob[index].nodo[i].thNodo,conf->tSigma,conf->fSigmaD);
		}
	}

//Habiendo calculado las salidas, obtiene el fitness (NO AJUSTADO) de las salidas actuales respecto a las salidas
//el error el es el promedio de errores ABSOLUTOS de todas las salidas respecto a las salidas deseadas.
	acum=0;
	j=0;
	for (i=(conf->numEntradas+conf->numBias);i<(inicOcultas);i++){//para nodos de salida
		if(conf->pob[index].nodo[i].nodeFunction==2){//TODO para nodos de salida//ESTA CONDICION SE PUEDE QUITAR
			if (j==nSalidas){//TODO ESTE BLOQUE SE PUEDE QUITAR //TODO
				fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 10 en funcion evaluarGenoma(%u,%u) número de salidas %u no corresponde con neuronas de salida en el genoma.",index,primero,nSalidas);
				return(0);
			}
			acum=acum+fabs(conf->pob[index].nodo[i].valor-salidas[j]);////TODO , mejorar CALCULO (se quitó normalización) DE ESTE ERROR  SI salida[j]=0
			j++;
		}
	}

	conf->pob[index].fitness += (acum/(float)conf->numSalidas);

	//Actualiza la lista de representantes (se busca el MENOR fitness debido a que todavía no se ha calculado el fitness real)
	if(conf->pob[index].fitness<conf->representantes[conf->pob[index].especie]){
		conf->representantes[conf->pob[index].especie]=index;
	}

	//fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>F%1.1f\n",conf->pob[index].fitness);
//	if(conf->pob[index].fitness<0)
//	{
//	fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>Ins=%1.1f,%1.1f, Nout=%1.1f, d=%1.1f fit=%1.1f, acum=",conf->pob[index].nodo[0].valor,conf->pob[index].nodo[1].valor,conf->pob[index].nodo[3].valor,salidas[0],conf->pob[index].fitness);
//		return(0);
//	}
// :)

	return(1);
}*/

unsigned crossover(unsigned indexpob1, unsigned indexpob2, unsigned indexpobOut ,float super, float promediarPob,float porcentEnableds,TConfig* conf) //OPTIMIZADA TODO: debería ir en especie
{
//Realiza el cruce entre dos genomas dados sus indexpob y lo coloca en un elemento de conf->pob (puede ser uno de los padres)
//Toma calcula el fitness de los dos genomas, se heredan los matching genes randomly,
//Los disjounsigned y excess se heredan solo del fittest,
//Parámetros: 	indexpob1, indexpob2 = genomas a cruszar
//				indexpobOut = indexpob donde se debe colocar el genoma resultante del cruce.
//				super = float entre 0 y 1, Probabilidad de heredar los exess y disjounsigned  ints del menos apto (aparte de los que se heredan normalmente del más apto)
//				promendiarProb = float entre 0 y 1 = probabilidad de que en caso de matching, se promedien los pesos en lugar de
    //seleccionarlos aleatoriamente entre los padres.
//Retona 0 si hubo error, 1 si ok, coloca en la variable global tempGenoma el genoma generado por crossover
    unsigned mejor=indexpob1;
    unsigned peor=indexpob2;
    unsigned i;
    unsigned j;
    Genoma* pGenomaMejor; //usados para acelerar calculoa
    Genoma* pGenomaPeor;

    //TODO QUITAR: verifica que los máximos innov num y máximos numNodo y conex correspondan con los que se encuentran en los genoma
    if ((verificarGenoma(indexpob1,conf)==0)||(verificarGenoma(indexpob2,conf)==0)||(verificarGenoma(indexpobOut,conf)==0))
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 46.0 en funcion crossover(%u,%u,%u,%1.1f,%1.1f) llamando a verificarGenomas()\n",indexpob1, indexpob2, indexpobOut , super, promediarPob);
        return(0);
    }//TODO QUITAR CUANDO NO HAYA ERRORES

    // coloca en mejor el indexpob del que tiene mayor fitness
    if (conf->pob[indexpob2].fitness>conf->pob[indexpob1].fitness)
    {
        mejor=indexpob2;
        peor=indexpob1;
    }
    // si fitness mayor= fitness menor, si super=0 entonces mejor=el que tenga menor totalconexiones, peor el que tenga mayor totalConexiones
    if (conf->pob[indexpob2].fitness==conf->pob[indexpob1].fitness)
    {
        if (conf->pob[indexpob2].totalConexiones<conf->pob[indexpob1].totalConexiones)
        {
            mejor=indexpob2;
            peor=indexpob1;
        }
    }
    // asigna punteros al mejor y peor genoma
    pGenomaMejor = &(conf->pob[mejor]);
    pGenomaPeor = &(conf->pob[peor]);

    // probando
   // fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>%u>%u,", mejor,indexpobOut);
   // copiarGenoma(mejor, indexpobOut,conf);

// TODO:Probando nuevo crossover que solo copia El mejor a hijo si hijo es =! peor  y luego comparando conex a conex (y nodo por fdt?)
//      en caso de genes con ==innovNumCon, se debe escoger uno de los dos pesos al azar para esa conex en el Hijo.
// TODO: PROBANDO POR ERROR DE SC

    // si hijo es != peor
    if (indexpobOut!=peor)
    {
        // copia mejor a hijo
        copiarGenoma(mejor,indexpobOut,conf); //   :)
    }
    // busca en el peor conexiones con =innovnum a alguna de m
    for (i=0;i<conf->pob[peor].totalConexiones;i++) // barre econex del peor
    {
        for (j=0;j<conf->pob[mejor].totalConexiones;j++) // barre conex del mejor
        {
            if (conf->pob[mejor].conex[j].innovNum==conf->pob[peor].conex[i].innovNum) // si los innovnums son iguales
            {
                if (((float)randL(conf))<=0.5){ // escoge al azar entre los dos.
                    // si el hijo no es el peor, asigna el peso del peor, porque originalmente tenía el del mejor.
                    if (indexpobOut!=peor)
                    {
                        conf->pob[indexpobOut].conex[j].peso=conf->pob[peor].conex[i].peso;
                        //solo si estaba disabled en el mejor la copia act/desact, para evitar no-viables.
                        if (conf->pob[indexpobOut].conex[j].enabled==0)
                            conf->pob[indexpobOut].conex[j].enabled=conf->pob[peor].conex[i].enabled;
                    }
                    // si el hijo es el peor, asigna el peso del mejor, porque originalmente tenía el del peor.
                    else
                    {
                        conf->pob[indexpobOut].conex[i].peso=conf->pob[mejor].conex[j].peso;
                        //solo si estaba disabled en el hijo(peor) la copia act/desact, para evitar no-viables.
                        if (conf->pob[indexpobOut].conex[i].enabled==0)
                            conf->pob[indexpobOut].conex[i].enabled=conf->pob[mejor].conex[j].enabled;
                    }
                }
            }
        }
    }


/** Probando nuevo crossover que solo copia El mejor a hijo si hijo es =! peor  y luego comparando conex a conex (y nodo por fdt?)
    //en caso de genes con ==innovNumCon, se debe escoger uno de los dos pesos al azar para esa conex en el Hijo.

    // copia com malloc+memcpy el pGenomaMejor->nodo[i]  (uno por uno creo) a tempGenoma.nodo
    if((tempGenoma.nodo=malloc(sizeof(GenNodoF)*pGenomaMejor->totalNodos))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 40 en funcion crossover(%u,%u,%u,%1.1f,%1.1f) llamando a malloc(%u)\n", indexpob1, indexpob2, indexpobOut , super, promediarPob,pGenomaMejor->totalNodos*(unsigned int)sizeof(GenNodoF) );
        return(0);
    }
    if ((tempGenoma.nodo=memcpy(tempGenoma.nodo,pGenomaMejor->nodo,sizeof(GenNodoF)*pGenomaMejor->totalNodos))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 41 en funcion crossover(%u,%u,%u,%1.1f,%1.1f) llamando a memcpy(), el arreglo de nodos del mejor es NULL\n", indexpob1, indexpob2, indexpobOut , super, promediarPob );
        return(0);
    }
    tempGenoma.totalNodos=pGenomaMejor->totalNodos;
    tempGenoma.maxInnovNumNodo=pGenomaMejor->maxInnovNumNodo;
    tempGenoma.maxInnovNumConex=pGenomaMejor->maxInnovNumConex;
    //para todos los nodos hace el contHijos=0;
    for (i=0; i<tempGenoma.totalNodos; i++)
    {
        tempGenoma.nodo[i].contHijos=0;
        tempGenoma.nodo[i].conexHijo=NULL; //al copiar inicializa este puntero en null para el nuevo elemento(Se configura en actualizarPNodos).
    }
    // actualiza totalconexiones y reserva memoria para nodos.
    tempGenoma.totalConexiones=pGenomaMejor->totalConexiones;
    tempGenoma.conex=(GenConexF *) malloc(sizeof(GenConexF)*tempGenoma.totalConexiones);
    if(tempGenoma.conex==NULL)
    {
       fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 42.1 en funcion crossover(%u,%u,%u,%1.1f,%1.1f) llamando a malloc(%u)\n",indexpob1, indexpob2, indexpobOut , super, promediarPob,tempGenoma.totalConexiones*(unsigned int)sizeof(GenConexF));
       return(0);
    }
    // para cada i desde 0, mientras i<pGenomaMejor->totalConexiones
    for (i=0; i<pGenomaMejor->totalConexiones; i++)
    {
        pTempGenomaConex = &(tempGenoma.conex[i]);
        // busca en conf->pob[peor] , pGenomaMejor->conex[i].innovNum, si lo encuentra, escoge leatoriamente entre mejor y peor y lo copia a tempgenoma
        if((j=buscarInnovConex(peor,pGenomaMejor->conex[i].innovNum,conf))!=UINT_MAX)
        {
            //TODO: verificar si este parametro funciona por encima de cero en redes grandes sinó quitarlo.
            if (((float)randL(conf))<=promediarPob	) //los promedia
            {
                *pTempGenomaConex=pGenomaMejor->conex[i];
                pTempGenomaConex->peso=(pGenomaMejor->conex[i].peso+pGenomaPeor->conex[j].peso)/2;

            }
            else //escoge aleatoriamente entre el mejor y el peor
            {
                //TODO: Probar colocando como parámetro si el porcentaje de probabilidad de heredar de mejor o peor y si influye en redes grandes
                if (((float)randL(conf))<=0.5)  //TODO: PUEDE HABER ERROR AQUI SE PODRÍA USAR GENOMA EN INDEXPOB DESPUES DE
                                                    //REPS
                    *pTempGenomaConex=pGenomaMejor->conex[i];
                else
                    *pTempGenomaConex=pGenomaPeor->conex[j];
            }
            if(pGenomaMejor->conex[i].enabled!=pGenomaPeor->conex[j].enabled)
            {
                if(((float)randL(conf))<porcentEnableds)
                    pTempGenomaConex->enabled=pGenomaPeor->conex[j].enabled;
                else
                    pTempGenomaConex->enabled=pGenomaMejor->conex[i].enabled;
            }
        }
        // si no lo encuentra, copia la conexión del mas apto a tempGenoma.
        else
        {
            // copia la conexión del mejor a tempgenoma
            *pTempGenomaConex=pGenomaMejor->conex[i];
        }
        //incrementa el contHijo del nodo padre de la conexión
        tempGenoma.nodo[pTempGenomaConex->indexOut].contHijos++;
    }
    //TODO: Verificar si este parámetro funciona en redes grandes, sinó quitarlo
    if (((float)randL(conf))<=super)
    {
        //Copia a tempGenoma los genes de nodo que no existan en tempGenoma.nodo.
        k=tempGenoma.totalNodos;
        for (i=0; i<pGenomaPeor->totalNodos; i++)
        {
            genFound=0;
            for (j=0; j<k; j++)
            {
                if(tempGenoma.nodo[j].innovNum==pGenomaPeor->nodo[i].innovNum)
                    genFound=1;
            }
            if (genFound==0)  //si no se encontró en tempGenoma, se le adiciona.
            {
                if ((tempGenoma.nodo=(GenNodoF *)realloc(tempGenoma.nodo,sizeof(GenNodoF)*(tempGenoma.totalNodos+1)))==NULL)
                {
                    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 43 en funcion crossover(%u,%u,%u,%1.1f,%1.1f) llamando a realloc(tempGenoma.nodo,%u)\n",indexpob1, indexpob2, indexpobOut , super, promediarPob,(tempGenoma.totalNodos+1)*(unsigned int)sizeof(GenNodoF));
                    free(tempGenoma.conex);
                    return(0);
                }
                tempGenoma.nodo[tempGenoma.totalNodos++]=pGenomaPeor->nodo[i];
                if (pGenomaPeor->nodo[i].innovNum>tempGenoma.maxInnovNumNodo)
                    tempGenoma.maxInnovNumNodo=pGenomaPeor->nodo[i].innovNum;
            }
        }
        //Copia a tempGenoma los genes disjounsigned o excess de pGenomaPeor->conex
        k=tempGenoma.totalConexiones;
        for (i=0; i<pGenomaPeor->totalConexiones; i++)
        {
            genFound=0;
            for (j=0; j<k; j++)
            {
                if(tempGenoma.conex[j].innovNum==pGenomaPeor->conex[i].innovNum)
                    genFound=1;
            }
            if (genFound==0)  //si no se encontró en tempGenoma, se le adiciona.
            {
                if ((tempGenoma.conex=(GenConexF *)realloc(tempGenoma.conex,sizeof(GenConexF)*(tempGenoma.totalConexiones+1)))==NULL)
                {
                    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 44 en funcion crossover(%u,%u,%u,%1.1f,%1.1f) llamando a realloc(tempGenoma.nodo,%u)\n",indexpob1, indexpob2, indexpobOut , super, promediarPob,(tempGenoma.totalConexiones+1)*(unsigned int)sizeof(GenConexF));
                    return(0);
                }
                tempGenoma.conex[tempGenoma.totalConexiones++]=pGenomaPeor->conex[i];
                if (pGenomaPeor->conex[i].innovNum>tempGenoma.maxInnovNumConex)
                    tempGenoma.maxInnovNumConex=pGenomaPeor->nodo[i].innovNum;
            }
        }
    }
    // coloca la especie igual a la del mejor y fitness en 0
    tempGenoma.especie=pGenomaMejor->especie;
    tempGenoma.fitness=0;
    // libera la mamoria usada por el arreglo de nodos y conexiones del nodo hijo a reemplazar.
    if (pGenomaOut->nodo!=NULL) free((void *)pGenomaOut->nodo);
    if (pGenomaOut->conex!=NULL) free((void *)pGenomaOut->conex);
    // copia el tempGenoma al indexpob genoma de la pob.
    *pGenomaOut=tempGenoma;
    // copia tempGenoma a pGenomaOut->nodo
    if((pGenomaOut->nodo=(GenNodoF *)malloc(sizeof(GenNodoF)*tempGenoma.totalNodos))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 45 en funcion crossover(%u,%u,%u,%1.1f,%1.1f) llamando a malloc(%u)\n",indexpob1, indexpob2, indexpobOut , super, promediarPob,(tempGenoma.totalNodos)*(unsigned int)sizeof(GenNodoF));
        free(tempGenoma.conex);
        return(0);
    }
    pGenomaOut->nodo=memcpy(pGenomaOut->nodo,tempGenoma.nodo,sizeof(GenNodoF)*tempGenoma.totalNodos);
    // copia tempGenoma a pGenomaOut->conex
    if((pGenomaOut->conex=(GenConexF *)malloc(sizeof(GenConexF)*tempGenoma.totalConexiones))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 46 en funcion crossover(%u,%u,%u,%1.1f,%1.1f) llamando a malloc(%u)\n",indexpob1, indexpob2, indexpobOut , super, promediarPob,(tempGenoma.totalConexiones)*(unsigned int)sizeof(GenConexF));
        free(tempGenoma.conex);
        return(0);
    }
    pGenomaOut->conex=memcpy(pGenomaOut->conex,tempGenoma.conex,sizeof(GenConexF)*tempGenoma.totalConexiones);
    //TODO: QUITAR esta Funcion si no hay errores
    if (verificarGenoma(indexpob1,conf)==0)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 46.5 en funcion crossover() llamando a verificarGenomas()\n");
        free(tempGenoma.conex);
        return(0);
    }
    if (verificarGenoma(indexpob2,conf)==0)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 46.6 en funcion crossover() llamando a verificarGenomas()\n");
        free(tempGenoma.conex);
        return(0);
    }
    if (verificarGenoma(indexpobOut,conf)==0)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 46.7 en funcion crossover() llamando a verificarGenomas()\n");
        free(tempGenoma.conex);
        return(0);
    }

    //calcula el max y min Th para el genoma nuevo.
    calcularLimThOneGenome(indexpobOut,conf);
    //libera la memoria usada por el arreglo de nodos y conex de tempgenoma
    if (tempGenoma.nodo!=NULL) free((void *)tempGenoma.nodo);
    if (tempGenoma.conex!=NULL) free((void *)tempGenoma.conex);
*/
    return(1);

}

unsigned evaluarGenoma(int indexPob, Genoma* pGenoma, unsigned primero, float *entradas, float *salidas, TConfig* conf)  //OPTIMIZADA
{
// Funcion evaluarGenoma() para un genoma i obtiene los valores y los fitness NO ajustados(1-error) para cada NODO de un genoma (incluyendo las salidas)
// también acumula fitness en pob[index].fitness para ser procesado luego por evaluarPob
// Parámetros:
//				index = indice de genoma por parámetro y
//				primero = especifica si es la primera vez (=1) que se evalúa el genoma para inicializar si no es la primera vez=0
//				entradas = puntero a un arreglo de nEntradas valores float que serán las entradas a evaluar
//				salidas = puntero a un arreglo de nSalidas valores float que serán las salidas deseadas, respecto a las cuales se obtendrá el error y por tant el fitness = 1-error.
//				nEntradas = número de elementos en el arreglo entradas
//				nSalidas = número de elementos en el arreglo de salidas
// sus valores en 0 excepto los de las entradas.
// salida: retorna el genoma en el indice index evaluado para la entrada con la variable valor de la estructura GenNodoF evaluada.
// retorna 0 si hubo error.
// contadores
    unsigned i=0;
    unsigned j=0;
    float acum=0;
    float tmp;
    GenNodoF* pNodo;
    unsigned numEntradas=conf->numEntradas;
    unsigned numEntradasBias = numEntradas + conf->numBias;
    unsigned inicOcultas=(numEntradasBias+conf->numSalidas);//posición inicial de nodos ocultos
    // si es primera vez (primero=1) hace cero todos los demás pGenoma->GenNodo[i].Valor y marca valornocalculado.
//PROBANDO PARA VER SI SE PUEDE LLEGAR A 90% CON ENTRADAS EN FSIGMA
    if (primero)
    {
        for (i=0; i<pGenoma->totalNodos; i++)
        {
            if ((pNodo = &(pGenoma->nodo[i]))->nodeFunction!=3)
            {
                pNodo->valor=0;
                pNodo->estadoC=0;
            }
        }
    }
    else  // si no es primero, solo marca como no actualizado el valor de cada nodo del genoma que no sea entrada o bias.
    {
//PROBANDO PARA VER SI SE PUEDE LLEGAR A 90% CON ENTRADAS EN FSIGMA
        //for (i=numEntradasBias; i<pGenoma->totalNodos; i++)
        for (i=0; i<pGenoma->totalNodos; i++)
        {
            if ((pNodo = &(pGenoma->nodo[i]))->nodeFunction!=3)
            {
                pNodo->estadoC=0;
            }
        }
    }
    // inicializa los contO[indexPob]=0;
    conf->contO[indexPob]=0;
    // coloca los elementos leidos en los valores de los nodos correspondientes de entrada
    for (i=0; i<numEntradas; i++)
    {
        //MOdificado en version0.66 para probar entrada lineal
        pGenoma->nodo[i].valor=entradas[i]; //FUNCIONA MEJOR EN REDES GRANDES QUE CON FSIGMA
        // TODO: para cada parámetro (pesos y thresholds) se debe guardar un registro de si el último cambio aumentó el fitness y el signo
        //del incremento, si en el último incremento hubo mejora, continuar con el mismo signo
        //sinó, usar el signo contrario para el próximo incremento aleatorio.
        //pNodo=&(pGenoma->nodo[i]);
        //pNodo->valor=fSigma(entradas[i]-pNodo->thNodo,conf->tSigma,conf->fSigmaD,conf);
    }
    // calcula los valores de las neuronas de salida (Las entradas siempre están al principio).
    // Para todos los que tengan GenNodo(i).function=0 (entradas) busca todos los que tengan GenNodo(j).nodoOut=i y los suma en acum (multiplicados por el peso de la conexión).
    for (i=numEntradasBias; i<inicOcultas; i++)
    {
        pGenoma->nodo[i].valor=calcularValorNodo(pGenoma,indexPob,&(pGenoma->nodo[i]),i,conf);
    }
    // habiendo calculado las salidas, obtiene el fitness (NO AJUSTADO) de las salidas actuales respecto a las salidas
    // el error el es el promedio de errores ABSOLUTOS de todas las salidas respecto a las salidas deseadas.
    acum=0;
    j=0;
    //TODO: este for y el anterior se pueden unir pero quitando acum del anterior for.

    for (i=numEntradasBias; i<inicOcultas; i++) //para nodos de salida
    {
        //TODO: hacer función de error(para que acum=acum+fError) de be tener como params el nodo[i].valor y la salida[j], el valor de retorno debe sor positivo entre 0 y 1
        tmp=fError(pGenoma->nodo[i].valor,salidas[j]);
		//fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>OutR=%3.3f, OutC=%3.3f, Error=%3.3f\n", salidas[j],pGenoma->nodo[i].valor, tmp );
        //if (index==conf->representantes[0]) fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>tmp=%3.3f, acum=%3.3f \n",pGenoma->nodo[i].valor, salidas[j]);
        if (tmp==-1)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>Error 345 en evaluarGenoma() llamando a fError, error negativo");
            return(0);
        }
        acum=acum+tmp;////TODO , mejorar CALCULO (se quitó normalización) DE ESTE ERROR  SI salida[j]=0
        j++;
    }
    // incrementa el fitness con el acumulado de error dividido entre el número de salidas (usado en evaluarpob para calculo de fitness total del genoma)
    pGenoma->fitness += acum;//TODO: verificar si funciona quitando la división por el número de salidas. antes era(acum/(float)conf->numSalidas)
    return(1);
}
//TODO: Falta adicionar a lista de hijos de padre en funciones de nuevoNodo y NuevaConex, verificar si las dos mutaciones funcionan bién con eso
//	también verificar el genoma inicial y marcar como valorCalculado=1 y numHijos=0 a las neuronas de entrada.
//	también verificar en crossover si los nodos del hijo quedan con la misma lista de nodos hijo y conexHijo.

unsigned genomaPerfecto(unsigned index, TConfig* conf)
{
//Coloca un genoma perfecto de xor en la posición deseada. retorna 0 si hay error, 1 swi ok
//parametros :index
    //genera genoma inicial (todos los th en 0.5)
    if (genomaInicial(index,conf->numEntradas,conf->numSalidas,conf->numBias,0,0,conf)==0)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 34 en Funcion genomaPerfecto() llamando a genomaInicial()\n");
        return(0);
    }
    //agrega nuevo nodo 4 entre 2 y 3
    if((nuevoNodo(index,2,conf))==0) ////TODO CURVAS DE I2 para seleccionar i
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 19 en funcion genomaPerfecto() llamando a funcion nuevoNodo()\n");
        return(0);
    }
    //hace peso=0 de última conex no se necesita bias
    conf->pob[index].conex[conf->pob[index].totalConexiones-1].peso=0;
    //agrega nueva conex entre 0 y 4 peso=1
    if(nuevaConex(index,0,4, 1.0, 0, 1,conf)==0)
    {
        //Si hubo error retona 0
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\Error 13 en funcion genomaPerecto() llamando a nuevaConex()\n");
        return(0);
    }
    //agrega conexión entre 1 y 4 peso =-1
    if(nuevaConex(index,1,4, -1.0, 0, 1,conf)==0)
    {
        //Si hubo error retona 0
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\Error 13 en funcion genomaPerecto() llamando a nuevaConex()\n");
        return(0);
    }
    //establece peso de -1 para la conex 0 (0->3)
    conf->pob[index].conex[0].peso=-1.0;
    //establece peso de 2 para la conexión 3 (de 4 a 3)
    conf->pob[index].conex[3].peso=2.0;
    //TODO: hacer función de preanálisis de entradas y salidas para establecer la ganancia de normalización de las entradas y de denormalización de salidas
    //TODO: la función ed preanálisis  también calcula los thresholds de las neuronas de entrada (para las demás es inicialmente 0.5)
    //TODO: verificar cuál es el mejor threshold inicial ara ocultas y salida? usando evaluaciones de fsigma (es el máximo o la mitad del unsigned  intervaloinicialmente?)
    return(1);
}

unsigned mutarAC(unsigned indexpob, unsigned maxIntentos, TConfig* conf)  //OPTIMIZADA
{
//Agrega una conexión al azar usando la funcion nuevaconexión
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
    unsigned i=0;
    unsigned j=0;
    unsigned k=0;
    unsigned l=0;
    unsigned innovIn=0;
    unsigned innovOut=0;
    float peso;
    //verifica que la conexión no exista:
    // para k=0 hasta maxIntentos:
    for (k=0; k<maxIntentos; k++)
    {
        // selecciona dos nodos random i,j menores a maxNodos El límite de neuronas es 32767 (RAND_MAX), se puede aumentar usando rand*rand/1073676289
        i=(unsigned  int)floor((((float)randL(conf))*(conf->pob[indexpob].totalNodos-1))+0.5);
        // NO se permiten conexiones a nodos bias como destino de la conexión
        do
        {
            j=(unsigned  int)floor((((float)randL(conf))*(conf->pob[indexpob].totalNodos-1))+0.5);
        }
        while(conf->pob[indexpob].nodo[j].nodeFunction==3);
        innovIn=conf->pob[indexpob].nodo[i].innovNum;
        innovOut=conf->pob[indexpob].nodo[j].innovNum;
        // si buscarInnovConexPorNodos == UINT_MAX
        if((l=buscarInnovConexPorNodos(indexpob,innovIn,innovOut,conf))==UINT_MAX)
        {
            // se debe además verificar que en casos de conexiones recurrentes, el nodo no sea de entrada o bias.
            // 0=entrada,1=oculto,2=salida,3=bias.
            if (i==j)
            {
                if ((conf->pob[indexpob].nodo[i].nodeFunction==1)||(conf->pob[indexpob].nodo[i].nodeFunction==2))
                {
                    // adiciona la nueva conexión con la funcion:
                    // unsigned nuevaConex(indexpob,innovIn,innovOut, peso, recurrente, 1);//retorna 0 si hay error
                    peso=(((float)randL(conf))*4)-2; // (-2,2)peso aleatorio para la nueva conex (-2,2)
                    if(nuevaConex(indexpob,i,j, peso, 1, 1,conf)==0)
                    {
                        //Si hubo error retona 0
                        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\Error 13 en funcion mutarAC(%u,%u) llamando a nuevaConex(%u,%u,%u,%u,%u,%u)\n",indexpob,maxIntentos,indexpob,innovIn,innovOut, 1, 1, 1);
                        return(0);
                    }
                    // retorna 1 si OK
                    /*	if (erificarGenoma(indexpob)==0){
                    		fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 14.4 en funcion MutarAC(%u,%u) llamando a verificarGenomas()\n",indexpob,maxIntentos);
                    		return(0);
                    	}*/
                    return(1);
                }
            }
            else
            {
                peso=(((float)randL(conf))*conf->pob[indexpob].nodo[j].contHijos*4)-(2*conf->pob[indexpob].nodo[j].contHijos); // (-1,1) por cada conthijo peso aleatorio para la nueva conex (-2,2)
                if(nuevaConex(indexpob,i,j, peso, 0, 1,conf)==0)
                {
                    //Si hubo error retona 0
                    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 14 en funcion mutarAC(%u,%u) llamando a nuevaConex(%u,%u,%u,%u,%u,%u)\n",indexpob,maxIntentos,indexpob,innovIn,innovOut, 1, 0, 1);
                    return(0);
                }
                // retorna 1 si OK
                return(1);
            }
        }
        else
        {
            if(conf->pob[indexpob].conex[l].enabled==0)
            {
                if(((float)randL(conf))<conf->porcentEnableds)
                {
                    conf->pob[indexpob].conex[l].enabled=1;
                    //incrementa conCont del nodoout de la conex
                    calcularLimThOneNode(indexpob,conf->pob[indexpob].conex[l].nodoOut,conf);
                    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>R(%u)",l);
                    return(1);
                }
            }
        }
    }
    return(2); //No se pudo encontrar una conexión no existente.
}

unsigned mutarAN(unsigned indexpob, TConfig* conf)  //OPTIMIZADA
{
//Agrega un nodo al azar usando la funcion nuevoNodo, entran en la selección todas las conexiones existentes.
//reerva memoria para el nuevo tamaño del genoma con realloc y luego adiciona el nuevo nodo.
//¿Es posible agregar un nodo bias? si, es necesario?
// Retorna 0 si hubo error , 1 si ok.
// Parámetros:	indexpob	= index de conf->pob que se desea mutar.
    float i=0;
    ////TODO:Sustituir todos los rand por (float)(float)rand()
    //Crea el nuevo nodo con la fución nuevoNodo(unsigned indexpob, unsigned indexElimInnovConex )
    if (conf->pob[indexpob].totalConexiones>conf->minConexMutPeso)
    {
        //Selecciona al azar un número i entre 0 y maxConex-1 (como rand solo va hasta 32k, se usa rand*rand)
        i=(float)randL(conf);
        i*=((float)conf->pob[indexpob].totalConexiones-1);
        if((nuevoNodo(indexpob,(unsigned  int)i,conf))==0) ////TODO CURVAS DE I2 para seleccionar i
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 19 en funcion mutarAN(%u) llamando a funcion nuevoNodo(%u,%u)\n",indexpob,indexpob,(unsigned  int)i);
            return(0);
        }
        return (1);
    }
    return (1);
}

unsigned nuevaConex(unsigned indexpob,unsigned indexIn, unsigned indexOut, float peso, short unsigned recurrente, short unsigned enabled, TConfig* conf) //OPTIMIZADA TODO: Va en gen.c no aquí.
{
    /******************************************/
    /*       funcion  nuevaConex()             */
    /******************************************/
//Adiciona o cobreescribe unaonexión a un genoma de la conf->población (param = indice del genoma en la conf->pob, todas
//las variables de struct GenConexF excepto innovNum que es una variable global para cada gen.).
//además, incrementa el contador de hijos del nodo indexOut
//Parámetros: indexpob, IndexIn, IndexOut,recurrente, peso, enabled.
//Retorna 0 si hubo error, 1 si ok.
////TODO MODIFICAR LA FUNCION EVALUAR GENOMA PARA EVALUAR POR INNOVNUM DE NODOS y CONEXIONES en lugar de indexes.
////TODO: Verificar si genoma adicionado ya existe en el genoma, si existe, le adiciona su peso. y tiene 50% de chance de ser enabled si
//			una de ellas (existente o nueva)es disabled(manejar conCont del nodo destino).
    //obtener memoria con realloc en conf->pob[indexpob].conex para un tamaño (totalConexiones+1)sizeof(GenNodoF)
    Genoma* pGenoma = &(conf->pob[indexpob]); //usado para acelerar operaciones en el genoma
    GenConexF* pConex; //puntero a la nueva conex creada, usado para acelerar operaciones en el nodo
    unsigned innovIn = pGenoma->nodo[indexIn].innovNum;
    unsigned innovOut = pGenoma->nodo[indexOut].innovNum;
    unsigned i;
    // verifica si la conex ya existe.
    for (i=0;i<conf->pob[indexpob].totalConexiones;i++)
    {
        if ((conf->pob[indexpob].conex[i].indexIn==indexIn)&&(conf->pob[indexpob].conex[i].indexOut==indexOut))
        {
            // si existe, sobreescribe los valores y retorna OK.
            conf->pob[indexpob].conex[i].peso=peso;
            conf->pob[indexpob].conex[i].enabled=enabled;
            conf->pob[indexpob].conex[i].recurrente = ((indexIn==indexOut)? 1: 0);
            return(1);
        }
    }
    // ubica memoria para la nueva conexión
    if ((pGenoma->conex=(GenConexF *)realloc(pGenoma->conex, sizeof(GenConexF)*(pGenoma->totalConexiones+1)))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 11 en funcion nuevaConex(%u,%u,%u,%1.1f,%u,%u) llamando a realloc(conf->pob[%u].conex,%u)\n",indexpob,indexIn, indexOut, peso, recurrente, enabled,indexpob,(pGenoma->totalConexiones+1)*(unsigned int)sizeof(GenConexF));
        return(0);
    }
    pConex=&(pGenoma->conex[pGenoma->totalConexiones]);
    //Inicializa los valores de la nueva conexión
    if ((pConex->innovNum=nuevaInnovCon(innovIn,innovOut,conf))==UINT_MAX)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 12 en funcion nuevaConex(%u,%u,%u,%1.1f,%u,%u) llamando a nuevaInnovCon(%u,%u)\n",indexpob, indexIn, indexOut, peso,recurrente,enabled,innovIn,innovOut);
        return(0);
    }
    //La sig. linea mantiene el innovNum máximo para el genoma indexPop
    if (pGenoma->maxInnovNumConex<pConex->innovNum) pGenoma->maxInnovNumConex=pConex->innovNum;
    pConex->nodoIn=innovIn;
    pConex->nodoOut=innovOut;
    pConex->indexOut=indexOut;
    pConex->indexIn=indexIn;
    pConex->peso=peso;
    pConex->recurrente= ((indexIn==indexOut)? 1: 0);
    pConex->enabled=enabled;
    if (enabled==1)
    {
        calcularLimThOneNode(indexpob,indexOut,conf);
    }
    pGenoma->nodo[indexOut].contHijos++;
    //Incrementa el total de nconexiones  en pGenoma->totalConexiones++;
    pGenoma->totalConexiones++;
    return(1);
}

unsigned nuevoNodo(unsigned indexpob, unsigned indexElimConex, TConfig* conf)  //OPTIMIZADA TODO: Va en gen.c no aquí.
{
// Adiciona un nuevo nodo al genoma unsigned  interrumpiendo la conexión indicada y conservando el peso original en la conex nuevo a out y
// colocando peso=1 para la conex in a nuevo y coloca en disabled la conexión original.
// Parámetros : 	indexpob = indice de la pob para el genoma al que se adicionará el nodo
// 					indexElimInnovConex = index del nodo en el arreglo de nodos.
// Retorna 0 si hubo algún error, 1 si ok.
// Coloca en disabled la conexión con número de innovación indexElimConex
    Genoma* pGenoma = &(conf->pob[indexpob]); //usado para acelerar operaciones en el genoma
    GenNodoF* pNodo; //puntero al nuevo nodo creado, usado para acelerar operaciones en el nodo
    pGenoma->conex[indexElimConex].enabled=0;
    // obtener memoria con realloc en pGenoma->nodo para un tamaño (totalNodos+1)sizeof(GenNodoF)
    if ((pGenoma->nodo=(GenNodoF *)realloc(pGenoma->nodo, sizeof(GenNodoF)*(pGenoma->totalNodos+1)))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 15 en fupGenoma->nodo[pGenoma->totalNodos].minTh=0;ncion nuevoNodo(%u,%u) llamando a realloc(pGenoma->nodo,%u).",indexpob,indexElimConex,(pGenoma->totalNodos+1)*(unsigned int)sizeof(GenNodoF));
        return(0);
    }
    //Asigna pNodo al nodo recién ubicado en memoria para acelerar calculos.
    pNodo = &(pGenoma->nodo[pGenoma->totalNodos]);
    // inicializa los valores del nuevo nodo. obteniendo el número de innovación para New.
    if ((pNodo->innovNum=nuevaInnovNodo(pGenoma->conex[indexElimConex].nodoIn,pGenoma->conex[indexElimConex].nodoOut,conf))==UINT_MAX)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 16 en funcion nuevoNodo(%u,%u) en funcion nuevaInnovNodo(%u,%u)",indexpob,indexElimConex,pGenoma->conex[indexElimConex].nodoIn,pGenoma->conex[indexElimConex].nodoOut);
        return(0);
    }
    pNodo->nodeFunction=1; //1=oculto
    pNodo->estadoC=0;
    pNodo->valor=0.5;
    pNodo->maxTh=1;
    pNodo->minTh=0;
    pNodo->contHijos=0;
    pNodo->conexHijo=NULL;
    // TODO: que valor de threshold debería tener el nuevo nodo?. Por ahora se coloca en el medio del límite de maxTh y minTh.
    pNodo->thNodo=conf->Fthreshold;
    if (conf->tSigma>1000)  //si tsigma es -1,1 //TODO: verificar sies necesario esto si se usa mutar th,
    {
        pNodo->minTh=-1.0;
        pNodo->valor=0;
    }
    // mantiene en para el genoma el máximo número de innovación de nodo.
    if (pGenoma->maxInnovNumNodo<pNodo->innovNum)
    {
        pGenoma->maxInnovNumNodo=pNodo->innovNum;
    }
    // inicializa los valores de la nueva conexión New->Out co peso  igual a la conexión eliminada y con el valor de enabled de la anterior?????//TODO
    if(nuevaConex(indexpob,pGenoma->totalNodos, pGenoma->conex[indexElimConex].indexOut, pGenoma->conex[indexElimConex].peso  , 0, 1,conf)==0)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 17 en funcion nuevoNodo(%u,%u) en funcion nuevaConex(%u,%u,%u,%1.1f,%u,%u)",indexpob,indexElimConex,indexpob,pGenoma->totalNodos, pGenoma->conex[indexElimConex].indexOut, pGenoma->conex[indexElimConex].peso  , 0, 1);
        return(0);
    }
    // inicializa los valores de la nueva conexión In->New con peso 1
    if(nuevaConex(indexpob,pGenoma->conex[indexElimConex].indexIn,pGenoma->totalNodos,1, 0, 1,conf)==0)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 18 en funcion nuevoNodo(%u,%u) llamando a funcion nuevaConex(%u,%u,%u,%u,%u,%u)",indexpob,indexElimConex,indexpob,pGenoma->conex[indexElimConex].indexIn,pGenoma->totalNodos,1 , 0, 1);
        return(0);
    }
    pGenoma->totalNodos++;
    return(1);
}

unsigned genomaInicial(unsigned index,unsigned nEntradas, unsigned nSalidas, unsigned nBias, unsigned primero, unsigned especie,TConfig* conf)  // OTIMIZADA TODO: optimizar con punteros
{
//Crea un nuevo genoma totalmente conectado asigna la especie 0 y lo ubica en el índice index de la población.
//parámetros: index=donde queda el genoma inicial, nEntradas, nSalidas, nBias.
//			primero si =1 se bora lista de innovaciones , representantes, conservacion, generaciones sin mejora, etc...
//retorna 0 si hubo error, 1 si la creación fué exitosa.
    unsigned i=0;
    unsigned j=0;
    unsigned k=0;
    ////TODO
    //Crear todos los nodos de entrada y salida (a mano).
    //tenemos puntero a estructura Genoma, pero tenemos variables de genoma sin inicializar.
    //Se hace calloc para el arreglo de nodos de tamaño bias+entradas+salidassizeof(nodo)
    if (conf->pob[index].nodo!=NULL) free((void *)conf->pob[index].nodo);
    if (!(conf->pob[index].nodo=( GenNodoF *) calloc(1,sizeof( GenNodoF)*(nEntradas+nSalidas+nBias))))
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 20 en funcion genomaInicial(%u,%u,%u) llamando a calloc(1,%u).\n",nEntradas,nSalidas,nBias,(nEntradas+nSalidas+nBias)*(unsigned int)sizeof( GenNodoF));
        return(0);
    }
    //Para cada  nodo de salida, crea una conexión a cada nodo de entrada y bias con peso 1.
    //Obtiene puntero para conexiones con tamaño nSalidas*(nentradas+nBias)sizeof(GenConexF)
    if (conf->pob[index].conex!=NULL) free((void *)conf->pob[index].conex);
    if ((conf->pob[index].conex=( GenConexF *) calloc(1,sizeof( GenConexF)*nSalidas*(nEntradas+nBias)))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 21 en funcion genomaInicial(%u,%u,%u) llamando a calloc(1,%u).\n",nEntradas,nSalidas,nBias,nSalidas*(nEntradas+nBias)*(unsigned int)sizeof( GenConexF));
        return(0);
    }
    if (primero==1)
    {
        //Se reserva memoria para el arreglo de innovacione para nodos.
        if (!(conf->listaInnovNodo=( TListaInnov *) malloc(sizeof( TListaInnov)*(nEntradas+nSalidas+nBias))))
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 22 en funcion genomaInicial(%u,%u,%u) llamando a malloc(%u).\n",nEntradas,nSalidas,nBias,(nEntradas+nSalidas+nBias)*(unsigned int)sizeof( TListaInnov));
            return(0);
        }
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>PrimerGenoma\n");
        //inicializa los punteros de arreglos de nodos de innovs :
        for(i=0; i<(nEntradas+nSalidas+nBias); i++)
        {
            conf->listaInnovNodo[i].numOut=1;
            if ((conf->listaInnovNodo[i].nodoOut=(TNodoOut*)malloc((unsigned  int)sizeof(TNodoOut)))==NULL)
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 22.1 en funcion genomaInicial(%u,%u,%u) llamando a malloc(%u).\n",nEntradas,nSalidas,nBias,(unsigned  int)sizeof(TNodoOut));
                return(0);
            }
            conf->listaInnovNodo[i].nodoOut[0].innovNum=i;
            conf->listaInnovNodo[i].nodoOut[0].nodoOut=i;
        }
        conf->contInnovNodo=(nEntradas+nSalidas+nBias);
        //Se reserva memoria para el arreglo de innovacione para conexiones.
        if (!(conf->listaInnovCon=( TListaInnov *) malloc(sizeof(TListaInnov)*(nEntradas+nSalidas+nBias))))
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 23 en funcion genomaInicial(%u,%u,%u) llamando a calloc(1,%u).\n",nEntradas,nSalidas,nBias,nSalidas*(nEntradas+nBias)*(unsigned int)sizeof(TListaInnov));
            return(0);
        }
        //inicializa los punteros de arreglos de conexiones de innovs:
        k=0;
        for (i=0; i<(nEntradas+nBias); i++)
        {
            conf->listaInnovCon[i].numOut=nSalidas;
            if((conf->listaInnovCon[i].nodoOut=(TNodoOut*) malloc(sizeof(TNodoOut)*nSalidas))==NULL)  //porque es totalmente conectado
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 23.1 en funcion genomaInicial(%u,%u,%u) llamando a malloc(%u).\n",nEntradas,nSalidas,nBias,(unsigned  int)sizeof(TNodoOut));
                return(0);
            }
            for (j=0; j<nSalidas; j++)
            {
                conf->listaInnovCon[i].nodoOut[j].nodoOut=(nEntradas+nBias)+j;
                conf->listaInnovCon[i].nodoOut[j].innovNum=k;
                k++;
            }
        }
        //completa los valores de los demás elementos de entrada de listaInnovCon
        for (i=(nEntradas+nBias); i<(nEntradas+nBias+nSalidas); i++)
        {
            conf->listaInnovCon[i].numOut=0;
            conf->listaInnovCon[i].nodoOut=NULL;
        }
        conf->contInnovCon=nSalidas*(nEntradas+nBias);
// IMPRIME LA LISTA DE INNOVACIONES DE NODO Y CONEX
//        imprimirListasInnov(conf);
        //reserva memoria para el arreglo de conf->representantes de tamaño 1 (trivial)
        if (conf->representantes!=NULL) free((void *)conf->representantes);
        if  ((conf->representantes=(unsigned *) calloc(1,(unsigned  int)sizeof(unsigned  int)))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 24 en funcion genomaInicial(%u,%u,%u) llamando a calloc(1,%u).\n",nEntradas,nSalidas,nBias,(unsigned  int)sizeof(unsigned  int));
            return(0);
        }
        //reserva memoria para el arreglo de especies en conservación de temaño conf->maxEspeciesConservacion
        if (conf->conservacionEsp!=NULL) free((void *)conf->conservacionEsp);
        if  ((conf->conservacionEsp=(unsigned *) calloc(conf->spEspecies,(unsigned  int)sizeof(unsigned  int)))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 24.4 en funcion genomaInicial(%u,%u,%u) llamando a calloc(maxEspeciesConservación,%u).\n",nEntradas,nSalidas,nBias,(unsigned  int)sizeof(unsigned  int));
            return(0);
        }
        //reserva memoria para el arreglo de numero de genomas por especie, tamaño=conf->numEspecies
        if (conf->conservacionEsp!=NULL) free((void *)conf->conservacionEsp);
        if  ((conf->numGenomasPorEspecie=(unsigned *) calloc(conf->spEspecies,(unsigned  int)sizeof(unsigned  int)))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 24.4 en funcion genomaInicial(%u,%u,%u) llamando a calloc(maxEspeciesConservación,%u).\n",nEntradas,nSalidas,nBias,(unsigned  int)sizeof(unsigned  int));
            return(0);
        }
        //reserva memoria para conf->contGeneracSinMejora
        if (conf->contGeneracSinMejora!=NULL) free((void *)conf->contGeneracSinMejora);
        if  ((conf->contGeneracSinMejora=(unsigned *) calloc(1,(unsigned  int)sizeof(unsigned  int)))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 24.5 en funcion genomaInicial(%u,%u,%u) llamando a calloc(1,%u).\n",nEntradas,nSalidas,nBias,(unsigned  int)sizeof(unsigned  int));
            return(0);
        }
    }
    //Se inicializa la funcion y el valor en los correspondientes valores iniciales
    //para entradas:
    for (i=0; i<nEntradas; i++)
    {
        conf->pob[index].nodo[i].nodeFunction=0; //0=entrada,1=oculto,2=salida,3=bias.
        conf->pob[index].nodo[i].valor=0; //Valor de salida de cada nodo(para computar la ann sin matriz de pesos).
        conf->pob[index].nodo[i].estadoC=1; // no se calcula, se tiene de las entradas.
        conf->pob[index].nodo[i].maxTh=1; //TODO: calcular rango de thresholds con preanálisis
        conf->pob[index].nodo[i].minTh=0;
        conf->pob[index].nodo[i].contHijos=0;
        conf->pob[index].nodo[i].conexHijo=NULL;
        // TODO: que valor de threshold debería tener el nuevo nodo?. Por ahora se coloca en el medio del límite de maxTh y minTh.
        conf->pob[index].nodo[i].thNodo=conf->Fthreshold;
        if (conf->tSigma>1000)  //si tsigma es -1,1
        {
            conf->pob[index].nodo[i].minTh=-1.0;
            conf->pob[index].nodo[i].valor=0;
        }
        //coloca valores de listainnov para los nodos
        if ((conf->pob[index].nodo[i].innovNum=nuevaInnovNodo(i,i,conf))==UINT_MAX)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 25.1 en funcion genomaInicial(%u,%u,%u) llamando a nuevaInnovNodo(%u,%u)\n",nEntradas,nSalidas,nBias,i,i);
            return(0);
        }
    }
    //para bias:
    for (i=nEntradas; i<(nEntradas+nBias); i++)
    {
        conf->pob[index].nodo[i].nodeFunction=3; //0=entrada,1=oculto,2=salida,3=bias.
        conf->pob[index].nodo[i].valor=1; //Ya que es bias, se inicializa en 1
        conf->pob[index].nodo[i].estadoC=1; //nunca se calcula
        conf->pob[index].nodo[i].maxTh=1; //TODO: calcular rango de thresholds con preanálisis
        conf->pob[index].nodo[i].minTh=0;
        conf->pob[index].nodo[i].contHijos=0;
        conf->pob[index].nodo[i].conexHijo=NULL;
        // TODO: que valor de threshold debería tener el nuevo nodo?. Por ahora se coloca en el medio del límite de maxTh y minTh.
        conf->pob[index].nodo[i].thNodo=conf->Fthreshold;
        if (conf->tSigma>1000)  //si tsigma es -1,1
        {
            conf->pob[index].nodo[i].minTh=-1.0;
        }
        if ((conf->pob[index].nodo[i].innovNum=nuevaInnovNodo(i,i,conf))==UINT_MAX)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 25.2 en funcion genomaInicial(%u,%u,%u) llamando a nuevaInnovNodo(%u,%u)\n",nEntradas,nSalidas,nBias,i,i);
            return(0);
        }
    }
    //para salidas:
    for (i=(nEntradas+nBias); i<(nEntradas+nBias+nSalidas); i++)
    {
        conf->pob[index].nodo[i].nodeFunction=2; //0=entrada,1=oculto,2=salida,3=bias.
        conf->pob[index].nodo[i].valor=0; //Valor de salida de cada nodo(para computar la ann sin matriz de pesos).
        conf->pob[index].nodo[i].estadoC=0;
        conf->pob[index].nodo[i].maxTh=1; //TODO: calcular rango de thresholds con preanálisis
        conf->pob[index].nodo[i].minTh=0;
        conf->pob[index].nodo[i].contHijos=nEntradas+nBias;
        conf->pob[index].nodo[i].conexHijo=NULL;
        // TODO: que valor de threshold debería tener el nuevo nodo?. Por ahora se coloca en el medio del límite de maxTh y minTh.
        conf->pob[index].nodo[i].thNodo=conf->Fthreshold;
        if (conf->tSigma>1000)  //si tsigma es -1,1
        {
            conf->pob[index].nodo[i].minTh=-1.0;
            conf->pob[index].nodo[i].valor=0;
        }

        if ((conf->pob[index].nodo[i].innovNum=nuevaInnovNodo(i,i,conf))==UINT_MAX)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 25.3 en funcion genomaInicial(%u,%u,%u) llamando a nuevaInnovNodo(%u,%u)\n",nEntradas,nSalidas,nBias,i,i);
            return(0);
        }
    }
    //Se inicializa en los valores necesarios.
    k=0;
    for (j=(nEntradas+nBias); j<(nEntradas+nBias+nSalidas); j++)
    {
        for (i=0; i<(nEntradas+nBias); i++)
        {
            conf->pob[index].conex[k].nodoIn = i;
            conf->pob[index].conex[k].nodoOut = j;
            conf->pob[index].conex[k].indexIn = i;
            conf->pob[index].conex[k].indexOut = j;
            conf->pob[index].conex[k].peso = 1;
            conf->pob[index].conex[k].recurrente = 0;
            conf->pob[index].conex[k].enabled = 1;
            ////TODO: funcion para obtener automáticamente el número de innovación correspondiente a la conexión
            //basado en la lista de innovaciones y en los nodos de entrada y salida yagregarlo a la lista sino existe.
            //¿es posible hacer otra lista de innovaciones de nodos que lleven en su estructura
            //el innovConNum de entradas y salidas? podría esto verificarse al momento de hacer el cruce únicamente?
            //¿Puede codificarse en el nombre del nodo, el número de conexiones adyacientes?

            //Al aparecer un nuevo nodo su número(nombre) de innovación puede obtenerse  si se compara con
            //una lista donde se lleva para cada nuevo nodo en la generación actual
            //el nodo de origen y destino de la conexión que se unsigned  interrumpió para crearlo.
            //y ese número se usa como identificador del nodo en la estructura GenConexF y en la lista de innovaciones.
            ////TODO:Crear estructuras y modificaciones para usar esto.
            //Para esta primera generación, como los nodos son conocidos, su número de innovación es el mismo número de nodo.
            //y por tanto,:
            if ((conf->pob[index].conex[k].innovNum = nuevaInnovCon(i,j,conf))==UINT_MAX)
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 26 en funcion genomaInicial(%u,%u,%u) en funcion nuevaInnovCon(%u,%u) \n",nEntradas,nSalidas,nBias,i,j);
                return(0);
            }
            //Coloca el número de innovación correspondiente a cada conexión
            k++;
        }
    }
    conf->pob[index].maxInnovNumConex=(nSalidas*(nEntradas+nBias))-1;
    conf->pob[index].maxInnovNumNodo=nSalidas+nEntradas+nBias-1;
    //Inicializa los valores restantes del genoma Inicial
    conf->pob[index].especie = especie; //establece la especie inicial.

    if (primero==1)
    {
        conf->numEspecies=1; //adiciona una especie especie (especie 0)
        conf->contGeneracSinMejora[especie]=0;//coloca 0 generaciones sin mejor aen fitness para la especie 0
        //reserva memoria para el nuevo arreglo de especies con tamaño conf->numEspeciessizeof(unsigned  int)
        //reserva memoria para el nuevo arreglo de conf->contGeneracSinMejora
        //coloca al indexpob como representante de la nueva especie
        conf->representantes[conf->numEspecies-1]=index;
        //coloca en 0 el número de generaciones sin mejora para la nueva especie
        conf->contGeneracSinMejora[conf->numEspecies-1]=0;

    }
    conf->pob[index].totalNodos = nEntradas+nSalidas+nBias;
    conf->pob[index].totalConexiones = nSalidas*(nEntradas+nBias);
    conf->pob[index].numHijos=-1;
    conf->pob[index].fitness = 0; // error inicial=100% también inicializa el contador de fitness
    calcularLimThOneGenome(index,conf);
    return(1);
}

void genomaMasLejano(int indexPob,int intentos, TConfig* conf)
// busca el genoma más lejano a todos los de su especie entre la población  <intentos> veces y lo coloca en indexPob
// también debe tener distancia menor a 1/2 de distacia a especie más cercana.
{
    int i,j;
    int miEspecie=conf->pob[indexPob].especie;
    unsigned offs=(2*conf->sizePob)+conf->spEspecies;
    float distEspCercana=999999; //si hay una sola especie, se usa este valor
    float distancia[1000];
    float minDist=9999999;
    float increm=(1-conf->porcentMutPeso)/intentos;
    float tmp;
    unsigned foundIndex=0;
    // calcula distancia a la especie más cercana
    distEspCercana=distEspecieCercana(indexPob,conf->c1, conf->c2, conf->c3, conf->eG_t,conf)/2;
    for (i=0;i<intentos;i++) //inicializa todas las mínimas distancias en 0: necesario al final de esta func
    {
        distancia[i]=0;
    }
    // copiarGenoma representante[miEspecie] a indexpob (SE DEBE HACER?)
    for (i=0;i<intentos;i++)
    {
        // copia representante[miEspecie] a los n intentos
        copiarGenoma(conf->representantes[miEspecie],offs+i,conf);
        // aplica perturbación cada vez mayor de pesos empezando en conf->porcentMutPeso terminando terminando en 1
        perturbarPeso(offs+i,increm*(float)i,conf->probMutPeso,0,0,0,0, conf);
        //randomizarPesos(offs+i,conf->pesoMinInicial,conf->pesoMaxInicial,conf);
        //inicializa minDist
        minDist=999999;
        // para cada intento para cada genoma de la pob:
        for (j=0;j<conf->sizePob;j++)
        {
            //solo si es de la misma especie que indexpob
            if (conf->pob[j].especie==miEspecie)
            {
                if (j!=indexPob)
                {
                    // medir distancia a cada uno de los de mi especie exc yo y guardar mínimo en distancia[i]
                    tmp=calcularDist(offs+i, j, conf->c1, conf->c2, conf->c3, conf->eG_t,conf);
                    if ((tmp<minDist)&&(tmp>0))
                    {
                        minDist=tmp;
                    }
                }
            }
        }
        //guarda la mínima distancia encontrada en la especie en distancia[i]
        distancia[i]=minDist;

    }

    tmp=0;
   // copiar máximo entre los temp[] a indexPob
    for(i=0;i<intentos;i++)
    {
        //si la distancia es máxima y menor a un medio de la distancia entre reps más cercanos.
        if ((distancia[i]>tmp)&&(calcularDist(offs+i, conf->representantes[miEspecie], conf->c1, conf->c2, conf->c3, conf->eG_t,conf)<distEspCercana)){
            tmp=distancia[i];
            foundIndex=i;
        }
    }
    copiarGenoma(offs+foundIndex,indexPob,conf);
}

//No quiero perderme ni un fotón de tus refléjos. :] Juli.


