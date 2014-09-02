/** Archivo de manejo de poblaci�n, incluye ciclo principal - C file
	usan/modifican la poblacion �nicamente (pueden usarse operaciones de especies, genomas o genes).
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
#ifndef GENOMA_H_INCLUDED
#include "genoma.h"
#define GENOMA_H_INCLUDED
#endif
#ifndef ESPECIE_H_INCLUDED
#include "especie.h"
#define ESPECIE_H_INCLUDED
#endif
#include "pob.h"
#include <cuda.h>

#define THREADSPB 96
#define MAXCONEX 1000

//kernel de cuda que toma datos de entrenamiento y lista de conexiones, genera vector de fitness.
__global__ void evalPobCUDA(int numEntradas, int numSalidas, int numBias, int numDatos, float* dataGTDf, int** conexIn, int** conexOut, float** conexPeso, int* tamListaConexPost, float** valorC, float** valorTr, float* fitness)
{
	int i,j,k;
    // vectores en registros para conexiones: valor[], conexIn[] conexOut[] conexPeso[]
    int r_conexIn[16];
    int r_conexOut[16];
    float r_conexPeso[16];
	// vector en registros para los nodos (HASTA 96 ENTRADAS,6 salidas,192 Nodos en total)
    float valor[16];
    // vector de salida de entrenamiento.
    float r_valorTr[4];
	// vector de tama�os de lista de conexiones.
	int r_tamListaConexPost;
    // calcula el id del thread
    int idx=blockDim.x*blockIdx.x+threadIdx.x;
	// calcula el n�mero de repeticiones para floor de todos los datos/buffSize
    int numReps=(numDatos*(numEntradas+numSalidas))/THREADSPB;
	// THREADSPB tambi�n es el n�mero de datos completos leidos de global a shared.
	// vector de datos de entrenamiento leidos en shared de tama�o THREADSPB
	__shared__ float s_dataGTDf[THREADSPB];
	// vector de fitness resultante para cada genoma por repetici�n entre grbaciones a global.
	__shared__ float s_fitness[THREADSPB];
    // vector de valor calculado para para calculo de fitness por correlaci�n.
	__shared__ float s_valorC[THREADSPB];
    // vector de valor de entrenamiento para para calculo de fitness por correlaci�n.
	__shared__ float s_valorTr[THREADSPB];

    // inicializa el fitness del genoma  en 0
    s_fitness[idx]=0;
    // calcula l�mites de for para optimizaci�n.
    int indSalidas=(numEntradas+numBias+numSalidas);
    int indBias=(numEntradas+numBias);
    int indInS=numEntradas*4;
    int g_indexGTD=0;
    int s_indexGTD=0;
    int numRepsShared=THREADSPB/(numEntradas+numSalidas);
    float A =-2.435;
    float acum=0;
    float vCm=0; //para las medias
    float vTm=0;
    int contDato=0;
    float sum0=0;
    float sum1=0;
    float sum2=0;
	int r_tamListaConexPostAnt;
    __syncthreads();
	// coloca en 0 valor[]
	for (i=0;i<16;i++)
    {
        valor[i]=0;
    }

	// vector de tama�os de lista de conexiones.
    __syncthreads();
    r_tamListaConexPost=tamListaConexPost[idx];
	r_tamListaConexPostAnt=r_tamListaConexPost-1;
    // lee las conexiones en los arreglos de registros (LENTO)
    __syncthreads();
	//for (i=0;i<r_tamListaConexPost;i++)
	if (r_tamListaConexPost>16) r_tamListaConexPost=16;
	for (i=0;i<r_tamListaConexPost;i++)
    {
        r_conexIn[i]=conexIn[idx][i];
        r_conexOut[i]=conexOut[idx][i];
        r_conexPeso[i]=conexPeso[idx][i];
    }

	// para el n�mero de repeticiones:
    for (i=0;i<numReps;i++)
    {
		
		// coloca los valores de entrada desde dataGTDf[indexGTD] en el buffer s_dataGTDf[]
        __syncthreads();
        s_dataGTDf[idx]=dataGTDf[g_indexGTD+idx];
       // __syncthreads();
        // para cada uno de los datos en s_dataGTDf[]
        s_indexGTD=0;

		for (j=0;j<numRepsShared;j++)
        {

            // coloca bias en 1 TODO: se puede colocar solo na vez al principio (VERIFICAR).
            valor[numEntradas]=1;
			// coloca entradas en sus respectivos nodos valor[].

			//			for (k=0;k<numEntradas;k++)

			for (k=0;k<numEntradas;k++)
            {
                __syncthreads();
                valor[k]=s_dataGTDf[s_indexGTD];
				s_indexGTD++;
            }
			// coloca las salidas de entrenamiento en el vector valorTr
            //for (k=0;k<numSalidas;k++)
			for (k=0;k<numSalidas;k++)
            {
                __syncthreads();
				r_valorTr[k]=s_dataGTDf[s_indexGTD];
                __syncthreads();
                s_valorTr[s_indexGTD]=r_valorTr[k];
                // acumula para calulo de medias
                vTm=vTm+r_valorTr[k];
				__syncthreads();
                s_indexGTD++;
            }
            // calcula el fSigma de las entradas
            for (k=0;k<numEntradas;k++)
            {
                valor[k]=2*(expf(A*((valor[k]-1)*(valor[k]-1))))-1;
            }

			// inicializa el acumulador.
            acum=0;
            for (k=0;k<r_tamListaConexPost;k++)
            {
                acum=acum+(r_conexPeso[k]*valor[r_conexIn[k]]);
                // si el conexOut no es el �ltimo
				// OPTIMIZABLE: r_tamListaConexPost-1 puede ser una cte o reg
                if (k<(r_tamListaConexPostAnt))
                {
                    // si el ConexOut siguiente es diferente al actual
                    if (r_conexOut[k]!=r_conexOut[k+1])
                    {
                        // calcula el sigma
                        valor[r_conexOut[k]]=2*(expf(A*((acum-1)*(acum-1))))-1;
                        //reinicializa el acumulador
                        acum=0;
                    }
                }
            }
            // saca el Fsigma del �ltimo nodo.// OPTIMIZABLE: r_tamListaConexPost-1 puede ser una cte o reg
            valor[r_conexOut[r_tamListaConexPostAnt]]=2*(expf(A*(acum-1)*(acum-1)))-1;
            // guarda los valores calculados en el arreglo valoresC[]
            s_valorC[idx]=valor[indBias];
            // acumula para calculo de fitness con correlaci�n
            vCm=vCm+valor[indBias];
            // incrementa el index del dato leido
            contDato++;

        }
        //copia s_valorC y s_valorTr a global
        for (j=0;j<numRepsShared;j++)
        {
            //copia s_valorC y s_valorTr a global
            __syncthreads();
            valorTr[idx][g_indexGTD+j]=s_valorTr[idx];
            __syncthreads();
            valorC[idx][g_indexGTD+j]=s_valorC[idx];
        }
        // incrementa index de dataGTDf
        __syncthreads();
        g_indexGTD+=THREADSPB;
	}
    // calcula el fitness como el coeficiente de correlaci�n de Pearson de los 2 vectores (1 si =es,-1 si inversos, 0 si diferentes)
    //fitness[idx]=correlacCUDA(valorC[idx],valorTr[idx]);
    vCm/=contDato;
    vTm/=contDato;
    // calcula sum0(0,nD,(Xi-Xm)*(Yi-Ym)),sum1(0,nD,sqr(Xi-Xm)) y sum2(0,nD,sqr(Yi-Ym)
	for (i=0;i<numDatos;i++)
    {
        // lee en s_valorC y s_valorTr los vectores globales
        __syncthreads();
        s_valorC[idx]=valorC[idx][i];
        __syncthreads();
        s_valorTr[idx]=valorTr[idx][i];
        __syncthreads();
        sum0+=((s_valorTr[idx]-vTm)*(s_valorC[idx]-vCm));
        __syncthreads();
        sum1+=((s_valorTr[idx]-vTm)*(s_valorTr[idx]-vTm));
        __syncthreads();
        sum2+=((s_valorC[idx]-vCm)*(s_valorC[idx]-vCm));
        __syncthreads();
        g_indexGTD+=THREADSPB;
    }
    // retorna la correlaci�n
    __syncthreads();
    fitness[idx]=(sum0/(sqrtf(sum1)*sqrtf(sum2)));
	fitness[idx]=0.5;
    __syncthreads();
}

int evalPob(int numGenomas, int numDatos, hdrGTDv1 headerGTD, hdrSNNv1* headerSNN, float** dataGTD, tConexDataF** listaConexData, float* fitness, int** ordenEval, int* tamListaConexPost, TConfig* conf)
// genera un vector de fitness al evaluar los datos GTD con la lista de conexiones de los genomas SNN
// asume que todos los SNN tienen iguales valores de encabezado excepto el n�mero de conex.
// retorna 0 si hubo error
{
    int i,j;
    int threadsPorBlock=THREADSPB;
    int blocksPorGrid=(1+threadsPorBlock-1)/threadsPorBlock;//?
	float fitness_h[THREADSPB];
	int tempCin[MAXCONEX];
	int tempCout[MAXCONEX];
	float tempCpeso[MAXCONEX];
    // inicializa vectores de evaluaci�n.
    for (i=0;i<numGenomas;i++)
    {
        // reserva memoria para ordenEval
        conf->ordenEval[i]=(int*)malloc(conf->pob[i].totalNodos*sizeof(int));
        if (!conf->ordenEval[i])
        {
            printf("\nError 66.0.7 en genOrdenEvalF1g llamando a malloc()");
            return(0);
        }
        // actualiza el headerSNN;
        genSSNhdr1(i,conf);
        // genera la lista de conex para este genoma
        genOrdenEvalF1g(i, conf);
        // reordena listaConexData para que est�n primero los nodos necesarios y elimina enableds (optimizaci�n)
        conf->tamListaConexPost[i]=ordenarListaConexF(conf->headerSNN[i], conf->listaConexData[i] ,conf->ordenEval[i],conf->tamOrdenEval[i]);
        // CUDA: reserva memoria en device para listaConexData[i][tamConexPost[i]]
        cudaMalloc(&(conf->cu_conexIn[i]),conf->tamListaConexPost[i]*sizeof(int));
		// CUDA: reserva memoria en device para listaConexData[i][tamConexPost[i]]
        cudaMalloc(&(conf->cu_conexOut[i]),conf->tamListaConexPost[i]*sizeof(int));
		// CUDA: reserva memoria en device para listaConexData[i][tamConexPost[i]]
        cudaMalloc(&(conf->cu_conexPeso[i]),conf->tamListaConexPost[i]*sizeof(float));
		// genera vector conexIn para transferir
		for (j=0;j<conf->tamListaConexPost[i];j++) 
		{
			tempCin[j]=conf->listaConexData[i][j].conexIn;
			tempCout[j]=conf->listaConexData[i][j].conexOut;
			tempCpeso[j]=conf->listaConexData[i][j].peso;
		}
        // CUDA:  transfiere listaConexData[i] ordenada a device
        cudaMemcpy(conf->cu_conexIn[i], tempCin, conf->tamListaConexPost[i]*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(conf->cu_conexOut[i], tempCout, conf->tamListaConexPost[i]*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(conf->cu_conexPeso[i], tempCpeso, conf->tamListaConexPost[i]*sizeof(float), cudaMemcpyHostToDevice);
    }
	printf("Llamando kernel de CUDA\n");
    // llama el kernel de CUDA para los 96 threads.
    evalPobCUDA<<<blocksPorGrid,threadsPorBlock>>>(conf->numEntradas,conf->numSalidas,conf->numBias,conf->numDatos,conf->cu_dataGTDf,conf->cu_conexIn,conf->cu_conexOut,conf->cu_conexPeso,conf->cu_tamListaConexPost,conf->cu_valorC, conf->cu_valorTr, conf->cu_fitness);
    // CUDA: transfiere el vector fitness[sizePob] de device a host
    cudaMemcpy(fitness_h, conf->cu_fitness, conf->sizePob*sizeof(float), cudaMemcpyDeviceToHost);
    // libera memoria de ordenEval y de device de listaConexData
    for (i=0;i<numGenomas;i++)
    {
		printf("F%3.3f",fitness_h[i]);
        free(conf->ordenEval[i]);
        // CUDA: libera memoria de  listaconexData ordenada
        cudaFree(conf->cu_conexIn[i]);
		cudaFree(conf->cu_conexOut[i]);
		cudaFree(conf->cu_conexPeso[i]);
    }
    printf("\n");
    return(1);
}

unsigned distribProc(TConfig* conf)
// lee los archivos NNP con el prefix fileNameDistrib y compara cada genoma con los
// representantes actuales, si el m�s cercano de distrib es mejor, lo deja.
// usado para procesamiento distribuido.
// retorna 0 si error, 1 si ok.
{
    int i, j, k;
    char tmpString[512];
    FILE* fileIn=NULL;
    int leidos;
    hdrNNPv1 headerNNP;
    unsigned buscado1;
    unsigned buscado2;
    unsigned especieCercana;
    // copia el fileNameDistrib a tmpstring;
    strcpy(tmpString, conf->fileNameDistrib);
    // para i= el n�mero cargarDistrib
    for (i=0; i<conf->cargarDistrib; i++)
    {
        // genera el nombre de archivo
        sprintf(tmpString, "%s%i.nnp", conf->fileNameDistrib, i);
        // abre el archivo
        fileIn=fopen(tmpString,"rb");
        if (fileIn==NULL)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 1.1 en funcion distribProc() llamando a fopen\n");
            return(1);
        }
        else
        {
            // lee encabezado NNP
            leidos = fread(&headerNNP,sizeof(hdrNNPv1),1,fileIn);
            if (leidos<1)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 1.2 en funcion distribProc() llamando a fread\n");
                return(1);
            }
            else
            {
                // busca primer genoma no representante
                k=0;
                while ((k==conf->representantes[conf->pob[k].especie])&&(k<conf->sizePob))
                    k++;
                buscado1=k;
                if (buscado1>=conf->sizePob)
                {
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>\nError 1.3 en funcion distribProc() genoma no encontrado\n");
                    return(0);
                }
                // guarda buscado1 en tmpgenoma
                copiarGenoma(buscado1,conf->tmpIndexPob,conf);
                // para cada SNN hasta numGenomas
                for (j=0; j<headerNNP.numGenomas; j++)
                {
                    // carga SNN en busca1
                    if (snnDataLoader(buscado1,0,fileIn,conf)==0)
                    {
                        fclose(conf->logFile);
                        conf->logFile=fopen(conf->fileNameLog,"a+");
                        fprintf(conf->logFile,"<br>\nError 1.4 en funcion distribProc() llamando a snnDataLoader(%d)\n",buscado1);
                        return(0);
                    }
                    // busca la especie mas cercana a pob[busca1]
                    especieCercana=especieMinDist(buscado1, conf->c1,conf->c2,conf->c3,conf->eG_t,conf);
                    // actualiza la especie
                    conf->pob[buscado1].especie=especieCercana;
                    // si headerSNN.fitness > representantes[especieCercana].fitness
                    if ((conf->pob[buscado1].fitness>conf->pob[conf->representantes[especieCercana]].fitness)&&(conf->pob[buscado1].fitness>conf->pob[conf->sizePob+especieCercana].fitness))
                    {
                        // busca2 el que NO sea representante y de  especieCercana
                        k=0;
                        while (((k==conf->representantes[conf->pob[k].especie])||(conf->pob[k].especie!=especieCercana)||(k==buscado1))&&(k<conf->sizePob))
                            k++;
                        buscado2=k;
                        // si no lo encontr�, muestra error
                        if (buscado2>=conf->sizePob)
                        {
                            fclose(conf->logFile);
                            conf->logFile=fopen(conf->fileNameLog,"a+");
                            fprintf(conf->logFile,"<br>\nError 1.5 en funcion distribProc()\n");
                            return(0);
                        }
                        // si busca1.fitness>busca2.fitness //nunca se sabe :)
                        if (conf->pob[buscado1].fitness>conf->pob[buscado2].fitness)
                        {
                            // copiarGenoma(busca1,busca2)
                            fclose(conf->logFile);
                            conf->logFile=fopen(conf->fileNameLog,"a+");
                            fprintf(conf->logFile,"D(%d,%d)%d",i,conf->pob[buscado1].especie,buscado2);
                            copiarGenoma(buscado1,buscado2,conf);
                        }
                        // para debugging
                    }
                }
                // restaura tmp a busca1
                copiarGenoma(conf->tmpIndexPob,buscado1,conf);
                // cierra el archivo de entrada.
                fclose(fileIn);
            }
        }
    }
    return(1);
}

unsigned cargarDesdeNNP(char* filename, unsigned especieDestino, unsigned cantidad,TConfig* conf)
// Carga desde un archivo NNP (Neural Network Population) los genomas que se colocan
// en la especie deseada.
// retorna 0 si hubo error, 1 si ok
{
    FILE *fileIn;
    size_t leidos=0;
    hdrNNPv1 headerNNP;
    int i,k;
    // abre el archivo NNP para lectura
    if ((fileIn=fopen(filename,"rb"))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 49 en funcion cargarDesdeNNP() llamando a fopen\n");
        return(0);
    }
    // lee el encabezado NNP
    leidos=fread(&headerNNP,sizeof(hdrNNPv1),1,fileIn);
    // verifica leidos
    if (leidos<1)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>Error 343 en cargarDesdeNNP() llamando a fread()");
        return(0);
    }
    // verifica fileID y versi�n de header NNP.
    if ((headerNNP.fileID[0]!='N')||(headerNNP.fileID[1]!='N')||(headerNNP.fileID[2]!='P')||(headerNNP.version!=1)||(headerNNP.numGenomas>32000))
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 344 en cargarDesdeNNP() error en encabezado NNPv1 %c %c %c %d ,%d",headerNNP.fileID[0],headerNNP.fileID[1],headerNNP.fileID[2],headerNNP.version,headerNNP.numGenomas);
        return(0);
    }
    // cuenta i=0 hasta numGenomas
    if (cantidad>conf->numEspecies)
        cantidad=conf->numEspecies;
    if (cantidad>headerNNP.numGenomas)
        cantidad=headerNNP.numGenomas;
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>NNP%d=",cantidad);
    for (i=0; i<cantidad; i++)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"%d,",i);
        // Lee SNN
        if (snnDataLoader(conf->representantes[i],i,fileIn,conf)==0)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>Error 343.1 en cargarDesdeNNP() llamando a snnDataLoader()");
            return(0);
        }
        // copia el nuevo genoma sobre todos los de su especie.
        for (k=0; k<conf->sizePob; k++)
        {
            if ((conf->pob[k].especie==i)&&(k!=conf->representantes[i]))
                copiarGenoma(conf->representantes[i],k,conf);
        }

    }
    // cierra el archivo NNP
    fclose(fileIn);
    return(1);
}

// todas las funciones de este archivo deben cumplir con que no modifican a los representantes (NMR)
unsigned cicloPrincipal(TConfig* conf)  // OPTIMIZADA    //TODO verificar por NMR todas las funcs en el ciclo, adicionar par�metros al llamar a Funcion perturbarPeso con pert no uniforme, //TODO Funcion perturbarTh.
{
    // Prerequisisto: Tener una conf->poblaci�n inicial usando primeraGen()
    // Realiza el ciclo principal de NEAT, realiza  evaluaci�n, especiaci�n,seleccion(y cruce) y mutaci�n AN + AC (con probabilidades de entrada)
    // hasta que se cumpla maxconf->iteracion o minFitness se alcance.
    // Par�metros:
    //			probMutAN = entre 0 y 1 prob de mutaci�n AN
    //			probMutAC = entre 0 y 1 prob de mutaci�n AC
    //			maxconf->iteracion = m�ximo n�mero de veces que debe correrse el ciclo
    //			minFitness = entre 0 y 1 m�nimo fitness (promediado durante el arch de entrenamiento) necesario para detener el ciclo
    //			maxMemoriaUsada = Max memoria que se puede utilizar en MBytes (si se supera, sale del ciclo)
    //			fileNameGTDv1 = Path para el archivo de datos de entrenamiento en formato GTDv1.
    //			filenameMejor = Path  para el mejor genoma encontrado al  omento de terminar el ciclo principal.
    //			repTrain = n�mero de repeticiones del archivo de entrenamiento(para cada genoma).
    //			maxBufferSize = tam��o m�ximo del buffer de lectura de archivos de entrada., memoria usada = maxBufferSize*4 bytes
    // retorna -1 si hay error, indexpob de genoma con mayor fitness si se cumple una de las condiciones de parada.
    unsigned mejor=0;
    long unsigned memoriaUsada=0;
    unsigned i=0;
    float fitMejor=0;
    int tmpUsarD=conf->usarDistrib;
    // imprime el mejor genoma o de la poblaci�n
    // regla de los \n : uno controla como termina, no como empieza. (solo se puede o no colocar \n al final de cada linea, nunca al principio)
    //						para los errores se hace lo contrario: ellos controlan como empiezan, pero se supone que todo lo que se imprima despu�s
    //						ser�n m�s errores, por tanto, se debe colocar el \n al principio. y un \n al final del error de cicloprincipal./	fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>unsigned  introduciendo GenomaPerfecto como elemento 15");
    /*if (genomaPerfecto(15,conf)==0){ //necesaria esta evaluaci�n antes de primera especiaci�n?
    	fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 59 en funcion cicloPrincipal() llamando a evaluarPob(1,%u,%u,%u)",conf->maxBufferSize,conf->numEntradas,conf->numSalidas);
    	return(UINT_MAX);
    }*/
    // coloca iteracion en 0;
    conf->iteracion=0;

    calcularD(conf);
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>evaluando genoma inicial \n");
    // coloca el n�mero de especies en 1 ya que genoma inicial le asigna la especie 0 en primeraGen();
    conf->numEspecies=1;
    // realiza primera evaluaci�n //TODO: probar si funciona mejor inicializando (primer par�metro)
    if (evaluarPob(1,1,conf->maxBufferSize,conf->fileNameGTDv1, conf)==0)  //necesaria esta evaluaci�n antes de primera especiaci�n?
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 59 en funcion cicloPrincipal() llamando a evaluarPob(1,%u,%u,%u)",conf->maxBufferSize,conf->numEntradas,conf->numSalidas);
        return(UINT_MAX);
    }
    // llama a calcularLimTh de gen.h para ccontrolar el max y min Th de las neuronas de toda la pob
    calcularLimTh(conf);
    // especiaci�n (necesita haberse avaluado al menos una vez el genoma.)
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>Primera especiaci�n.\n");
    if (especiacion(conf)==0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 60 en funcion cicloPrincipal() llamando a especiacion(1,%u,%u,%u)",conf->maxBufferSize,conf->numEntradas,conf->numSalidas);
        return(UINT_MAX);
    }


    // calcula el mejor y su fitness
    if ((mejor=buscarMejorFitness(conf))==UINT_MAX)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 59_144 en cicloPrincipal() llamando a buscarMejorFitness()");
        return(UINT_MAX);
    }
    fitMejor=conf->pob[mejor].fitness;
    // mientras no se cumplan las tres condiciones de parada
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>Comenzando iteraci�n,\n");
    while((conf->pob[mejor].fitness<conf->minFitness)&&(conf->iteracion<conf->maxIteraciones)&&(memoriaUsada<conf->maxMemoriaUsada))
    {
		printf("\nWorkz!!");
        // TODO: en ninguna etapa se deben modificar los representantes excepto en la de perturbar peso (VERIFICAR)
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>Iteraci�n %4u =",conf->iteracion);
        // Imprime los mejores de cada especie
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile," M(%2u) =%3u(%9.9f)e:%2u ",conf->pob[mejor].especie,mejor,fitMejor,conf->numEspecies);
        // TODO: si hay m�s de 3 especies (debido a bug), se realiza sc, sinn� no!!!
        // realiza selecci�n y cruces Funcion seleccionCrossover() pag 54 - 55 phd y otras fuentes,maneja conservaci�n, eliminaci�n y vector de

//TODO: dividir sc en selecci�n y reproducci�n.
        if (conf->numEspecies>2)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"sc");
            if (seleccionCrossover(conf)==0)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 61 en funcion cicloPrincipal() llamando a  seleccionCrossover");
                return(UINT_MAX);
            }
        }
        // realiza mutaciones AN y AC		VERIFICADA
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"m");
        for (i=0; i<conf->sizePob; i++)
        {
            // realiza mutaci�n AN con rconf->pobabilidad probMutAN
            // probabilidad de mutaci�n
            if (((float)randL(conf))<=conf->probMutAN)
            {
                // if (conf->representantes[conf->pob[i].especie]!=i){//Muta si no es el representante
                if (mutarAN(i,conf)==0)
                {
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>\nError 62 en funcion cicloPrincipal() llamando a mutarAN(%u)",i);
                    return(UINT_MAX);
                }
            }
            //}
            // realiza mutaci�n AC con probabilidad porbMutAC
            // probabilidad de mutaci�n
            if (((float)randL(conf))<=conf->probMutAC)
            {
                // if (conf->representantes[conf->pob[i].especie]!=i)//Muta si no es el representante
                if (mutarAC(i,conf->maxIntentosMutarAC,conf)==0)
                {
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>\nError 63 en funcion cicloPrincipal() llamando a mutarAC(%u)",i);
                    return(UINT_MAX);
                }
            }
        }
        // realiza perturbaci�n de pesos 		VERIFICADA
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"w");
// Hace las perturbaciones de peso en paralelo
//#pragma omp parallel for
        for (i=0; i<conf->sizePob; i++)
        {
            // if (conf->representantes[conf->pob[i].especie]!=i){
            // perturbarPesoYth(i,conf->porcentMutPeso,conf->porcentMutNTh,conf->probMutTh,conf->probMutPeso, conf);
            perturbarPeso(i,conf->porcentMutPeso,conf->probMutPeso,0,0,conf->tipoPert,0, conf);
            //}
        }
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"N");
        // introduce poblaci�n desde archivo NNP en diferentes especies
        // FALTA: poner como param cargar al inicio o cargar al final
        // if ((conf->cargarNNP>0)&&(conf->numEspecies==conf->spEspecies))
        if ((conf->cargarNNP>0)&&(conf->numEspecies==conf->cargarNNP))
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nCargando Genomas desde NNP\n");
            if (cargarDesdeNNP(conf->fileNameNNPv1Load,0,conf->cargarNNP,conf)==0)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 355 en cicloPrincipal() llamando a cargarDesdeNNP()");
                return(UINT_MAX);
            }
            conf->cargarNNP=0;
        }
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"D");
        // implementaci�n de procesamiento distribuido(justo antes de la eval)
        if ((conf->cargarDistrib>0)&&(conf->numEspecies==conf->spEspecies)&&(tmpUsarD>0))
        {
            // TODO: FALTA: param para seleccionar si descargar antes  o despu�s.
            // realiza el procesamiento distribuido en los archivo actuales
            if (distribProc(conf)==0)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 355.1 en cicloPrincipal() llamando a distribProc()");
                return(UINT_MAX);
            }
            tmpUsarD--;
        }

/*        // deja descargando en paralelo archivos de procesamiento distrib mientras hace evaluaci�n
        if ((tmpUsarD==0)&&(conf->cargarDistrib>0))
        {
            system(conf->distriCmd);
            tmpUsarD=conf->usarDistrib;
        }
*/
        // evaluaci�n
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"e");
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"C");
        if (evaluarPob(1,0,conf->maxBufferSize,conf->fileNameGTDv1, conf)==0)  //// TODO EL ANTERIORMANTE DICHO PARAMETRO para conf->primero para las evaluaciones y hacer inicializaciones a 0 de los valores cuando se crean las neuronas.
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 59 en funcion cicloPrincipal() llamando a evaluarPob(1,%u,%u,%u)",conf->maxBufferSize,conf->numEntradas,conf->numSalidas);
            return(UINT_MAX);
        }
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"s");
        // especiaci�n (necesita haberse evaluado al menos una vez el genoma.)
        if (especiacion(conf)==0)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 60 en funcion cicloPrincipal() llamando a especiacion(1,%u,%u,%u)",conf->maxBufferSize,conf->numEntradas,conf->numSalidas);
            return(UINT_MAX);
        }
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"Sv");
        // guarda el mejor Genoma
        guardarGenomaSNN(conf->sizePob+(conf->pob[buscarMejorFitness(conf)].especie),conf->fileNameSNNv1,conf);
        // guarda representantes en NNP
        guardarRepresentantesNNP(conf->fileNameNNPv1, conf);
        // si cargarNNP==0 guarda en filenameNNPv1Load (para poder continuar de nuevo si se interrumpe)
        if (conf->cargarNNP==0)
            guardarRepresentantesNNP(conf->fileNameNNPv1Load, conf);
        guardarRepresentantesNNP(conf->fileNameNNPv1, conf);
//Todo: quitar cuando se ahay solucionado problama de delay en especiaci�n
        if (verificarListasInnov(conf)==0)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>Error 60.1 en cicloPrincipal llamando a verificarListaInnov()\n");
            return(UINT_MAX);
        }
        conf->iteracion++; // incrementa conf->iteracion
        //realiza competencia (si el porcentCompetencia es negativo, no hace nada)
        if (conf->porcentCompetencia>0)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"C");
            i=competencia(conf);
        }
        // verifica  y  busca el mejor y su fitness
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"v\n");
        // calcula m�ximo fitness entre los conf->representantes de la conf->pob obtiene el indexpob en mejor
        if ((mejor=buscarMejorFitness(conf))==UINT_MAX)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 59_144 en cicloPrincipal() llamando a buscarMejorFitness()");
            return(UINT_MAX);
        }
        fitMejor=conf->pob[mejor].fitness;
        memoriaUsada=calcularMemoriaUsada(conf->sizePob,conf)/(1048576); // calcula memoria usada en MB
        /*for (j=0;j<conf->numEspecies;j++)
        	fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>E%u[%u]=%3d(%2.2f), ",j,contarGPEsp(j, conf),conf->representantes[j],conf->pob[conf->representantes[j]].fitness);
        //Se llama a calcularLimTh(conf); de gen.h para ccontrolar el max y min Th de las neuronas de toda la pob
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\n");*/
        calcularLimTh(conf);
//while((conf->pob[mejor].fitness<conf->minFitness)&&(conf->iteracion<conf->maxconf->iteracion)&&(memoriaUsada<conf->maxMemoriaUsada)){
        if (conf->pob[mejor].fitness>=conf->minFitness)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>Salida de ciclo principal debido a MinFitness alcanzado.\n");
        }
        if (conf->iteracion>=conf->maxIteraciones)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>Salida de ciclo principal debido a m�ximo n�mero de conf->iteracion.\n");
        }
        if (memoriaUsada>=conf->maxMemoriaUsada)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>Salida de ciclo principal debido a m�xima memoria usada=%lu\n",memoriaUsada);
        }

    }

    return(mejor);
}

//Nueva evaluarPob, versi�n antigua al final.
unsigned evaluarPob(unsigned inicializar,unsigned primero,unsigned maxBufferSize, char *fileNameGTDv1, TConfig* conf)  // OPTIMIZADA, NMR
// lee los archivos SNN y GTD en memoria dentro de los buffers de entrada de la funci�n evalPob
// , luego se llama a la funci�n que eval�a toda la poblaci�n  y finalmente se actualiza el
// fitness de cada genoma de la pob desde el vector de retorno de evalPob();.
// evalPob(int numGenomas, int numDatos, hdrGTDv1 headerGTD, hdrSNNv1* headerSNN, float** dataGTD, tConexDataF** listaConexData, float* fitness)
// retora 0 si error,1 si OK
{
    int numGenomas=conf->sizePob;
    int i,leidos;
    FILE* fileIn;
    // si primero==1
    if (primero==1)
    {
        // abre archivo GTDv1
        fileIn=fopen(conf->fileNameGTDv1,"rb");
        if (!fileIn)
        {
            printf("\nError 57.1 en evaluarPob() llamando a fopen()");
            return(0);
        }
        // lee headerGTD
        leidos=fread(&conf->headerGTD,sizeof(hdrGTDv1),1,fileIn);
        // busca el filesize
        if (leidos<1)
        {
            printf("\nError 58 en evaluarPob() llamando a fread()");
            return(0);
        }
        // averigua el fileSize y calcula conf->numDatos
        fseek(fileIn, 0, SEEK_END);
        conf->numDatos = (ftell(fileIn) - sizeof(hdrGTDv1))/(conf->headerGTD.tamRegistros*(conf->headerGTD.numEntradas+conf->headerGTD.numSalidas));
        //printf("\nFileSize=%d",ftell(fileIn));
        rewind(fileIn);
        leidos=fread(&conf->headerGTD,sizeof(hdrGTDv1),1,fileIn); // para ubicarlo en posici�n de lectura de datos
        // reserva memoria en headerSNN[numGenomas]
        conf->headerSNN=(hdrSNNv1*)malloc(numGenomas*sizeof(hdrSNNv1));
        if (!conf->headerSNN)
        {
            printf("\nError 60 en evaluarPob() llamando a malloc()");
            return(0);
        }
        //reserva memoria para el vector de salidas calculadas para usar en c�lculo de fitness en evalGenom()
        conf->valoresC=(float*)malloc(conf->numDatos*sizeof(hdrSNNv1));
        if (!conf->valoresC)
        {
            printf("\nError 60 en evaluarPob() llamando a malloc()");
            return(0);
        }
        // reserva memoria para conf->ordenEval[i<numGenomas][numNodos[i]]
        conf->ordenEval=(int**)malloc(numGenomas*sizeof(int*));
        if (!conf->ordenEval)
        {
            printf("\nError 60.5 en evaluarPob() llamando a malloc()");
            return(0);
        }
        // reserva memoria para conf->tamOrdenEval[numGenomas]
        conf->tamOrdenEval=(int*)malloc(numGenomas*sizeof(int));
        if (!conf->tamOrdenEval)
        {
            printf("\nError 60.6 en evaluarPob() llamando a malloc()");
            return(0);
        }
        // tamListaConexPost[numGenomas]
        conf->tamListaConexPost=(int*)malloc(numGenomas*sizeof(int));
        if (!conf->tamListaConexPost)
        {
            printf("\nError 60.7 en evaluarPob() llamando a malloc()");
            return(0);
        }
        for (i=0;i<numGenomas;i++)
        {
            conf->ordenEval[i]=(int*)malloc(conf->maxNodos*sizeof(int));
            if (!conf->ordenEval[i])
            {
                printf("\nError 60.7 en evaluarPob() llamando a malloc()");
                return(0);
            }
        }
        // si no son datos tipo float
        if (conf->headerGTD.tamRegistros==4)
        {
            // reserva memoria en fitness[numGenomas]
            conf->fitness=(float*)malloc(numGenomas*sizeof(float));
            if (!conf->fitness)
            {
                printf("\nError 61 en evaluarPob() llamando a malloc()");
                return(0);
            }
            // reserva memoria en listaConexData[i<numGenomas][numConex[i]]
            conf->listaConexData=(tConexDataF**)malloc(numGenomas*sizeof(tConexDataD*));
            if (!conf->listaConexData)
            {
                printf("\nError 62 en evaluarPob() llamando a malloc()");
                return(0);
            }
            for (i=0;i<numGenomas;i++)
            {
                conf->listaConexData[i]=(tConexDataF*)malloc(conf->maxConex*sizeof(tConexDataF));
                if (!conf->listaConexData[i])
                {
                    printf("\nError 63.11 en evaluarPob() llamando a malloc()");
                    exit(0);
                }
            }
            // reserva memoria en conf->dataGTD[conf->numDatos][numEntradas+numSalidas]
            // para todo el buffer, ojo! obtiene toda la memoria para el index 0 para que sean consecutivos
            // durante la lectura
            conf->dataGTDf=(float**)malloc(conf->numDatos*sizeof(float*));
printf("NumDatos=%d\n",conf->numDatos);
            if (!conf->dataGTDf)
            {
                printf("\nError 64 en evaluarPob() llamando a malloc()");
                return(0);
            }
            conf->dataGTDf[0]=(float*)malloc((conf->headerGTD.numEntradas+conf->headerGTD.numSalidas)*conf->numDatos*sizeof(float));
            if (!conf->dataGTDf[0])
            {
                printf("\nError 65 en evaluarPob() llamando a malloc((%d+%d)*%d)",conf->headerGTD.numEntradas,conf->headerGTD.numSalidas,conf->numDatos*sizeof(float));
                return(0);
            }
            // coloca las filas de la matriz del buffer GTD apuntando al inicio de cada grupo de datos.
            for (i=1;i<conf->numDatos;i++)
            {
                //TODO: OJO!!!! Verificar si es sizeof(float)*(nIn+nOut)*i o (nIn+nOut)*i
                conf->dataGTDf[i]=conf->dataGTDf[0]+((conf->headerGTD.numEntradas+conf->headerGTD.numSalidas)*i);
            }
            // lee de GTD (numEntradas+numSalidas)*sizeof(float) en conf->dataGTD[i];
            leidos=fread(conf->dataGTDf[0],sizeof(float)*(conf->headerGTD.numEntradas+conf->headerGTD.numSalidas),conf->numDatos,fileIn);
            // verifica lectura.
            if (leidos<conf->numDatos)
            {
                printf("\nError 66 en evaluarPob() llamando a fread()");
                return(0);
            }
           // for (i=0;i<conf->numDatos;i++)
            //{
             //   printf("E0=%3.3f E1=%3.3f E2=%3.3f E3=%3.3f S0=%3.3f\n",conf->dataGTDf[i][0],conf->dataGTDf[i][1],conf->dataGTDf[i][2],conf->dataGTDf[i][4],conf->dataGTDf[i][5]);
            //}
        }
        // CUDA: reserva memoria en device CUDA para los datos de entrenamiento.
        cudaMalloc(&(conf->cu_dataGTDf),(conf->headerGTD.numEntradas+conf->headerGTD.numSalidas)*conf->numDatos*sizeof(float));
        // coloca las filas de la matriz del buffer GTD apuntando al inicio de cada grupo de datos.
/*        for (i=1;i<conf->numDatos;i++)
        {
//TODO: OJO: SE DEBE HACER EN EL KERNEL.
            //TODO: OJO!!!! Verificar si es sizeof(float)*(nIn+nOut)*i o (nIn+nOut)*i
            conf->dataGTDf[i]=conf->dataGTDf[0]+((conf->headerGTD.numEntradas+conf->headerGTD.numSalidas)*i);
        }
*/
        // CUDA: reserva memoria en device para listaConexData[sizePob]
        cudaMalloc(&(conf->cu_conexIn),conf->sizePob*sizeof(int*));
        // CUDA: reserva memoria en device para listaConexData[sizePob]
        cudaMalloc(&(conf->cu_conexOut),conf->sizePob*sizeof(int*));
        // CUDA: reserva memoria en device para listaConexData[sizePob]
        cudaMalloc(&(conf->cu_conexPeso),conf->sizePob*sizeof(float*));
        // CUDA: reserva memoria en device para tamListaConexData[sizePob]
        cudaMalloc(&(conf->cu_tamListaConexPost),conf->sizePob*sizeof(int));
        // CUDA: reserva memoria en device para fitness[sizePob]
        cudaMalloc(&conf->cu_fitness,conf->sizePob*sizeof(float));
        // CUDA: reserva memoria en device para valorC[sizePob] usado para calculo de fitness con correlaci�n.
        cudaMalloc(&(conf->cu_valorC),conf->sizePob*sizeof(float*));
        for (i=0;i<conf->sizePob;i++)
            cudaMalloc(&(conf->cu_valorC[i]),conf->numDatos*sizeof(float));
        // CUDA: reserva memoria en device para valorC[sizePob] usado para calculo de fitness con correlaci�n.
        cudaMalloc(&(conf->cu_valorTr),conf->sizePob*sizeof(float*));
        for (i=0;i<conf->sizePob;i++)
            cudaMalloc(&(conf->cu_valorTr[i]),conf->numDatos*sizeof(float));
        // CUDA: transfiere los datos de entrenamiento.dataGTDf.
        printf("Copiando datos de entrenamiento a dispositivo CUDA.\n");
        cudaMemcpy(conf->cu_dataGTDf, conf->dataGTDf[0], sizeof(float)*(conf->headerGTD.numEntradas+conf->headerGTD.numSalidas)*conf->numDatos, cudaMemcpyHostToDevice);
    }
    // si es o no primero
    // para i=0;i<numGenomas;i++
    for (i=0;i<numGenomas;i++)
    {
        // genera headerSNN[i]
        conf->headerSNN[i].actThreshold=conf->Fthreshold;
        conf->headerSNN[i].numEntradas=conf->numEntradas;
        conf->headerSNN[i].numSalidas=conf->numSalidas;
        conf->headerSNN[i].numBias=conf->numBias;
        conf->headerSNN[i].numHiddens=conf->pob[i].totalNodos-(conf->numEntradas+conf->numSalidas+conf->numBias);
        conf->headerSNN[i].numConex=conf->pob[i].totalConexiones;
        conf->headerSNN[i].sigmaFactor=conf->A;
        conf->headerSNN[i].tamRegistros=(conf->useFloat==1?4:8);
        conf->headerSNN[i].usarSigned=conf->usarSigned;
    }
    // llama a evalPob() (PARALELIZABLE(secci�n dentro de ella))
    evalPob(numGenomas, conf->numDatos, conf->headerGTD, conf->headerSNN, conf->dataGTDf, conf->listaConexData, conf->fitness, conf->ordenEval, conf->tamListaConexPost , conf);
    // coloca el fitness en cada genoma, para i=0;i<numGenomas;i++
    for (i=0;i<numGenomas;i++)
         conf->pob[i].fitness=conf->fitness[i];
       // printf("F%3.3f",conf->pob[i].fitness=conf->fitness[i]);
	
    return(1);
}

/*
unsigned evaluarPob(unsigned inicializar,unsigned primero,unsigned maxBufferSize, char *fileNameGTDv1, TConfig* conf)  // OPTIMIZADA, NMR
{
// eval�a toda la poblaci�n y deja el valor post fSigma en cada nodo.
// y calcula el fitness basado en el que se va acumulando con cada evaluaci�n de cada genoma.
// Los archivos de entrada y salida deben estar previamente abiertos para lectura binaria br
// retorna 0 si hay error, 1 si ok.
// //TODO: URGENTE PARA FX par�metro repeticiones para pasar los archivos de entrada repetidas veces por las redes neuronales al realizar las evaluaciones.
    unsigned i;
    unsigned j;
    int leidos=0;
    double pasadas=0; // usado para normalizar el fitness
    FILE* fileIn;
    //coloca los valores de cada neurona en 0 para comenzar la evaluaci�n. si incializar=1;
    //tambi�n Actualiza los punteros a conexHijo de cada nodo necesarios para evaluarGenoma
    if (primero==1)
    {
        //reserva memoria para sistema de ordenEval[][]
        // reserva memoria para conf->ordenEval[i<numGenomas][numNodos[i]]
        conf->ordenEval=(int**)malloc(conf->realSizePob*sizeof(int*));
        if (!conf->ordenEval)
        {
            printf("\nError 61.5 en .evaluarPob() llamando a malloc()");
            return(0);
        }
        // reserva memoria para m�ximo conf->maxNodos nodos por genoma
        for (i=0;i<conf->realSizePob;i++)
        {
            conf->ordenEval[i]=(int*)malloc(conf->maxNodos*sizeof(int));
            if (!conf->ordenEval[i])
            {
                printf("\nError 61.6 en .evaluarPob() llamando a malloc()");
                return(0);
            }
        }
        // reserva memoria para conf->contO[numGenomas] usado como contador para dentro de la funci�n recursiva
        conf->contO=(int*)malloc(conf->realSizePob*sizeof(int));
        if (!conf->contO)
        {
            printf("\nError 61.8 en .evaluarPob() llamando a malloc()");
            return(0);
        }
        // abre archivo GTDv1
        fileIn=fopen(conf->fileNameGTDv1,"rb");
        if (!fileIn)
        {
            printf("\nError 57.1 en evaluarPob() llamando a fopen()");
            return(0);
        }
        // lee headerGTD
        leidos=fread(&conf->headerGTD,sizeof(hdrGTDv1),1,fileIn);
        // busca el filesize
        if (leidos<1)
        {
            printf("\nError 58 en evaluarPob() llamando a fread()");
            return(0);
        }
        // averigua el fileSize y calcula conf->numDatos
        fseek(fileIn, 0, SEEK_END);
        conf->numDatos = (ftell(fileIn) - sizeof(hdrGTDv1))/(conf->headerGTD.tamRegistros*(conf->headerGTD.numEntradas+conf->headerGTD.numSalidas));
        if (conf->numDatos<1)
        {
            printf("\nError 65 en evaluarPob()llamando a ftell()");
            return(0);
        }
        //printf("\nFileSize=%d",ftell(fileIn));
        rewind(fileIn);
        // ubica en posici�n de lectura de datos despu�s del header.
        leidos=fread(&conf->headerGTD,sizeof(hdrGTDv1),1,fileIn);
        // reserva memoria en conf->dataGTD[conf->numDatos][numEntradas+numSalidas]
        // para todo el buffer, ojo! obtiene toda la memoria para el index 0 para que sean consecutivos
        // durante la lectura
        conf->dataGTDf=(float**)malloc(conf->numDatos*sizeof(float*));
        if (!conf->dataGTDf)
        {
            printf("\nError 64 en evaluarPob() llamando a malloc()");
            return(0);
        }
        conf->dataGTDf[0]=(float*)malloc((conf->headerGTD.numEntradas+conf->headerGTD.numSalidas)*conf->numDatos*sizeof(float));
        if (!conf->dataGTDf[0])
        {
            printf("\nError 65 en evaluarPob() llamando a malloc((%d+%d)*%d)",conf->headerGTD.numEntradas,conf->headerGTD.numSalidas,conf->numDatos*sizeof(float));
            return(0);
        }
        // coloca las filas de la matriz del buffer GTD apuntando al inicio de cada grupo de datos.
        for (i=1;i<conf->numDatos;i++)
        {
            //TODO: OJO!!!! Verificar si es sizeof(float)*(nIn+nOut)*i o (nIn+nOut)*i
            conf->dataGTDf[i]=conf->dataGTDf[0]+((conf->headerGTD.numEntradas+conf->headerGTD.numSalidas)*i);
        }
        // lee de GTD (numEntradas+numSalidas)*sizeof(float) en conf->dataGTD[i];
        leidos=fread(conf->dataGTDf[0],sizeof(float)*(conf->headerGTD.numEntradas+conf->headerGTD.numSalidas),conf->numDatos,fileIn);
        // verifica lectura.
        if (leidos<conf->numDatos)
        {
            printf("\nError 66 en evaluarPob() llamando a fread()");
            return(0);
        }
        fclose(fileIn);
    }
    if (inicializar==1)
    {
        for (i=0; i<conf->sizePob; i++)
        {
            //Actualiza los punteros a conecHijo de los nodos del genoma.
            if (actualizarPNodos(i,conf)==0)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 1150 en fnci�n evaluarPob() llamando a actualizarPNodos()\n");
                return(0);
            }
            for (j=0; j<conf->pob[i].totalNodos; j++)
            {
                if (conf->pob[i].nodo[j].nodeFunction!=3)
                    conf->pob[i].nodo[j].valor=0;
            }
        }
    }
    //Inicializa el acumulador de fitness
    for (i=0; i<conf->sizePob; i++)
    {
        conf->pob[i].fitness=0;
    }
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"0");
    //hace repeticiones deben ser 0 o m�s
    for (i=0; i<conf->numDatos; i++)
    {
        pasadas = pasadas+1;
        //eval�a entradas y salidas para todos los genomas
        //TODO: ES m�s r�pido si a evaluargenoma se le pasa todo el buffer en lugar de dato a dato
        //con pragme 0,1,2=1:38s,  SIN=48s
        //#pragma omp parallel for
        for (j=0; j<conf->sizePob; j++)
        {
            if(evaluarGenoma(j,&(conf->pob[j]),primero, &(conf->dataGTDf[i][0]) ,&(conf->dataGTDf[i][conf->headerGTD.numEntradas]),conf)==0) ////TODO: HACER PARAMETRO GLOBAL conf->PRIMERO PARA DIFERENTES APLICACIONES
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 59 en funcion evaluarPob() llamando a evaluarGenoma(genoma=%u,i=%u )\n",j,i);
            }
        }
    }
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"1");
    // calcula los fitness a partir del acumulado de cada genoma
    if (pasadas>0)
        for (i = 0; i < conf->sizePob; i++)
        {
            //conf->pob[i].fitness /= (conf->tSigma > 1000 ? 2*(conf->numSalidas*pasadas) : conf->numSalidas*pasadas);
            conf->pob[i].fitness /= (conf->numSalidas*pasadas);
            conf->pob[i].fitness = 1.0 - conf->pob[i].fitness;
        }
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"2");
    // TODO: colocar como par�metro si copiar en cada evaluaci�n el mejor genoma de una especie en luagar del peor de la especie para que pueda mutar pesos,etc...
    // cierra archivos de entrada y salida
    fclose(conf->fIn);
    // libera los punteros usados por los punteros a las conexhijos usados durante eval (deben ser rearmados nuevamente con actualizarPNodos).
    for (i=0; i<conf->sizePob; i++)
    {
        for (j=0; j<conf->pob[i].totalNodos; j++)
        {
            if ((conf->pob[i].nodo[j].conexHijo!=NULL)&&(conf->pob[i].nodo[j].contHijos>0)&&(conf->pob[i].nodo[j].nodeFunction!=3))
            {
                free(conf->pob[i].nodo[j].conexHijo);
            }
        }
    }
    //libera memoria de las listas inputs y outputs
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"3");
    ////TODO: hacer Funcion bufferize para leer las entradas y salidas en arreglos globales para que se lean una sola vez del disco durante todo el ciclo ppal
    return(1);
}
*/
unsigned seleccionCrossover(TConfig* conf)  // OPTIMIZADA, NMR (No modifica representantes) TODO: falta mirar si formula de numHijos funcionam bi�n
{
    // realiza selecci�n y cruce en toda la poblaci�n
    // Inicializa el arreglo conf->numGenomasPorEspecie (conf->fInal de pag 394(426) de AI game prog) y lo modifica para reducir a la mitad el n�mero de especies
    // en conservaci�n. Y reparte el resto entre las especies que n� est�n en conservaci�n.
    // Requiere haber evaluado fitness y haber realizado especiacion.
    // Retorna 0 si hay error, 1 si ok
    // Par�metros:	porcentRedConserv = (entre 0 y 1) porcentaje de genomas que se quitan al n�mero calclulado para especies en conserv. recom=0.3
    // porcentElim = (entre 0 y 1) porcentaje de genomas que se eliminan en cada generaci�n, el restante porcentaje se reproduce por cruce.
    // unsigned  intentosPareja= n�mero de unsigned  intentos para buscar pareja aleatoriamente entre los padres antes de "matrimonio forzoso :) "
    // super = float entre 0 y 1, Probabilidad de heredar los exess y disjounsigned  ints del menos apto (aparte de los que se heredan normalmente del m�s apto)
    // promendiarProb = float entre 0 y 1 = probabilidad de que en caso de matching, se promedien los pesos en lugar de
    // seleccionarlos aleatoriamente entre los padres.
    unsigned i;
    unsigned j;
    unsigned dif=0;
    unsigned k=0;
    float avgTotal=0; // para almacenar la suma de los promedios de fitness de todas las especies
    unsigned verificacion =0; // usado para verificar si la sumatoria de conf->numGenomasPorEspecie==conf->sizePob
    unsigned variacion=0; // IGUAL QUE VERIFICACI�N
    unsigned tmp;
    unsigned huboSwap=1;// para reliazr ordenamiento de listaOrden
    unsigned sentido=1;// usado para ordenar listaOrden
    unsigned x;
    unsigned Ge=30005; // n�mero de genomas por especie (total de hijos de toda la especie)
    unsigned Gp=7;// n�mero de genomas padre.
    float m;
    unsigned y;
    unsigned acum=0;
    float correc=0.0001;
    unsigned sum=0;
    unsigned huboReprod=1; // usado para verificar reproducci�n de todos los padres
    unsigned posLMadre=0; // posici�n en la lista Ordenada de la madre
    unsigned posHijo=0; // indexpob del hijo
    unsigned especie1=0;
    unsigned especieCercana=0;

    // reserva memoria para *conf->fitnessAvgPorEspecie
    if ((conf->fitnessAvgPorEspecie=(float*)calloc(conf->spEspecies,(unsigned  int)sizeof(float)))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 51 en funcion seleccionCruce() llamando a calloc(%u,%u)",conf->numEspecies,(unsigned  int)sizeof(float));
        return(0);
    }
    // inicializa  a a0
    for (i=0; i<conf->numEspecies; i++)
    {
        conf->fitnessAvgPorEspecie[i]=0;
    }
    // reserva memoria para conf->actNumGenomasPorEspecie de tama�o conf->numEspecies
    if ((conf->actNumGenomasPorEspecie=(unsigned  int*)calloc(conf->spEspecies,(unsigned  int)sizeof(unsigned  int)))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 52 en funcion seleccionCruce() llamando a calloc(%u,%u)",conf->numEspecies,(unsigned  int)sizeof(unsigned  int));
        return(0);
    }
    // inicializa  a 0
    for (i=0; i<conf->numEspecies; i++)
    {
        conf->actNumGenomasPorEspecie[i]=0;
    }
    // reserva memoria para la matriz de indexpob por especie[f][c]
    // memoria para filas:
    if((conf->listaOrdenFitness=(unsigned  int**)calloc(conf->spEspecies,(unsigned  int)sizeof(unsigned  int*)))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 53 en funcion seleccionCruce() llamando a malloc(%u)",conf->numEspecies*(unsigned int)sizeof(unsigned  int*));
        return(0);
    }
    // memoria para columnas y coloca el n�mero de hijos de todos los genomas en 0
    // coloca el n�mero de hijos para todos los genomas de la conf->poblaci�n en 0 antes de calcular su n�mero de hijos.
    for (i=0; i<conf->numEspecies; i++)
    {
        if((conf->listaOrdenFitness[i]=(unsigned  int*)calloc(1,sizeof(unsigned  int)*conf->sizePob))==NULL)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 54 en funcion seleccionCruce() en asignacion conf->listaOrdenFitness[%u]=(unsigned  int*)malloc(%u)",i,conf->sizePob*(unsigned int)sizeof(unsigned  int));
            return(0);
        }
    }
    // inicializa actNumGenomas por epecie a 0
    // se calcula con formula el conf->numGenomasPorEspecie para cada especie, adem�s se llena la matriz de orden y conf->actNumGenomasPorEspecie
    for (i=0; i<conf->sizePob; i++) //se barre conf->pob
    {
        conf->listaOrdenFitness[conf->pob[i].especie][conf->actNumGenomasPorEspecie[conf->pob[i].especie]]=i; //lena el indexpob correspondiente en la lista de orden
        conf->fitnessAvgPorEspecie[conf->pob[i].especie]+=conf->pob[i].fitness; //acumula fitness por especie
        conf->actNumGenomasPorEspecie[conf->pob[i].especie]++; //se actualiza conf->actNumGenomasPorEspecie
    }
    // para el calculo se divide conf->fitnessAvgPorEspecie entre conf->actNumGenomasPorEspecie correspondiente para obtener el conf->fitnessAvgPorEspecie
    for (i=0; i<conf->numEspecies; i++) //se barren las especies
    {
        conf->fitnessAvgPorEspecie[i]/=(float)conf->actNumGenomasPorEspecie[i];
        avgTotal+=conf->fitnessAvgPorEspecie[i];
    }
    // verifica si alguno de los genomas tiene el m�nimo de fitness requerido para tener minPorcentGenomasPorEspecie y si no lo asigna y recalcula el avgTotal
    for (i=0; i<conf->numEspecies; i++)
    {
        if ((conf->fitnessAvgPorEspecie[i]/avgTotal)<conf->minPorcentGenomasPorEspecie)
        {
            conf->fitnessAvgPorEspecie[i]=conf->minPorcentGenomasPorEspecie;
        }
    }

    /** PROBANDO!!!! CAMBIADO El fitness sharing para distribuci�n de poblaci�n, ahora depende del fitness solo del rep no del prom de todos
        // recalcula el n�mero de genomas por especie.
        avgTotal=0;
        for (i=0; i<conf->numEspecies; i++) //se barren las especies
        {
            avgTotal+=conf->fitnessAvgPorEspecie[i];
        }
            //calcula el fitness para cada especie dependiendo de el fitness de la especie y el fitness total


    //TODO:formula cambiada buscando bug de ciclo infinito         conf->numGenomasPorEspecie[i]=(unsigned  int)floor(((conf->fitnessAvgPorEspecie[i]/avgTotal)*(float)conf->sizePob)+0.5); //calcula el n�mero de genomas para cada especie

        for (i=0; i<conf->numEspecies; i++) //se barren las especies
        {
            conf->numGenomasPorEspecie[i]=(unsigned  int)floor((conf->fitnessAvgPorEspecie[i]/avgTotal)*(float)conf->sizePob)+0.5; //calcula el n�mero de genomas para cada especie
            verificacion+=conf->numGenomasPorEspecie[i];
        }
    */
//INICIO DE PRUEBA
    // recalcula el n�mero de genomas por especie.
    avgTotal=0;
    for (i=0; i<conf->numEspecies; i++) //se barren las especies
    {
        avgTotal+= conf->pob[conf->representantes[i]].fitness;
    }
    //calcula el fitness para cada especie dependiendo de el fitness de la especie y el fitness total


//TODO:formula cambiada buscando bug de ciclo infinito         conf->numGenomasPorEspecie[i]=(unsigned  int)floor(((conf->fitnessAvgPorEspecie[i]/avgTotal)*(float)conf->sizePob)+0.5); //calcula el n�mero de genomas para cada especie

    for (i=0; i<conf->numEspecies; i++) //se barren las especies
//TODO: probando para debugging de sigerr, antes no se le restaba numespecies a sizePob y solo se le sumaba 0.5 en lubgar de 1.5
    {
        conf->numGenomasPorEspecie[i]=(unsigned  int)floor((conf->pob[conf->representantes[i]].fitness/avgTotal)*(float)(conf->sizePob-(3*conf->numEspecies)))+3.5; //3 porque 1 champios,2 distrib, 3 worker//calcula el n�mero de genomas para cada especie
        verificacion+=conf->numGenomasPorEspecie[i];
    }
//FIN DE PRUEBA
    //verifica que la sumatoria de numgenomas por especie coincida con conf->sizePob, si es superior, resta genomas barriendo a cada especie en conf->numGenomasPorEspecie
    i=0;
    if (verificacion>conf->sizePob) //si hay m�s que los debidos
    {
        variacion=verificacion-conf->sizePob;
        if (conf->numGenomasPorEspecie[i]>1)
            while(variacion--)  //Disminuye uniformemente el conf->numGenomasPorEspecie variaci�n veces.
            {
                conf->numGenomasPorEspecie[i++]--;
                if(i==conf->numEspecies)
                    i=0;
            }
    }
    if (verificacion<conf->sizePob) //si hay menos que los debidos
    {
        variacion=conf->sizePob-verificacion;
        while(variacion--)  //Incrementa uniformemente el conf->numGenomasPorEspecie variaci�n veces.
        {
            conf->numGenomasPorEspecie[i++]++;
            if(i==conf->numEspecies)
                i=0;
        }
    }

    // imprime n�mero de genomas por especie
    /*   for (j=0; j<conf->numEspecies; j++)
       {
           fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>E%i=%u,",j,conf->numGenomasPorEspecie[j]);
       } */
    // se ordena la conf->listaOrdenFitness para cada especie con prioridad parael fitness y segunda prio la simplicidad (numConex)
    // poco sofisticado, pero funciona.

    for (i=0; i<conf->numEspecies; i++) //barre las especies
    {
        huboSwap=1;//para que entre la primera vez
        while(huboSwap)
        {
            huboSwap=0;
            // TODO: optimizar con punteros y variables toda esta secci�n
            if(sentido==1) //para ordenamiento descendente desde el primer elemento
            {
                for (j=0; j<(conf->actNumGenomasPorEspecie[i]-1); j++) //Barre los index de la especie
//TODO: quitar cuando corrijamos bug hang por max uint (desbordamiento)
                    if ((conf->actNumGenomasPorEspecie[i]-1)<=conf->sizePob)
                    {
                        if (conf->pob[conf->listaOrdenFitness[i][j]].fitness<conf->pob[conf->listaOrdenFitness[i][j+1]].fitness)
                        {
                            huboSwap=1;
                            swap(&(conf->listaOrdenFitness[i][j]),&(conf->listaOrdenFitness[i][j+1]));
                        }
                        if (conf->pob[conf->listaOrdenFitness[i][j]].fitness==conf->pob[conf->listaOrdenFitness[i][j+1]].fitness)
                        {
                            if(conf->pob[conf->listaOrdenFitness[i][j]].totalConexiones>conf->pob[conf->listaOrdenFitness[i][j+1]].totalConexiones)
                            {
                                huboSwap=1;
                                swap(&(conf->listaOrdenFitness[i][j]),&(conf->listaOrdenFitness[i][j+1]));
                            }
                        }
                    }
            }
            else //para ordenamiento ascendente desde el �ltimo elemento
            {
                for (j=(conf->actNumGenomasPorEspecie[i]-1); j>=1; j--) //Barre los index de la especie
//TODO: quitar cuando corrijamos bug hang por max uint (desbordamiento)
                    if (((conf->actNumGenomasPorEspecie[i]-1)<conf->sizePob)&&(j<=conf->sizePob))
                    {

                        if (conf->pob[conf->listaOrdenFitness[i][j-1]].fitness<conf->pob[conf->listaOrdenFitness[i][j]].fitness)
                        {
                            huboSwap=1;
                            swap(&(conf->listaOrdenFitness[i][j-1]),&(conf->listaOrdenFitness[i][j]));
                        }

                        if (conf->pob[conf->listaOrdenFitness[i][j-1]].fitness==conf->pob[conf->listaOrdenFitness[i][j]].fitness)
                        {
                            if(conf->pob[conf->listaOrdenFitness[i][j-1]].totalConexiones>conf->pob[conf->listaOrdenFitness[i][j]].totalConexiones)
                            {
                                huboSwap=1;
                                swap(&(conf->listaOrdenFitness[i][j-1]),&(conf->listaOrdenFitness[i][j]));
                            }
                        }

                    }

            }
            sentido*=-1; //invierte el sentido de la b�squeda
        }

    }


    //reorganiza listaorden: SOBREESCRIBE genomas faltantes en para que numgenomasPE y actnum sean iguales //TODO: verificar si esto es necesario.
    for (i=0; i<conf->numEspecies; i++)
    {
        if(conf->actNumGenomasPorEspecie[i]<conf->numGenomasPorEspecie[i])
        {
            dif=conf->numGenomasPorEspecie[i]-conf->actNumGenomasPorEspecie[i];
            while(dif)
            {
                for (j=0; j<conf->numEspecies; j++)
                {
                    if(conf->actNumGenomasPorEspecie[j]>conf->numGenomasPorEspecie[j])
                    {
                        dif--;
                        conf->actNumGenomasPorEspecie[j]--;
                        conf->listaOrdenFitness[i][conf->actNumGenomasPorEspecie[i]]=conf->listaOrdenFitness[j][conf->actNumGenomasPorEspecie[j]];
                        copiarGenoma(conf->listaOrdenFitness[i][0],conf->listaOrdenFitness[i][conf->actNumGenomasPorEspecie[i]], conf);
                        conf->actNumGenomasPorEspecie[i]++;
                        if (dif==0)
                            j=conf->numEspecies;
                    }
                }
            }
        }
    }

    //TODO: quitar esta parte SIGUIENTE  cuando funcione bi�n
    //verifica que la suma de numGenomasPorEspecie sea igual a sizePob
    acum=0;
    for (i=0; i<conf->numEspecies; i++)
    {
        acum+=conf->actNumGenomasPorEspecie[i];

    }
    if (acum!=conf->sizePob)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError an SeleccionCrossover() tama�o de la pob (%u) no coincide con suma de numGenomasxespecie(%u).",conf->sizePob,acum);
        return(0);
    }
    // AsignarNumHijos: se coloca descargando del total de hijos por especie calculado, el numHijos para cada genoma
    // Con porcentElim calcula el n�mero de genomas por especie que se deben reproducir y les asigna su respectivo numHijos.
    // corresponde a un segmento de linea con pendiente negativa entre x=0 y x=genomasPadre y con unsigned  integral=conf->numGenomasPorEspecie.
    // se hace barriendo la conf->listaOrdenFitness para cada especie y asignando el numHijos descontando un tanto de totalHijosPorEspecie-1 (debido a que el came�n se conserva y el numhijos se disminuye en 1) y si =0,
    //TODO: comprobar si esta formula funciona y leer y remover si es neceario el comment de arriba
    for (i=0; i<conf->numEspecies; i++) //barre las especies
    {
        Ge=conf->numGenomasPorEspecie[i]-1;// era -1 al finaln�mero de genomas para reemplazar (todos los de la especie menos el champion)
        if (Ge==-1)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 23 en seleccionCrossover(), numGenomasPorEspecie[%d]=%d\n",i,conf->numGenomasPorEspecie[i]);
            return(0);
        }
        sum=0;
        Gp=(unsigned  int)floor((float)Ge*(1-conf->porcentElim)+0.5);//n�mero padres por especie;
        for (x=1; x<=Gp; x++) //calcula el decremento de hijos por genoma.
        {
            sum+=x;
        }
        m=(float)Ge/(float)sum;// decremento de hijos por genoma
        acum=0;// para realizar verificaci�n de n�mero correcto de genomas.
        for (x=Gp; x>0; x--) //calcula el n�mero de hijos por genoma
        {
            y=(unsigned  int)floor(((m*(float)x)+0.5)+correc); //n�mero de hijos para genoma Gp-x de especie i
            conf->pob[conf->listaOrdenFitness[i][Gp-x]].numHijos=y;
            // fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>y(%u) = %u = %3.5f, ",Gp-x,y,(m*(float)x+0.5)+b-0.5+correc);
            acum+=y;
        }
        acum=Ge-acum;//calcula cuantos genomas sobraron o faltaron y los descuenta del champion.
        if(acum!=0) conf->pob[conf->listaOrdenFitness[i][0]].numHijos+=acum;
    }

    //imprimirSeleccion(conf);

// VERIFICAR QUE LOS HIJOS SUMEN FITNESS CORRECTO
    //for (x=0;x<Gp)
    //Se hace reproducci�n :))))  :)
    //para esto se barre la conf->listaOrdenFitness por cada especie y se realiza cruce entre padre y una madre aleatoria  (n unsigned  intentos) que tenga
    //el numHijos>0 si no se encuentra aleatorio, se barre de izq a der en busqueda.
    //Si el numHijos del PADRE es 1, se reemplaza el indexpob del padre por el hijo y se coloca numHijos=-1.
    //sin�, se coloca el resultado reemplazando a un genoma que tenga numHijos=0 EN TODA LA conf->pob (si despu�s del barrido no se encuentra (pej �ltimogenoma))
    //se toma como madre el genoma campe�n (index 0 de listaorden) y se reemplaza el del padre.
    //con cada reproducci�n se disminuye el numHijos del padre. si es =0 se pasa al pr�ximo genoma de listaOrden.
    //Si el genoma es un campe�n de especie, y su numHijos=1 se reproduce, se bua un numHijos=0 en conf->pob y se coloca ahi a su hijo, el numHijos
    //del campe�n se hace -1 para evitar que sea reemplazado por otros cruces y se conserve unsigned  intacto.
    //si no se encuentra madre con numHijos>0, se reproduce con champion y hijo sustituye a padre.
    for (i=0; i<conf->numEspecies; i++) //barre especies
    {
        if (conf->contGeneracSinMejora[i]<conf->maxGenParaNoCruce)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"%d",i);
            Gp = (unsigned  int)floor(((float)conf->numGenomasPorEspecie[i]-1)*(1-conf->porcentElim));//n�mero padres por especie;
            huboReprod=1; //inicializa para cada especie el huboReprod.
            // numHijos=0 para todos los dem�s genomas de listaOrdenFitness a partir de Gp
            for (k=Gp; k<conf->numGenomasPorEspecie[i]; k++)
                conf->pob[conf->listaOrdenFitness[i][k]].numHijos=0;


            while(huboReprod) //si hubo reproducci�n en el ciclo anterior
            {
                huboReprod=0;
                for (j=0; j<Gp; j++) //barre padres
                {
                    if (conf->pob[conf->listaOrdenFitness[i][j]].numHijos>0) //si el padre elegido tiene numHijos>0
                    {
                        tmp=conf->intentosPareja;
                        posLMadre=(unsigned  int)floor((Gp*((float)randL(conf)))+0.5000001);//posMadre=busca aleatoriamente entre 0 y Gp una madre se coloca por si unsigned  intentos=0
                        while (tmp--)  //mientras haya unsigned  intentos disponibles
                        {
                            if (Gp>1) // para evitar ciclo infinito y para evitar que los dos padres sean el mismo genoma
                                while((posLMadre=(unsigned  int)floor((Gp*((float)randL(conf)))+0.5000001))==j);//mientras ((posMadre=aleatoriio(0,Gp))==PosPadre);
                            if (conf->pob[conf->listaOrdenFitness[i][posLMadre]].numHijos>0)//si madre.numHijos>0 tmp=0
                                tmp=0;
                        }
                        if (conf->pob[conf->listaOrdenFitness[i][posLMadre]].numHijos==0) //si madre.numHijos==0
                        {
                            for (k=0; k<Gp; k++) //para k=0;k<Gp;k++ si no encontr� madre aleatoriamente, la busca secuencialmente
                            {
                                if (conf->pob[conf->listaOrdenFitness[i][k]].numHijos>0)//si si posK.numHijos>0
                                    if (j!=k)//si posK!=j,posMadre=k
                                        posLMadre=k;
                            }
                        }
                        if (conf->pob[conf->listaOrdenFitness[i][posLMadre]].numHijos>0) //si madre.numHijos>0
                        {
                            posHijo=UINT_MAX;//posHijo=-1
                            if (conf->pob[conf->listaOrdenFitness[i][j]].numHijos==1) //si padre.numHijos==1
                            {
                                if (j>0) //si posPadre>0 posHijo=listaOrden[i][j]
                                    posHijo = conf->listaOrdenFitness[i][j];
                                else // coloca el n�mero de hijos en -1 para el representante
                                    conf->pob[conf->listaOrdenFitness[i][j]].numHijos=-1;//sino padre.numHijos=-1 porque ee el champion
                            }

                            // adicionado para  competencia: busca posHijo y verifica que el index de un representante no sea escogido como hijo
                            // si posHijo==-1 si a�n no se tiene posici�n para el hijo busca en especie
                            if (posHijo==UINT_MAX)
                            {
                                if ((conf->numGenomasPorEspecie[i]-1)<conf->sizePob)
                                    for (k=(conf->numGenomasPorEspecie[i]-1); k>1; k--)
                                    {
                                        //para no sobreescribir el champion, se hace hasta k>0 no=0
                                        if(conf->pob[conf->listaOrdenFitness[i][k]].numHijos==0)
                                        {
                                            posHijo=conf->listaOrdenFitness[i][k];
                                            k=1; //sale del for
                                        }
                                    }
                            }

                            // si posHijo==-1 si a�n no se tiene posici�n para el hijo, busca en pob
                            if (posHijo==UINT_MAX)
                            {
                                for (k=0; k<conf->sizePob; k++) //busca primero entre los de su especie, luego busca entre toda la conf->poblaci�n un genoma que tenga numHijos=0 y guarda su indexpob en posHijo
                                    if(conf->pob[k].numHijos==0)
                                    {
                                        if(conf->representantes[conf->pob[k].especie]!=k)
                                        {
                                            posHijo=k;
                                            k=conf->sizePob; //sale del for
                                        }
                                    }
                            }

                            if (posHijo!=UINT_MAX) //si posHijo!=-1 es decir, si ya se tiene posici�n para el hijo.
                            {
                                if(crossover(conf->listaOrdenFitness[i][j],conf->listaOrdenFitness[i][posLMadre],posHijo,conf->super,conf->promediarPeso,conf->porcentEnableds, conf)==0) //crossover(indexpobPadre,indexpobMadre,indexpobHijo,super,promediarPob);
                                {
                                    fclose(conf->logFile);
                                    conf->logFile=fopen(conf->fileNameLog,"a+");
                                    fprintf(conf->logFile,"<br>\nError 55 en funcion seleccionCruce() llamando a crossover(%u,%u,%u,%1.1f,%1.1f)",conf->listaOrdenFitness[i][j],conf->listaOrdenFitness[i][posLMadre],posHijo,conf->super,conf->promediarPeso);
                                    return(0);
                                }

                                huboReprod=1;//huboReproducci�n=1
                                conf->pob[conf->listaOrdenFitness[i][j]].numHijos--;//padre.numHijos--

                                conf->pob[posHijo].numHijos=-1;//genomaindexpobHijo.numhijos=-1
                            }

                        }
                        else //sino (no se encontr� en especie madre con numhijos>0)
                        {
                            posLMadre=0;//posMadre=champion
                            if (j>0)
                            {
                                posHijo=conf->listaOrdenFitness[i][j];//posHijo=posPadre
                                //crossover(indexpobPadre,indexpobMadre,indexpobHijo,super,promediarPob);
                                if(crossover(conf->listaOrdenFitness[i][j],conf->listaOrdenFitness[i][posLMadre],posHijo,conf->super,conf->promediarPeso,conf->porcentEnableds, conf)==0) //crossover(indexpobPadre,indexpobMadre,indexpobHijo,super,promediarPob);
                                {
                                    fclose(conf->logFile);
                                    conf->logFile=fopen(conf->fileNameLog,"a+");
                                    fprintf(conf->logFile,"<br>\nError 56 en funcion seleccionCruce() llamando a crossover(%u,%u,%u,%1.1f,%1.1f)",conf->listaOrdenFitness[i][j],conf->listaOrdenFitness[i][posLMadre],posHijo,conf->super,conf->promediarPeso);
                                    return(0);
                                }

                                huboReprod=1;//huboReproducci�n=1
                                conf->pob[conf->listaOrdenFitness[i][j]].numHijos--;//padre.numHijos--
                                conf->pob[posHijo].numHijos=-1;//genomaindexpobHijo.numhijos=-1
                            }//MODIFICADO PARA EVITAR REEMPLAZAR A CHAMPION
                        }
                    }
                }
            }
        }
    }
    // realiza cruce inter-especies.
    if ((((float)randL(conf))<(conf->probInterSp*conf->sizePob))&&(conf->numEspecies>1))
    {
        // escoge una especie al azar
        especie1 = (conf->numEspecies-1)*((float)randL(conf));
        // busca su especie m�s cercana
        especieCercana = especieMinDist(conf->representantes[especie1],conf->c1,conf->c2,conf->c3,conf->eG_t,conf);
        // busca en pob un genoma que no sea el rep de espCercana
        i = 0;
        while (((i < conf->sizePob)&&(conf->pob[i].especie!=especieCercana))||(i==conf->representantes[especieCercana]))
            i++;
        if (i>=conf->sizePob)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 55.9 en seleccionCruce() no hay suficientes integrantes de la especie.");
        }
        else
        {
            // hace cruce de rep[esp1] y buscadoesp2 y deja al hijo en buscadoesp2
            if(crossover(conf->representantes[especie1],i,i,conf->super,conf->promediarPeso,conf->porcentEnableds, conf)==0) //crossover(indexpobPadre,indexpobMadre,indexpobHijo,super,promediarPob);
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 55.92 en funcion seleccionCruce() llamando a crossover(%u,%u,%u,%1.1f,%1.1f)",conf->listaOrdenFitness[i][j],conf->listaOrdenFitness[i][posLMadre],posHijo,conf->super,conf->promediarPeso);
                return(0);
            }
            // busca en pob un genoma que no sea el rep de especie1
            i = 0;
            while (((i < conf->sizePob)&&(conf->pob[i].especie!=especie1))||(i==conf->representantes[especie1]))
                i++;
            if (i>=conf->sizePob)
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 55.93 en seleccionCruce() durante cruce intersp.");
                return(0);
            }
            // hace cruce de rep[espCercana] y buscadoEsp1 y deha al hijo en buscadoEsp1
            if(crossover(conf->representantes[especieCercana], i, i, conf->super, conf->promediarPeso,conf->porcentEnableds, conf)==0) //crossover(indexpobPadre,indexpobMadre,indexpobHijo,super,promediarPob);
            {
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>\nError 55.94 en funcion seleccionCruce() llamando a crossover(%u,%u,%u,%1.1f,%1.1f)",conf->listaOrdenFitness[i][j],conf->listaOrdenFitness[i][posLMadre],posHijo,conf->super,conf->promediarPeso);
                return(0);
            }
            // imprimeISP
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"I");
        }
    }
    // libera memoria usada por todos los punteros
    for (i=0; i<conf->numEspecies; i++)
    {
        if (conf->listaOrdenFitness[i]!=NULL) free(conf->listaOrdenFitness[i]);
    }
    if (conf->listaOrdenFitness!=NULL) free(conf->listaOrdenFitness);
    if (conf->fitnessAvgPorEspecie!=NULL) free(conf->fitnessAvgPorEspecie);
    if (conf->actNumGenomasPorEspecie!=NULL) free(conf->actNumGenomasPorEspecie);
    // si no hubo errores, retorna 1
    return(1);
}

float especiacion(TConfig* conf)  // OPTIMIZADA, NMR //TODO: optimizar punteros
{
//realiza la asignaci�n de especies para toda la pob, el exterminio y tambi�n actualiza la liste de representantes de especies y el contador de generaciones sin mejora.
//copia el representante de cada especie a la posici�n sizePob+especie de la pob, sie el fitness es superior
//en caso contrario, copia el representante guardado en lugar del representante actual para conservar los cambios.
//para permitir su mutaci�n y reproducci�n normal en la poblaci�n sin perder info del mejor fitness.
//Se debe realizar DESPUES de la evaluaci�n.
//Altera el valor de threshold en un m�ximo porcentaje especificado para alcanzar el n�mero de especies requerido.
//si el n�mero de especies requerido ya se alcanz�,  no crea nuevas especies sino que asigna a caga genoma la especie m�s cercana.
//actualiza en cada ejecuci�n la lista de conf->representantes y el n�mero de generaciones sin mejora de fitness por especie.
//tambi�n actualiza la lista de especies en conservaci�n.
//Par�metros:	sPob = n�mero de genomas en la conf->poblaci�n
//				numEspDeseadas = n�mero de especies deseadas
//				numConserv = n�mero de especies en conservaci�n.
//				threshold = threshold inicial
//				porcentVarTh = entre 0 y 1 = porcentaje de variaci�n del threshold si no se ha alcanzado el m�ximo n�mero de especies.
//				c1= constante de proporcionalidad en distencia para n�mero de genes excess entre los padres
//				c2= constante de proporcionalidad en distancia para n�mero de genes disjounsigned entre los padres
//				c3= constante de proporcionalidad en distancia para promedio de diferencias de pesos en matching genes de los padres
//				eG_t=(creo que no es necesario REVISAR) n�mero de genes necesarios para considerar el genoma suficientemente grande y hacer n=1
//retorna 0 si hay error, 1 si ok
    unsigned i;
    unsigned j=0;
    unsigned mejor=0;
    //TODO : FALTA: Garantizar randomizando pesos que la distancia al m�s cercano sea mayora la m�nima
    //      se puede utilizar u ciclo hasta m�ximo de veces o >th y retornar el m�s lejano encontrado
    //      y si es mayor del m�nimo aumenta el threshold, sin� lo disminuyen para mentener
    //      el threshold de distancia �ptimo para cualquier n�mero de especies.
    // quitar este decremento,  solo copiar el genoma a cada (ver primeraGEn y exterminio)
    // si conf->numEspecies<conf->numEspecies deseadas decrementa el threshold en porcentVarTh (para que aparezcan como nuevas al ser diferentes a las actuales)
    if (conf->numEspecies < conf->spEspecies)
    {
        conf->threshold = conf->threshold * (1-conf->porcentVarTh);
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br> Eth=%2.2f",conf->threshold);
    }
    // TODO: Verificar si esto es necesario, sino, quitarlo.
    if (conf->numEspecies > conf->spEspecies) 	//AGREGAO HOY TAMBI�N: //TODO:manejar sin modificaci�n del threshold cuando conf->numEspecies>spEspecies para conservar fijo el n�mero de especies si se extermina una especie y se distribuye su poblaci�n de alguna manera/
    {
        conf->threshold = conf->threshold * (1+conf->porcentVarTh);
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nIncrementando threshold de compatibilidad de especies, nuevo valor = %3.3f",conf->threshold);
    }
    //TODO: verificar si se puede hacer el n�mero de especies completamente din�mico controlando el n�mero de especies con el threshold y agrgando o extinguiendo una especie.
    // incrementa el n�mero de generaciones sin mejora para todas las especies
    for (i = 0; i < conf->numEspecies; i++)
        conf->contGeneracSinMejora[i]++;
    // recalcula los representantes
    for (i = 0; i < conf->sizePob; i++)
    {
        if (conf->pob[i].fitness>conf->pob[conf->representantes[conf->pob[i].especie]].fitness)
        {
            //conf->representantes[conf->pob[i].especie]=i
            conf->representantes[conf->pob[i].especie]=i;
            conf->contGeneracSinMejora[conf->pob[i].especie]=0;
        }
    }

    // para toda la pob us asignarEspecies() para actualizar la especie a la mas cercana para el genoma, excepto para los reps.
    for (i=0; i<conf->sizePob; i++)
    {
        //si conf->numEspecies<conf->spEspecies
        if (conf->numEspecies<conf->spEspecies)
        {
            //asignarEspecie(i)
            //si i no es un representante de especie:
            if (conf->representantes[conf->pob[i].especie]!=i)
                conf->pob[i].especie=asignarEspecie(i,conf->threshold,conf->c1,conf->c2,conf->c3,conf->eG_t,conf);
        }
        else
        {
            //sino
            //conf->pob[i].especie = epecieMindist
            //si i no es un representante de especie:
            if (conf->representantes[conf->pob[i].especie]!=i)
                conf->pob[i].especie=especieMinDist(i,conf->c1,conf->c2,conf->c3,conf->eG_t,conf);
            //si conf->pob[i].fitness>conf->pob[conf->representantes[conf->pob[i].especie]].fitness
        }
    }

//	for (i=0;i<conf->numEspecies;i++)
//		conf->conservacionEsp[i]=i;
    //Ordena por fitness el arreglo conf->representantes en el arreglo conservaci�n, donde el index 0 es el de mayor fitness
    //para obtener en conservaci�nEsp en el index 0 el index de especie que tenga m�s fitness y decrementalmente hasta conf->numEspecies-1
    //si dos conf->representantes tienen el m�smo fitness, la mejor posici�n es la del que tenga menor totalConexiones
    /*
    	while(huboSwap){ //repite hasta que el arreglo conf->conservacionEsp est� ordenado (por fitness)
    		//recorre haciendo burbuja down de 0 a conf->numEspecies-1
    		huboSwap=0;
    		if (sentido==1){
    			for (i=0;i<(conf->numEspecies-1);i++){
    				if (conf->pob[conf->representantes[conf->conservacionEsp[i]]].fitness<conf->pob[conf->representantes[conf->conservacionEsp[i+1]]].fitness){
    					huboSwap=1;
    					swap(&(conf->conservacionEsp[i]),&(conf->conservacionEsp[i+1]));
    				}
    				if (conf->pob[conf->representantes[conf->conservacionEsp[i]]].fitness==conf->pob[conf->representantes[conf->conservacionEsp[i+1]]].fitness){
    					if(conf->pob[conf->representantes[conf->conservacionEsp[i]]].totalConexiones>conf->pob[conf->representantes[conf->conservacionEsp[i+1]]].totalConexiones){
    						huboSwap=1;
    						swap(&(conf->conservacionEsp[i]),&(conf->conservacionEsp[i+1]));
    					}
    				}
    			}
    		}
    		else{
    			//recorre haciendo burbuja up de conf->numEspecies-1 hasta 0
    			for (i=(conf->numEspecies-2);i>=0;i--){
    				if (conf->pob[conf->representantes[conf->conservacionEsp[i]]].fitness<conf->pob[conf->representantes[conf->conservacionEsp[i+1]]].fitness){
    					huboSwap=1;
    					swap(&(conf->conservacionEsp[i]),&(conf->conservacionEsp[i+1]));
    				}
    				if (conf->pob[conf->representantes[conf->conservacionEsp[i]]].fitness==conf->pob[conf->representantes[conf->conservacionEsp[i+1]]].fitness){
    					if(conf->pob[conf->representantes[conf->conservacionEsp[i]]].totalConexiones>conf->pob[conf->representantes[conf->conservacionEsp[i+1]]].totalConexiones){
    						huboSwap=1;
    						swap(&(conf->conservacionEsp[i]),&(conf->conservacionEsp[i+1]));
    					}
    				}
    			}

    		}
    		sentido=sentido*-1;

    	}
    		//coloca en 0 el conf->contGeneracSinMejora[] para las especies que se encuentran en conservaci�n, para que cuando salgan de conservaci�n
    		//tengan mejores probabilidades de mejorar antes de ser eliminadas.
    		//fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nEspecies en Conservacion = ");
    		for (i=0;i<conf->numConservacion;i++){
    			conf->contGeneracSinMejora[conf->conservacionEsp[i]]=0;
    			fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>%u(%2.5f), ",conf->conservacionEsp[i],conf->pob[conf->representantes[conf->conservacionEsp[i]]].fitness);
    		}
    		//Si las especies superan maxGenSinMejora coloca contGeneracSinMejora en 0 y coloca todos los genomas dde la espeie i como genomas Iniciales
    		*/
    // verifica si alguna especie debe se exterminada excepto la que tiene el mejor fitness.
    if ((mejor=buscarMejorFitness(conf)) == UINT_MAX)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 34.45 en Funcion especiacion() llamando a buscarMejorFitness()");
        return(0);
    }
    for (i=0; i<conf->numEspecies; i++)
    {
        if (conf->pob[mejor].especie!=i)
        {
            if (conf->contGeneracSinMejora[i]>conf->maxGeneracSinMejora)
            {
                // coloca generacSinMejora en 0
                fclose(conf->logFile);
                conf->logFile=fopen(conf->fileNameLog,"a+");
                fprintf(conf->logFile,"<br>X%2u ",i);
                conf->contGeneracSinMejora[i]=0;
                // sustituye el representante de la especie por un genoma inicial
                if (genomaInicial(conf->representantes[i],conf->numEntradas,conf->numSalidas,conf->numBias,0,i,conf)==0)
                {
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>\nError 34.5 en Funcion especiacion() llamando a genomaInicial(%u,%u,%u)",conf->numEntradas,conf->numSalidas,conf->numBias);
                    return(0);
                }
                //TODO : FALTA: Garantizar randomizando pesos que la distancia al m�s cercano sea mayora la m�nima
                //      se puede utilizar u ciclo hasta m�ximo de veces o >th y retornar el m�s lejano encontrado
                //      y si es mayor del m�nimo aumenta el threshold, sin� lo disminuyen para mentener
                //      el threshold de distancia �ptimo para cualquier n�mero de especies.
                if (randomizarPesos(conf->representantes[i],conf->pesoMinInicial,conf->pesoMaxInicial,conf)==0)
                {
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>\nError 35 en Funcion primeraGen() llamando a randomizarPesos(%u,%1.1f,%1.1f)\n",j,conf->pesoMinInicial,conf->pesoMaxInicial);
                    return(0);
                }
                // muta el genoma inicial AN+AC mutacionesPorExterminio Veces.
                for (j=0; j<conf->mutacionesPorExterminio; j++)
                {
                    // aplica mutaci�n AN al genoma copiado
                    if (mutarAN(conf->representantes[i],conf)==0)
                    {
                        fclose(conf->logFile);
                        conf->logFile=fopen(conf->fileNameLog,"a+");
                        fprintf(conf->logFile,"<br>\nError 34.6 en funcion especiacion() llamando a mutarAN(%u)\n",i);
                        return(0);
                    }
                    // aplica mutaci�n AC al genoma copiado
                    if (mutarAC(conf->representantes[i],conf->maxIntentosMutarAC,conf)==0)
                    {
                        fclose(conf->logFile);
                        conf->logFile=fopen(conf->fileNameLog,"a+");
                        fprintf(conf->logFile,"<br>\nError 34.7 en funcion especiacion() llamando a mutarAC(%u)\n",i);
                        return(0);
                    }
                }
                // coloca la especie = a la desaparecida
                conf->pob[conf->representantes[i]].especie=i;
                // copia el genoma representante de la nueva especie en todos los dem�s genomas con especie == i
                for (j=0; j<conf->sizePob; j++)
                {
                    if((conf->pob[j].especie==i)&&(j!=conf->representantes[i]))  //si no es el represeentante
                    {
                        //copia el genoma del representante a cada genoma de la poblaci�n perteneciente a la especie i

                        copiarGenoma(conf->representantes[i],j,conf);

                        //perturba pesos aleatoriamente para los genomas copiados
                        //TODO: probando mutar pesos totalmente despu�s de extinsi�n en lugar de perturbar los del genoma rep.
                        perturbarPeso(i,conf->porcentMutPeso,conf->probMutPeso,0,0,0,0, conf);
                        /*       if (randomizarPesos(j,conf->pesoMinInicial,conf->pesoMaxInicial,conf)==0){
                                   fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 35 en Funcion primeraGen() llamando a randomizarPesos(%u,%1.1f,%1.1f)\n",j,conf->pesoMinInicial,conf->pesoMaxInicial);
                                   return(0);
                               }   */

                    }
                }

                //guarda una nueva copia del genoma representante.
                copiarGenoma(conf->representantes[i],conf->sizePob+i,conf);
                // TODO: es necesiario que si hubo exterminio de una especie, realiza evaluaci�n de la poblaci�n para alistarse para el cruce.
                // EvaluarEspecie
                // TODO: quitar si se resolvi� error de fitness>1
/*                if (evaluarEspecie(1,i,conf->maxBufferSize,conf->fileNameGTDv1, conf)==0) ////TODO EL ANTERIORMANTE DICHO PARAMETRO para conf->primero para las evaluaciones y hacer inicializaciones a 0 de los valores cuando se crean las neuronas.
                {
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>\nError 59.55 en funcion especiacion() llamando a evaluarEspecie()");
                    return(0);
                }
                */
            }
        }

        // Actualiza representante de nueva especie
        // Compara entre el genoma guardado y el representante ANTERIOR de cada especie, deja el mejor guardado y en el index del rep.
        if (conf->pob[conf->representantes[i]].fitness < conf->pob[conf->sizePob+i].fitness)
        {
            copiarGenoma(conf->sizePob+i,conf->representantes[i],conf);
        }
        // Busca entre toda la pob los nuevos representantes.
        for (j = 0; j < conf->sizePob; j++)
        {
            if (conf->pob[j].especie == i)
            {
                if (conf->pob[j].fitness > conf->pob[conf->representantes[i]].fitness)
                {
                    conf->representantes[i]=j;
                }
            }
        }
        // Compara entre el genoma guardado y el representante ACTUAL de cada especie, deja el mejor guardado y en el index del rep.
        // Espto se hace por si el index del representante var�a entre una generaci�n y otra mejorando al anterior rep mutado(que debe haber disminuido fitness) y as� no sobreescribir el representantte nuevo.
        if (conf->pob[conf->representantes[i]].fitness > conf->pob[conf->sizePob+i].fitness)
        {
            copiarGenoma(conf->representantes[i],conf->sizePob+i,conf);
        }
    }
    // fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>,numEspecies=%u",conf->numEspecies);
    return(1);
}

unsigned primeraGen(unsigned tamPob, unsigned nEntradas, unsigned nSalidas, unsigned nBias, float minPeso, float maxPeso, unsigned maxMutacionesAN,unsigned maxMutacionesAC,unsigned maxIntentosMutAC, short unsigned useMutarAC, short unsigned useMutarAN,unsigned useRandomization, TConfig* conf)  //OPTIMIZADA
{
//genera una poblaci�n inicial de genomas de nIn entradas, nOut salidas y nBias bias a  partir de la mutaci�n
// (AN + AC) de un genoma inicial (en indice 0) totalmente conectado y que es el representante de la especie 0
//Adem�s muta aleatoriamente los pesos de llas conexiones de todos los genomas.
//retorna 0 si hubo alg�n error en el posicionamiento en memoria de las estructuras de los genomas.
//necesita: funcion (//TODO) mutarAN(), mutar(AC), genomaInicial(), nuevoNodo(), nuevaConex(), especieMinDist(),asignarEspecie()
////TODO: se deben hacer dos arreglos para los genomas de los conf->representantes de cada especie de la generaci�n actual y para
//los conf->representantes de las especies de la generaci�n anterior.
//En cada generaci�n se debe comparar cada genoma con los conf->representantes de la generaci�n anterior para determinar la especie a la que pertenecen
//y se toma el genoma con menor error como representante de la especie en el arreglo actual.
//Se debe llevar un record para cada representante de cada especie del n�mero de generaciones que lleva sin mejorar, para poder eliminaar de esta
//manera especies que se quedaron estancadas(EXCEPTO LA MEJOR(o n mejores?)).
//El n�mero de hijos que puede producir cada especie depende del fitness de sus individuos comparado con el promedio de fitness total como se
//muestra en la pag 394 del libro de AI game programming. (great tool)
//Para la primera generaci�n: Averiguar
//Algoritmo para especiar toda la conf->poblaci�n despu�s de mutaci�n y cruce, en pag 54 de disertaci�n de PhD
    unsigned i;
    unsigned j;
    unsigned numMutacionesAN;
    unsigned numMutacionesAC;
    unsigned mejor=conf->maxIntentosDist;
//TODO : FALTA: Garantizar randomizando pesos que la distancia al m�s cercano sea mayora la m�nima
//      se puede utilizar u ciclo hasta m�ximo de veces o >th y retornar el m�s lejano encontrado
//      y si es mayor del m�nimo aumenta el threshold, sin� lo disminuyen para mentener
//      el threshold de distancia �ptimo para cualquier n�mero de especies.
// para esto llenar la poblaci�n de genomas, randomizar pesos de todos escoge como rep de esp 0 a 0
// luego busca entre los que no son reps de ninguna especie el que tenga mayor m�nima distancia hasta todos los representantes
// y lo hace representante de la siguiente especie hasta spEspecies. En cada iteraci�n se DEBE hacer una randomizaci�n de los
// no-representantes (pues son m�s cercanos que el nuevo rep). DEJAR el thresold de distancia inicializado en la m�nima (porque puede haber distancias muy grandes si hay mutaciones)
// encontrada entre cada rep y los dem�s, luego hacer control sobre el th durante cada exterminio para mantenerse en el m�ximo posible
// realizando conf->iterSelRand iteraciones de randomizaci�n de pesos (y mutaciones Iniciales), manteniendo el genoma de
// MAYOR distancia a cada representante y ubicandolo como representante de su especie
    //busca el mayor entre maxIntentosDist ymaxIntentosDistInicial
    if (conf->maxIntentosDist<conf->maxIntentosDistInicial)
        mejor=conf->maxIntentosDistInicial;
    // reserva memoria para vector de genomas conf->pob lo crea de tama�o tamPob+spEspecies para poder guardar en los �ltimos elementos los representantes de cada especie(adicionado en versi�n 0.59).
    // se adiciona 1 para usar el �ltimo como genoma
    conf->realSizePob=mejor+(2*tamPob)+conf->spEspecies+1;
    if (inicializarPob(conf->realSizePob,nEntradas,nSalidas,nBias,conf)==0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 33 en Funcion primeraGen() llamando a incializarconf->pob(%u,%u,%u,%u)\n",tamPob,nEntradas,nSalidas,nBias);
        return(0);
    }
    //marca la posici�n del genoma temporal (usado para procesamiento distribuido.)
    conf->tmpIndexPob=mejor+(2*tamPob)+conf->spEspecies;
    //Crea el genoma inicial en conf->pob[0]
    if (genomaInicial(0, nEntradas, nSalidas, nBias, 1, 0, conf)==0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 34 en Funcion primeraGen() llamando a genomaInicial(%u,%u,%u)\n",nEntradas,nSalidas,nBias);
        return(0);
    }
    //randomiza pesos de genoma inicial
    if (useRandomization==1)
        if (randomizarPesos(0,minPeso,maxPeso,conf)==0)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 35 en Funcion primeraGen() llamando a randomizarPesos(%u,%1.1f,%1.1f)\n",0,minPeso,maxPeso);
            return(0);
        }

    //imprimirGenoma(0,conf);
    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>\nClonando genoma inicial...\n");
    //Para los dem�s genomas copia del original, randomiza y hace las dos mutaciones (DEBE incluir los genomas guardados de reps de cada especie).
    for (i=1; i<((tamPob*2)+conf->spEspecies); i++)
    {
        //Copia el genoma inicial(0) a los dem�s desde 1 hasta tamPob-1
        if (copiarGenoma(0,i,conf)==0)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 36 en Funcion primeraGen() llamando a copiarGenoma(0,%u)\n",i);
            return(0);
        }
    }
//PROBANDO MAXDISTANCIA, quitar si ok.
    //genomaMasLejano(1,conf->maxIntentosDistInicial,conf);
    //fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>dist=%7.7f, espCerc=%7.7f",calcularDist(1, conf->representantes[conf->pob[1].especie], conf->c1, conf->c2, conf->c3, conf->eG_t,conf)
    //       ,distEspecieCercana( 1,conf->c1, conf->c2, conf->c3, conf->eG_t,conf)/2);

    fclose(conf->logFile);
    conf->logFile=fopen(conf->fileNameLog,"a+");
    fprintf(conf->logFile,"<br>Aplicando mutaciones iniciales...\n");
    //hace maxDist, randomiza y hace las dos mutaciones (DEBE incluir los genomas guardados de reps de cada especie).
    for (i=1; i<((tamPob*2)+conf->spEspecies); i++)
    {

        // si se est� usando mzxDisteancia, se calcula el nuevo nodo m�s distante a todos los existentes
        // debe ir antes de todas las otras mutaciones y es excluyente con userandomization
        if (conf->maxIntentosDistInicial>0)
        {
            genomaMasLejano(i, conf->maxIntentosDistInicial,conf);
        }
        else
        {
            //randomiza pesos de conexiones de cada genoma copiado
            if (useRandomization==1)
            {
                if ((randomizarPesos(i,minPeso,maxPeso,conf))==0)
                {
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>\nError 37 en Funcion primeraGen() llamando a randomizarPesos(%u,%1.1f,%1.1f)\n",i,minPeso,maxPeso);
                    return(0);
                }
            }
        }

        //Aplica mutaci�n AN
        if (useMutarAN==1)
        {
            numMutacionesAN=(unsigned )floor((((float)randL(conf))*maxMutacionesAN) + 0.5);
            for (j=0; j<numMutacionesAN; j++)
                if (mutarAN(i,conf)==0)
                {
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>\nError 38 en funcion primeraGen() llamando a mutarAN(%u)\n",i);
                    return(0);
                }
        }

        //Aplica mutaci�n AC
        if (useMutarAC==1)
        {
            numMutacionesAC=(unsigned )floor((((float)randL(conf))*maxMutacionesAC) + 0.5);
            for (j=0; j<numMutacionesAC; j++)
                if (mutarAC(i,maxIntentosMutAC,conf)==0)
                {
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>\nError 39 en funcion primeraGen() llamando a mutarAC(%u,%u)\n",i,maxIntentosMutAC);
                    return(0);
                }
        }
    }
    return(1);
}

unsigned inicializarPob(unsigned tamPob,unsigned nEntradas,unsigned nSalidas,unsigned nBias, TConfig* conf)  // OPTIMIZADA, NMR,
{
//Inicializa el vector de genomas (conf->poblacci�n)
//tambi�n obtiene memoria para arreglo de nodos y conexiones de tama�o nEntradas*nSalidas*nBias
//Obtiene memoria para cada nueva estructura Genoma.
//Retorna 0 si hay error
    unsigned i=0;
    // Reserva memoria para randList[conf->tamRandList]
    conf->randList=(tRandList*) malloc(conf->tamRandList*sizeof(tRandList));
    printf("RandListInit\n");
    // Inicializa los valores aleatorios entre 0 y 1 en randList
    for(i=0;i<(conf->tamRandList-1);i++)
    {
        conf->randList[i].valor=(float)rand()/RAND_MAX;
        conf->randList[i].next=&(conf->randList[i+1]);//siguiente
    }
    // arregla el siguiente del �ltimo para que sea el primero.
    conf->randList[conf->tamRandList-1].valor=(float)rand()/RAND_MAX;
    conf->randList[conf->tamRandList-1].next=&(conf->randList[0]);//siguiente
    // apunta el puntero global de posici�n en lista de aleatorios al primero.}
    conf->punteroRand=&(conf->randList[0]);
    i=0;
    // reserva memoria para POB
    printf("PobListInit\n");
    if (conf->pob!=NULL) free((void *)conf->pob);
    if ((conf->pob = ( Genoma *) calloc(tamPob,(unsigned  int)sizeof(Genoma)))==NULL)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError 28 en funcion inicializarPob(%u,%u,%u,%u) llamando a calloc(%u,%u)\n",tamPob,nEntradas,nSalidas,nBias,tamPob,(unsigned  int)sizeof(Genoma));
        return(0);
    }
    for (i=0; i<tamPob; i++)
    {
        if (conf->pob[i].nodo!=NULL) free((void *)conf->pob[i].nodo);
        if ((conf->pob[i].nodo=(GenNodoF*)calloc(1,sizeof(GenNodoF)*(nEntradas+nBias+nSalidas)))==NULL)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 29 en funcion inicializarPob(%u,%u,%u,%u) llamando a calloc(1,%u)\n",tamPob,nEntradas,nSalidas,nBias,(nEntradas+nBias+nSalidas)*(unsigned int)sizeof(GenNodoF));
            return(0);
        }
        if (conf->pob[i].conex!=NULL) free((void *)conf->pob[i].conex);
        if ((conf->pob[i].conex=(GenConexF*)calloc(1,sizeof(GenConexF)*((nEntradas+nBias)*nSalidas)))==NULL)
        {
            fclose(conf->logFile);
            conf->logFile=fopen(conf->fileNameLog,"a+");
            fprintf(conf->logFile,"<br>\nError 30 en funcion inicializarPob(%u,%u,%u,%u) llamando a calloc(1,%u)\n",tamPob,nEntradas,nSalidas,nBias,(nEntradas+nBias+nSalidas)*(unsigned int)sizeof(GenConexF));
            return(0);
        }

    }
    return(1);
}

int competencia(TConfig* conf)
// Realiza competencia, que copmara con una copia guardada en indexPob+numEspecies+sizePob
// usa variable conf->porcentCompetencia como prob de ganancia del backup si es mejor que el actual,
// si es peor, simplemente se reemplaza.
// si es negativa, no hace competencia.
// solo compara si los correspondientes indexpob son de la misma especie( para evitar acaparamiento)
// se debe ejecutar despu�s de especiacion().
// retorna 0 si hay error, 1 si ok
{
    int i;
    unsigned bkp;
    if (conf->porcentCompetencia<0)
    {
        fclose(conf->logFile);
        conf->logFile=fopen(conf->fileNameLog,"a+");
        fprintf(conf->logFile,"<br>\nError: porcentCompetencia negativo");
        return(0);
    }
    //para todos los genomas de la poblaci�n:
    for(i=0; i<conf->sizePob; i++)
    {
        bkp=conf->sizePob+i+conf->spEspecies; // posici�n del backup
        // si i y backup son de la misma especie
        if (conf->pob[i].especie==conf->pob[bkp].especie)
        {
            // si backup.fitness>i.fitness
            if (conf->pob[bkp].fitness>conf->pob[i].fitness)
            {
                // prob de porcentCompetencia
                if (((float)randL(conf))<conf->porcentCompetencia)
                {
                    // copia backup a i
                    fclose(conf->logFile);
                    conf->logFile=fopen(conf->fileNameLog,"a+");
                    fprintf(conf->logFile,"<br>S=%u,",i);
                    copiarGenoma(bkp,i,conf);
                }
                else //TODO: PROBAR QUITANDO EL ELSE
                {
                    // copia i malo a bkp (como porque puede tenerse que pasar por etapa de malo durante mutaci�n y luego mejorar)
                    copiarGenoma(i,bkp,conf);
                }

            }
            else
            {
                copiarGenoma(i,bkp,conf);
            }

        }
        else //si son de diferentes especies en el mismo index, copia el nuevo a backup
        {
            copiarGenoma(i,bkp,conf);
        }
    }
    return(1);
}
