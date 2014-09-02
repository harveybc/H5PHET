/** H5PHET - SNN Evaluator.
	Evalúa una red neuronal especificada en un archivo SNN (Simple Neural Network) con un
	archivo de entrenamiento GTD (Generic Training Data) y muestra el error obtenido por la red.
	Parte de H5PHET.
	Por Harvey Bastidas.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// estructura del encabezado de GTDv1 (Generic Training Data).
// formato GDTv1: [encabezado][d1_entrada_0]..[d1_entrada_n][d1_salida_0]..[d1_salida_n].._entrada_0]..[dn_entrada_0]..[dn_entrada_n][dn_salida_0]..[dn_salida_n]
// estructura del encabezado de SNNv1 (Simple Neural Network).
// formato SNNv1: [encabezado][unsigned conexIn[numConex]][unsigned conexOut[numconex]][char enabled[numconex]][float/double peso[numConex]]

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
// formato SNNv1: [encabezado][unsigned conexIn[numConex]][unsigned conexOut[numconex]][float/double peso[numConex]]
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
 	double lastFitness;
} hdrSNNv1;

unsigned* conexIn;
unsigned* conexOut;

double* pesoD;
float* pesoF;
char* enabled;
hdrSNNv1 headerSNN;
hdrGTDv1 headerGTD;
float** dataGTDf; // buffer usado para almacenar los datos de entrenamiento en formato GTD
float* Vc; // buffer para valores calculados, usado para calculo de fitness con correlación Pearson.

/******************************************************************************************/


float correlac(float* Vc,float** dataGTDf, hdrGTDv1 headerGTD, int numDatos)
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
    for (i=0;i<numDatos;i++)
    {
        xm+=dataGTDf[i][headerGTD.numEntradas];
        ym+=Vc[i];
    }
    xm/=numDatos;
    ym/=numDatos;
    // calcula sum0(0,nD,(Xi-Xm)*(Yi-Ym)),sum1(0,nD,sqr(Xi-Xm)) y sum2(0,nD,sqr(Yi-Ym)
    for (i=0;i<numDatos;i++)
    {
        sum0+=((dataGTDf[i][headerGTD.numEntradas]-xm)*(Vc[i]-ym));
        sum1+=((dataGTDf[i][headerGTD.numEntradas]-xm)*(dataGTDf[i][headerGTD.numEntradas]-xm));
        sum2+=((Vc[i]-ym)*(Vc[i]-ym));
    }
    // retorna la correlación
    return(sum0/(sqrt(sum1)*sqrt(sum2)));
}


/**************************************   fError        ********************************************************/

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

/**************************************   fSigma        ********************************************************/

float fSigma(float fX, float A, unsigned usarSigned)
{
	if (usarSigned==1)
	{

		return 2*(exp(-A*(fX * fX)))-1;
	}
	else
	{
		return(exp(-A*(fX*fX)));
	}
}
/**************************************   fSigma        ********************************************************/

int buscarConex(int cxIn, int cxOut)
{
    int i;
    for(i=0;i<headerSNN.numConex;i++)
        if ((conexIn[i]==cxIn)&&(conexOut[i]==cxOut))
            return(i);
    printf("Error 0.001: la conex no existe.");
    exit(0);
}

/*************************************    calcularValor          ***************************************************************/


///**********************************  MAIN **************************************///

int main()
{
    float* valor; //hasta 10k neuronas
    //char valorCalculado[100];
    float error=0;;
    float fitness=0;;
    int numDatos=0;
    int i,k;
    float acum;
    unsigned j;
    int leidos;
    FILE *fileIn;
    unsigned totalNodos;
    unsigned tmpSize=0;
    int tamListaConex;
    int* listaConex;
    // abre el archivo de salida para escritura
    if ((fileIn=fopen("/home/harveybc/Dropbox/h5phet/mejor.snn\0","rb"))==NULL)
    {
        printf("\nError 48 llamando a fopen()\n");
        return(0);
    }
    // lee el encabezado int numEntradas,numBias,numSalidas,numHiddens,numConex
    leidos = fread(&headerSNN, sizeof(hdrSNNv1), 1, fileIn);
    // calcula el total de nodos.
    totalNodos = headerSNN.numEntradas+headerSNN.numBias+headerSNN.numSalidas+headerSNN.numHiddens;
    // ubica memoria para los arreglos de conexionesIO y pesos.
    tmpSize = headerSNN.numConex*sizeof(unsigned);
    if ((conexIn = (unsigned*) malloc(tmpSize*sizeof(unsigned)))==NULL)
        return(0);
    if ((conexOut = (unsigned*) malloc(tmpSize*sizeof(unsigned)))==NULL)
        return(0);
    if ((enabled = (char*) malloc(tmpSize*sizeof(char)))==NULL)
        return(0);
   // si tipo de dato de pesos ees float
    tmpSize = headerSNN.numConex*sizeof(float);
    if (headerSNN.tamRegistros==4)
    {
        pesoF = malloc(tmpSize*sizeof(float));
    }
    printf("Encabezado:\numEntradas=%i,numBias=%i,numSalidas=%i,numHiddens=%i,numConex=%i,sigmaFactor=%3.3f,usarSigned=%d,actTreshold=%3.3f\n",headerSNN.numEntradas,headerSNN.numBias,headerSNN.numSalidas,headerSNN.numHiddens,headerSNN.numConex,headerSNN.sigmaFactor,headerSNN.usarSigned,headerSNN.actThreshold);
    // lee los arreglos en orden: int conexIn[numConex],conexOut[numConex],double peso[numConex]
    leidos = fread(conexIn, sizeof(unsigned), headerSNN.numConex, fileIn);
    leidos += fread(conexOut, sizeof(unsigned), headerSNN.numConex, fileIn);
    leidos += fread(enabled, sizeof(char), headerSNN.numConex, fileIn);
    leidos+=fread(&tamListaConex,sizeof(int),1,fileIn);
    printf("\ntamListaConex=%d",tamListaConex);
    listaConex=(int*)malloc(tamListaConex*sizeof(int));
    leidos+=fread(listaConex,sizeof(int),tamListaConex,fileIn);
    // lee el arreglo de pesos si el tipo de dato es float
    if (leidos!=((3*headerSNN.numConex)+1+tamListaConex))
        exit(0);
    if (headerSNN.tamRegistros==4)
    {
        if (pesoF)
        {
            printf("\nPesos Leidos\n");
            leidos = fread(pesoF, sizeof(float), headerSNN.numConex, fileIn);
            if (leidos!=headerSNN.numConex)
            {
                printf("\nError 102.05 No se pudo leer el arreglo de pesos");
            }
        }
        else
        {
            printf("\nError 102.1 en main() el buffer pesoF es NULL");
            return(0);
        }
    }
    else
    {
        printf("\nError 51: No se pudieron leer los pesos, tamRegistros=%i\n",headerSNN.tamRegistros);
        return(0);
    }
    // imprime cada conexión del genoma leído
    for (i=0;i<headerSNN.numConex;i++)
	{
        if (headerSNN.tamRegistros==4)
        {
            printf("C(%i)%i-%i[%3.3f],",i,conexIn[i],conexOut[i],pesoF[i]);
        }
	}
    //imprime el orden de evaluación
    for (i=0;i<tamListaConex;i++)
	{
        printf("\n---->C[%d]\n",listaConex[i]);
	}

    // imprime error si el número de pesos leídos es diferente al especificado en el header
    if (leidos!=headerSNN.numConex)
    {
        printf("\nError leyendo archivo SNN\n");
        return(0);
    }
    //cierra el archivo SNN.
    if (fclose(fileIn)!=0)
    {
        printf("\nError 49 llamando a fclose() de archivo SNN\n");
        return(0);
    }
    // Abre archivo de entradas y salidas para lectura
    printf("Abriendo archivo de entrenamiento GTD.\n");
    if ((fileIn=fopen("/home/harveybc/Dropbox/h5phet/trdata.gtd\0","rb"))==NULL)
    {
        printf("\nError 50 llamando a fopen()\n");
        return(0);
    }
//Inicio adicion
    // lee headerGTD
    leidos=fread(&headerGTD,sizeof(hdrGTDv1),1,fileIn);
    // busca el filesize
    if (leidos<1)
    {
        printf("\nError 58 en evaluarPob() llamando a fread()");
        return(0);
    }
    // averigua el fileSize y calcula numDatos
    fseek(fileIn, 0, SEEK_END);
    numDatos = (ftell(fileIn) - sizeof(hdrGTDv1))/(headerGTD.tamRegistros*(headerGTD.numEntradas+headerGTD.numSalidas));
    if (numDatos<1)
    {
        printf("\nError 65 en evaluarPob()llamando a ftell()");
        return(0);
    }
    //printf("\nFileSize=%d",ftell(fileIn));
    rewind(fileIn);
    // ubica en posición de lectura de datos después del header.
    leidos=fread(&headerGTD,sizeof(hdrGTDv1),1,fileIn);
    // reserva memoria en conf->dataGTD[numDatos][numEntradas+numSalidas]
    // para todo el buffer, ojo! obtiene toda la memoria para el index 0 para que sean consecutivos
    // durante la lectura
    dataGTDf=(float**)malloc(numDatos*sizeof(float*));
    if (!dataGTDf)
    {
        printf("\nError 64 en evaluarPob() llamando a malloc()");
        return(0);
    }
    dataGTDf[0]=(float*)malloc((headerGTD.numEntradas+headerGTD.numSalidas)*numDatos*sizeof(float));
    if (!dataGTDf[0])
    {
        printf("\nError 65 en evaluarPob() llamando a malloc()");
        return(0);
    }
    // coloca las filas de la matriz del buffer GTD apuntando al inicio de cada grupo de datos.
    for (i=1;i<numDatos;i++)
    {
        //TODO: OJO!!!! Verificar si es sizeof(float)*(nIn+nOut)*i o (nIn+nOut)*i
        dataGTDf[i]=dataGTDf[0]+((headerGTD.numEntradas+headerGTD.numSalidas)*i);
    }
    // lee de GTD (numEntradas+numSalidas)*sizeof(float) en conf->dataGTD[i];
    leidos=fread(dataGTDf[0],sizeof(float)*(headerGTD.numEntradas+headerGTD.numSalidas),numDatos,fileIn);
    // verifica lectura.
    if (leidos<numDatos)
    {
        printf("\nError 66 en evaluarPob() llamando a fread()");
        return(0);
    }
    fclose(fileIn);    // confirma que los datos de GTD y SNN sean compatibles.
// fin adicion
    if ((headerGTD.numEntradas!=headerSNN.numEntradas)||(headerGTD.numSalidas!=headerSNN.numSalidas)||(headerGTD.tamRegistros!=headerSNN.tamRegistros)||(headerGTD.usarSigned!=headerSNN.usarSigned))
    {
    	printf("Error 51.1 los encabezados de SNN y GTD no coinciden.\n");
    	printf("GTD: In = %u , Out = %u , tamReg = %u , usarSigned = %u\n",headerGTD.numEntradas,headerGTD.numSalidas,headerGTD.tamRegistros,headerGTD.usarSigned);
		printf("SNN: In = %u , Out = %u , tamReg = %u , usarSigned = %u, sigmaFactor=%3.3f\n",headerSNN.numEntradas,headerSNN.numSalidas,headerSNN.tamRegistros,headerSNN.usarSigned,headerSNN.sigmaFactor);
    	return(0);
    }
/*  insertado para usar correlación como fitness */

        Vc=(float*)malloc(numDatos*sizeof(float));

/** insertado para hacer evaluar igual a neat */

        totalNodos=headerSNN.numEntradas+headerSNN.numBias+headerSNN.numSalidas+headerSNN.numHiddens;
        // reserva memoria para valor
        valor=(float*)malloc(totalNodos*sizeof(float));
        // inicializa valor y valorAnt
        for (k=0;k<totalNodos;k++)
        {
            valor[k]=0;
        }
        //coloca bias en 1 ://TODO: falta para varias bias, es necesario? NO
        valor[headerSNN.numEntradas]=1;
        fitness=0;
        for (j=0;j<numDatos;j++)
        {
            // coloca los valores de entrada de valor[] en valorAnt[] y copia los nuevos a valor[]
           // memcpy(valorAnt,valor,indInS);
            memcpy(valor,&(dataGTDf[j][0]),(headerSNN.numEntradas*sizeof(float)));
  /*          valor[0]=dataGTDf[j][0];
            valor[1]=dataGTDf[j][1];
            valor[2]=dataGTDf[j][2];
            valor[3]=dataGTDf[j][3];*/
            // coloca como no - calculados todos los valores.

            //calcula fsigma para las entradas
           // printf("\n");
            for (k=0;k<headerGTD.numEntradas;k++)
            {
               // printf("E[%d]=%3.3f ",k,valor[k]);
                valor[k]=2*((float)exp(-(float)headerSNN.sigmaFactor*((valor[k]-(float)headerSNN.actThreshold)*(valor[k]-(float)headerSNN.actThreshold))))-1;;
               //printf("v[%d]=%3.3f",k,valor[k]);
            }

            // coloca bias con valor 1 y calculado.
            //valorCalculado[headerGTD.numEntradas]=1;
            valor[headerGTD.numEntradas]=1;
            // calcula el nuevo valor[] con el genoma i y dataGTD
            acum=0;
            for (k=0;k<tamListaConex;k++)
            //for (k=0;k<headerSNN[i].numConex;k++)
            {
                acum+=(pesoF[listaConex[k]]*valor[conexIn[listaConex[k]]]);

                // si el conexOut no es el último
                if (k<(tamListaConex-1))
                {
                    // si el ConexOut siguiente es diferente al actual
                    if (conexOut[listaConex[k]]!=conexOut[listaConex[k+1]])
                    {
                        // calcula el sigma y actualiza valorAnt
                        valor[conexOut[listaConex[k]]]=2*((float)exp(-(float)headerSNN.sigmaFactor*((acum-(float)headerSNN.actThreshold )*(acum-(float)headerSNN.actThreshold))))-1;
                        //printf("S");
                        //reinicializa el acumulador
                        acum=0;
                    }
                }
                //printf("C%d(p%3.3f*v%3.3f=%3.3f)",k,pesoF[listaConex[k]],valor[conexIn[listaConex[k]]],(pesoF[listaConex[k]]*valor[conexIn[listaConex[k]]]));
            }
            // saca el Fsigma del último nodo.
            //printf("\nva=%3.3f",acum);
            valor[conexOut[listaConex[tamListaConex-1]]]=(2*((float)exp(-(float)headerSNN.sigmaFactor*(acum-(float)headerSNN.actThreshold)*(acum-(float)headerSNN.actThreshold))))-1;
            //printf("v=%3.3f",valor[conexOut[listaConex[tamListaConex-1]]]);
            Vc[j]=valor[conexOut[listaConex[tamListaConex-1]]];
            // copia valor a valorAnt
            // suma los errores de las neuronas de salida
            error=0;
            for (k=(headerSNN.numEntradas+headerSNN.numBias);k<(headerSNN.numEntradas+headerSNN.numBias+headerSNN.numSalidas);k++)
            {
                error+=fError(valor[k],dataGTDf[j][headerGTD.numEntradas]);
               // printf("S0=%3.3f\n",dataGTDf[j][headerGTD.numEntradas]);
//                printf("E0=%1.1f,%1.1f E1=%1.1f,%1.1f E2=%1.1f,%1.1f E3=%1.1f,%1.1f Vt=%3.3f,Vc=%3.3f\n",valor[0],dataGTDf[j][0],valor[1],dataGTDf[j][1],valor[2],dataGTDf[j][2],valor[3],dataGTDf[j][3],dataGTDf[j][headerGTD.numEntradas],valor[k]);
            }
            // acumula en fitness el error encontrado
            //fitness+=error;
            printf("Vt=%3.3f , Vc=%3.3f\n",dataGTDf[j][headerGTD.numEntradas],valor[(headerSNN.numEntradas+headerSNN.numBias)]);
        }
        // calcula el promedio de los errores
        fitness/=numDatos;
        // hace fitness = 1-error promedio
        fitness=correlac(Vc,dataGTDf,headerGTD,numDatos-1);
    printf("\nFitness = %11.11f",fitness);

    // libera memoria de los arreglos.
    free(conexIn);
    free(conexOut);
    free(enabled);

    if (headerSNN.tamRegistros==4)
    {
		free(pesoF);
    }
    return(0);
}

