#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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

double correlac(int x, int y, double** buffer, double** correlacM, double* medias, int numDatos)
// retorna el coeficiente de correlación de dos arreglos de datos x e y de tamaño numDatos.
// la matriz de correlaciones correlacM debe inicializarse cn valores >10 la primera vez.
// parámetros:
//      x,y = vectores a comparar
//      xm,ym = media de x y media de y, se deben proporcionar por optimización de pre-calculo de medias
//      numDatos = tamaño de los vectores
{
    int i=0;

    double sum0=0;
    double sum1=0;
    double sum2=0;
    // si x y y son iguales, retorna 1
    if (x==y)
    {
        correlacM[x][y]=1;
        correlacM[y][x]=1;
        return(1);
    }
    // si ya fué calculada.
    if ((correlacM[x][y]<10)||(correlacM[y][x]<10))
    {
        return(correlacM[y][x]);
    }

    // calcula sum0(0,nD,(Xi-Xm)*(Yi-Ym)),sum1(0,nD,sqr(Xi-Xm)) y sum2(0,nD,sqr(Yi-Ym)
    for (i=0;i<numDatos;i++)
    {
        sum0+=((buffer[i][x]-medias[x])*(buffer[i][y]-medias[y]));
        sum1+=((buffer[i][x]-medias[x])*(buffer[i][x]-medias[x]));
        sum2+=((buffer[i][y]-medias[y])*(buffer[i][y]-medias[y]));
    }
    // retorna sum0/(sqrt(sum1)*sqrt(sum2))
    if (sum1==0)
        printf("\nError en correlacentrada %d ",x);
    if (sum2==0)
        printf("\nError en entrada %d ",y);
        // actualiza la matriz de correlaciones
    correlacM[y][x]=sum0/(sqrt(sum1)*sqrt(sum2));
    correlacM[x][y]=correlacM[y][x];
    return(correlacM[y][x]);
}

double ponderacMultiOut(int index, int* s, double** buffer, double** correlacM, double* medias, hdrGTDv1 headerGTD,int numDatos, int metodoPonderac) // o poción multijugos
// calcula la multip de la correl de index(#entrada) con con las salidas sobre
// la multip de correlac de index con las entradas SELECCIONADAS).
// parámetros:
//      x,y = vectores a comparar
//      xm,ym = media de x y media de y, se deben proporcionar por optimización de pre-calculo de medias
//      numDatos = tamaño de los vectores
//      métodoPonderac = 0=son sumatoria de correlac de entradas en denominador o 1= con multiplicatoria.
{
    unsigned i;
    double mUp=(metodoPonderac==0?0:1);
    double mDn=(metodoPonderac==0?0:1);  //QUEDA? FALTA: MOFIFICADO PARA PRUEBA,ORIGINALMENTE 1 y en Dn es multip no sum

    //verifica param index
    if (index>(headerGTD.numEntradas+headerGTD.numSalidas))
    {
        printf("\nError 12 en ponderacMultiOut() index>numIns+numOuts");
    }
    // calcula mUp= multip(correlac(index,Salidas))
    for (i=0;i<headerGTD.numSalidas;i++)
    {
        if (metodoPonderac==0)
            mUp=mUp+fabs(correlac(s[index], headerGTD.numEntradas+i, buffer, correlacM, medias, numDatos));
        else
            mUp=mUp*fabs(correlac(s[index], headerGTD.numEntradas+i, buffer, correlacM, medias, numDatos));
    }
    // calcula mDn= multip(correlac(index,entradasSeleccionadas)) es decir entre 0 y index
    for (i=0;i<index;i++)
        if (metodoPonderac==0)
            mDn=mDn+fabs(correlac(s[index], i, buffer, correlacM, medias, numDatos));
        else
            mDn=mDn*fabs(correlac(s[index], i, buffer, correlacM, medias, numDatos));
// TODO: PROBANDO: originalmente era mUp/mDn
    return(mUp*mUp/mDn); //siempre es positivo
}

int main()
{
    int numEntradasSelec=24; //Número de entradas con las que se generará el GTD V1.
    int metodoPonderac=0; //0=son sumatoria de correlac de entradas en denominador o 1= con multiplicatoria.
    int numDatos=150000;
    char fileNameIn[300];
    char fileNameOut[300];
    char fileNameLog[300];
    int i, j, result, tmpI ;
    double maxD;
    int leidos;
    FILE* fileIn;
    FILE* fileOut;
    FILE* fileLog;
    hdrGTDv1 header;
    strcpy(fileNameIn,"d:\\h5phet\\trdata.gtd\0");
    strcpy(fileNameOut,"d:\\h5phet\\trdata_cut.gtd\0");
    strcpy(fileNameLog,"d:\\h5phet\\log_selec.html\0");
    int numEntradasAnt;
    int* s;

    double** correlacM;
    double** buffer;
    double* medias;
    double* valorPond;
    float* tmpFloatArr;
    double* tmpDoubleArr;
    double tmp,mUp,mDn;
    double sumaTotal=0;
    // abre GTDv1In para lectura
    printf("\n%s\n",fileNameIn);
    fileIn = fopen(fileNameIn,"rb");
    if (fileIn==NULL)
    {
        printf("Error 1: main() llamando a fopen()\n");
        return(0);
    }
    // lee encabezado GTDv1
    leidos = fread(&header,sizeof(hdrGTDv1),1,fileIn);
    if (leidos<1)
    {
        printf("Error 2: main() llamando a fread()\n");
        return(0);
    }
    if ((header.numEntradas>32000)||(header.numSalidas>32000))
    {
        printf("Error 3: main() numEntradas>32k o numSalidas>32k\n");
        return(0);
    }
    printf("\nH5PHET - Input selector V0.1\nnumEntradas=%d , numSalidas=%d\n, tamRegistros=%d, usarSigned=%d\n",header.numEntradas, header.numSalidas, header.tamRegistros, header.usarSigned);
    // reserva memoria para el arreglo corelac** [numEntradas+numSalidas][numEntradas+numSalidas]
    correlacM = (double**)malloc((header.numEntradas+header.numSalidas)*sizeof(double*));
    for (i=0;i<(header.numEntradas+header.numSalidas);i++)
        correlacM[i] = (double*)malloc((header.numEntradas+header.numSalidas)*sizeof(double));
    // reserva memoria para el arreglo medias* [numEntradas+numSalidas]
    medias = (double*)malloc((header.numEntradas+header.numSalidas)*sizeof(double));
    // reserva memoria para el buffer** [numEntradas+numSalidas][numDatos]
    buffer = (double**)malloc(numDatos*sizeof(double*));
    for (i=0;i<numDatos;i++)
        buffer[i] = (double*)malloc((header.numEntradas+header.numSalidas)*sizeof(double));
    // reserva memoria para valorPond[numEntradas]
    valorPond = (double*)malloc(header.numEntradas*sizeof(double));
    // reserva memoria para el arreglo s[numEntradas]
    s = (int*)malloc(header.numEntradas*sizeof(int));
    // verifica si algúm malloc retornó NULL
    if ((!s)||(!correlacM)||(!buffer)||(!medias)||(!valorPond))
    {
        printf("\nError 3.5: main() llamando a malloc()\n");
        return(0);
    }
    // inicializa la matriz correlac=20
    for (i=0;i<(header.numEntradas+header.numSalidas);i++)
        for(j=0;j<(header.numEntradas+header.numSalidas);j++)
            correlacM[i][j]=20;// para diferenciar las ya calculadas
    // inicializa el arreglo s[i]=i
    for (i=0;i<header.numEntradas;i++)
        s[i]=i;
    // inicializa el vector medias en 0
    for (i=0;i<(header.numEntradas+header.numSalidas);i++)
        medias[i]=0;
    // lee numDatos en el buffer
    if (header.tamRegistros==4)
    {
        tmpFloatArr = malloc((header.numEntradas+header.numSalidas)*sizeof(float));
        //verifica si alguno de los malloc retornó NULL
        if (!tmpFloatArr)
        {
            printf("\nError 3.7: main() llamando a malloc()\n");
            return(0);
        }
        for (i=0;i<numDatos;i++)
        {
            leidos=fread(tmpFloatArr,sizeof(float),header.numEntradas+header.numSalidas,fileIn);
            if (leidos<(header.numEntradas+header.numSalidas))
            {
                printf("\nError 4: main() llamando a fread(), eof antes de numDatos?\n");
                return(0);
            }
            // convierte cada leido en double y lo guarda en buffer[0<numIn+numOut]
            for (j=0;j<(header.numEntradas+header.numSalidas);j++)
            {
                buffer[i][j]=(double)tmpFloatArr[j];
            }
        }
    }
    else //para datos double
    {

        for (i=0;i<numDatos;i++)
        {
            leidos=fread(buffer[i],header.numEntradas+header.numSalidas,sizeof(double),fileIn);
            if (leidos<(header.numEntradas+header.numSalidas))
            {
                printf("\nError 5: main() llamando a fread(), eof antes de numDatos?\n");
                return(0);
            }
        }
    }
    // cierra archivo
    fclose(fileIn);
    // calcula el vector de medias
    for(i=0;i<numDatos;i++)
        for (j=0;j<(header.numEntradas+header.numSalidas);j++)
            medias[j]+=buffer[i][j];
    for (j=0;j<(header.numEntradas+header.numSalidas);j++)
    {
        medias[j]/=numDatos;
    }
// :)    Si esto funciona bién me pido la mejor Hamburguesa Disponible   :)
    // calcula el s[0]=ponderacMultiOut()
    maxD=0;
    result=0;
    for (i=0;i<header.numEntradas;i++)
    {
        mUp=(metodoPonderac==0?0:1);
        mDn=(metodoPonderac==0?0:1);
        // calcula mUp= multip(correlac(index,Salidas)
        for (j=0;j<header.numSalidas;j++)
        {
            if (metodoPonderac==0)
                mUp=mUp+fabs(correlac(i, header.numEntradas+j, buffer, correlacM, medias, numDatos));
            else
                mUp=mUp*fabs(correlac(i, header.numEntradas+j, buffer, correlacM, medias, numDatos));
        }
        // calcula mDn= multip(correlac(index,entradasSeleccionadas)) es decir entre 0 y index
        for (j=0;j<header.numEntradas;j++)
            if (metodoPonderac==0)
                mDn=mDn+fabs(correlac(i, j, buffer, correlacM, medias, numDatos));
            else
                mDn=mDn*fabs(correlac(i, j, buffer, correlacM, medias, numDatos));
//TODO: probando: antes era fabs(mUp)/fabs(mDn) y ahora fabs(mUp*mUp)/fabs(mDn)
        if (mDn!=0)
        {
            printf("e%3d = %5.5E\n",i,fabs(mUp*mUp)/fabs(mDn));
            sumaTotal+=(fabs(mUp*mUp)/fabs(mDn));
            if ((fabs(mUp*mUp)/fabs(mDn))>maxD)
            {
                maxD=fabs(mUp*mUp)/fabs(mDn);

                result=i;
            }
        }

/*
       //probando
        if (mDn!=0)
        {
            printf("e%3d = %5.5E\n",i,fabs(mUp));
            sumaTotal+=(fabs(mUp));
            if ((fabs(mUp))>maxD)
            {
                maxD=fabs(mUp);

                result=i;
            }
        }

    }
*/
    // coloca en s[0] el resultado (lo intercambia con el valor anterior)
    s[0]=result;
    s[result]=0;
    valorPond[s[0]]=maxD;
    printf("\nSUM = %7.7f\n  0 = %3d (%7.7E)\n",sumaTotal,s[0],valorPond[s[0]]);
    sumaTotal=0;
    // para los seleccionados.;i++
    for (i=1;i<numEntradasSelec;i++)
    {
        //	max=0
        maxD=0;
        //	result=i
        result=i;
        valorPond[i]=0;
        //	para j=i;j<numEntradas;j++
        for (j=i;j<header.numEntradas;j++)
        {
            //  tmp=ponderacMultiOut(s[j])
            tmp=ponderacMultiOut(j, s,buffer, correlacM, medias, header, numDatos, metodoPonderac);

           //printf("e%d=%5.5E\n",s[j],fabs(tmp));
            //	si max<=fabs(tmp)
            if (fabs(tmp)>maxD)
            {
                //	max=fabs(tmp);
                maxD=fabs(tmp);
                //	result=s[j];
                result=j;
                //	valorPond[i]=tmp;
                valorPond[i]=tmp;
            }
        }
        //	si result!=i
        if (result!=i)
        {
            tmpI=s[i];
            s[i]=s[result];
            s[result]=tmpI;
        }
        else
             valorPond[i]=ponderacMultiOut(result, s,buffer, correlacM, medias, header, numDatos,metodoPonderac);

        //  imprime("%3d = %3d (%7.7f)",i,s[i],valorPond[i])
        sumaTotal+=(fabs(valorPond[i]));
        printf("%3d = %3d (%7.7E)\n",i,s[i],valorPond[i]);
    }
    printf("\nSUM_NuevaPond = %7.7f\n",sumaTotal);
    printf("Guardando archivo\n");
    // abre GTDv1Out para escritura.
    fileOut=fopen(fileNameOut,"wb");
    if (!fileOut)
    {
        printf("\nError 7: main() llamando a fopen()");
        return(0);
    }
    // compone nuevo encabezado GTDv1
    numEntradasAnt=header.numEntradas;
    header.numEntradas=numEntradasSelec;
    // escribe encabezado GTDv1 y verifica.
    leidos=fwrite(&header, sizeof(hdrGTDv1),1,fileOut);
    if (leidos<1)
    {
        printf("\nError 8: main() llamando a fwrite()");
        return(0);
    }
    // para i=0;i<numDatos;i++
    tmpDoubleArr = (double*) malloc((header.numEntradas+header.numSalidas)*sizeof(double));
    for(i=0;i<numDatos;i++)
    {
        // si tamRegistros==4
        if (header.tamRegistros==4)
        {
            // para j=0;j<numEntradasSelec;j++
            //compone el vector en tmpVector de floats
            for (j=0;j<numEntradasSelec;j++)
            {
                tmpFloatArr[j] = (float) buffer[i][s[j]];
            }
            for (j=0;j<header.numSalidas;j++)
            {
                tmpFloatArr[numEntradasSelec+j] = (float) buffer[i][numEntradasAnt+j];
            }
            //escribe el vector.
            leidos=fwrite(tmpFloatArr,sizeof(float),numEntradasSelec+header.numSalidas,fileOut);
            if (leidos<(numEntradasSelec+header.numSalidas))
            {
                printf("\nError 9: main() llamando a fwrite()");
                return(0);
            }
        }
        // sino escribe igual pero doubles.
        else
        {

            // para j=0;j<numEntradasSelec;j++
            //compone el vector en tmpVector de floats
            for (j=0;j<numEntradasSelec;j++)
            {
                tmpDoubleArr[j] = buffer[i][s[j]];
            }
            for (j=0;j<header.numSalidas;j++)
            {
                tmpDoubleArr[numEntradasSelec+j] = buffer[i][numEntradasAnt+j];
            }
            //escribe el vector.
            leidos=fwrite(tmpDoubleArr,sizeof(double),numEntradasSelec+header.numSalidas,fileOut);
            if (leidos<(numEntradasSelec+header.numSalidas))
            {
                printf("\nError 9: main() llamando a fwrite()");
                return(0);
            }
        }
    }
    // cierra fOut
    fclose(fileOut);
    printf("\n");
    return 0;
}


/*
Stay in one spot, another day of monotony
Has gotten me to the point, I'm like a snail
I've got to formulate a plot or I end up in jail or shot
Success is my only motherfucking option, failure's not
Mom, I love you, but this trailer's got to go
I cannot grow old in Salem's lot
So here I go, is my shot.
Feet fail me not cause maybe the only opportunity that I got

You can do anything you set your mind to, man
-Eminem
*/

