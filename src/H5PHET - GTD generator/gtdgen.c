/** H5PHET - GTD Generator.
	Crea archivos de datos de entrenamiento en formato GTD para realizar pruebas con H5PHET - NEAT.
	Parte de H5PHET.
	Por Harvey Bastidas.
*/

#include <stdio.h> //para printf
#include <stdlib.h> //para malloc
#include <time.h> //para time
#include <math.h> //para floor
#include <errno.h>

// estructura del encabezado de GTDv1
typedef struct
{
 	char fileID[3];
 	char version;
 	unsigned usarSigned;
 	unsigned tamRegistros;
 	unsigned numEntradas;
 	unsigned numSalidas;
} hdrGTDv1;

/******************************************* inicRandomSeed ********************************************/
void inicRandomSeed() //inicializa el random seed
{
	int i,k;
	int j=1;
	// inicializa random seed
	// (genera un número entre 0 y 1e3)
    i=(unsigned int) time(NULL);
    j=j*(unsigned int)(1+floor((i-10*floor(i/10))+0.5));
    j=j*(unsigned int)(floor((i-1000*floor(i/1000))+0.5));
    j=j+(unsigned int)floor((i-10*floor(i/10))+0.5)+1;
    j=j+(unsigned int)floor((i-1000*floor(i/1000))+0.5);
	// inicializa el pseudo-random seed con la hora actual.
    srand(j);
    k=1;
    for (i=0;i<j;i++)
        k=rand();
    // re-Inicializa el pseudo-random seed con el último rand obtenido.
    srand(k);
}

/******************************************* escribirEncabezadoGTDv1 ********************************************/
size_t escribirEncabezadoGTDv1(FILE* fileOut,unsigned char usarSigned,unsigned char tamRegistros,unsigned numEntradas,unsigned numSalidas)
// escribe el encabezado de Generic Training Data:
// 		3 unsigned char = GTD, char version=1 , unsigned char usarSigned, unsigned char tamRegistros, unsigned numEntradas, unsignedNumSalidas.
// Parámetros:
// 	    fileOut = handler de archivo de salida previamente abierto para escritura.
//		usarSigned = 0=NO, 1=SI
//		tamRegistros = 4=float, 8=double
//		numEntradas = número de entradas
//		numSalidas = número de salidas
// Retorna 0 si hay error, != 0 si OK
{
	hdrGTDv1 header;
	int headerLen;
	header.fileID[0]='G';
	header.fileID[1]='T';
	header.fileID[2]='D';
	header.version=1;
	header.usarSigned=usarSigned;
	header.tamRegistros=tamRegistros;
	header.numEntradas=numEntradas;
	header.numSalidas=numSalidas;
	// escribe GTD y la versión
	headerLen = fwrite(&header, sizeof(hdrGTDv1), 1, fileOut);
	if (headerLen !=1)
	{
		printf("\nError en funcion escribirEncabezadoGTDv1() tamaño de encabezado escrito incorrecto (%i).\n",headerLen);
		return(0);
	}
	return(1);
}

/******************************************* probarGTD() ********************************************/
void probarGTD(char* filenameGTDv1)
// lee el archivo GTD especificado e imprime los 20 primeros registros, retona 0 si error.
{
	hdrGTDv1 header;
	float listaFloats[10160]; //SI se usa un valor de 1M o más se genera una excepción en windows
	int tmp;
//	double listaDoubles[10000000];  //máximo 10M datos para esta prueba
	FILE* fOut;
	int i;
	// abre el archivo
    if ((fOut=fopen(filenameGTDv1,"rb"))==NULL)
    {
		printf("\nError abriendo el archivo de entradas %s para escritura\n",filenameGTDv1);
		return;
    }
    // lee encabezado
	if (fread(&header, sizeof(hdrGTDv1), 1, fOut)!=1)
	{
		printf("Error leyendo encabezado GTD de %s\n",filenameGTDv1);
		return;
	}
    // falta : comparar 3 chars iniciales con GTD;
	// si versión != 1 notifica error.
    if (header.version!=1)
    {
    	printf("Error: versión de GTD diferente de 1\n");
    	return;
    }
	// lee 20 primeros datos.
	if (header.tamRegistros==4) //si es 4=float
	{
		if ((header.numEntradas<255)&&(header.numSalidas<255))//apra evitar Denial of Service por datos corruptos de numEntradas y numSalidas.
		{
			if ((tmp=fread(listaFloats, (header.numEntradas+header.numSalidas)*sizeof(float), 20, fOut))!=20)
			{
				printf("Error leyendo valores de  GTD  %s  n´mero de registros <20 \n",filenameGTDv1);
				return;
			}
			else
			{
				// imprime los 20 registros
				for(i=0;i<20;i++)
					if ((header.numEntradas<255)&&(header.numSalidas<255))//apra evitar Denial of Service por datos corruptos de numEntradas y numSalidas.
						printf("dato %d :  %3.3f XOR %3.3f = %3.3f\n",i,listaFloats[i*(header.numEntradas+header.numSalidas)],listaFloats[(i*(header.numEntradas+header.numSalidas))+1],listaFloats[(i*(header.numEntradas+header.numSalidas))+2]);
			}
		}
	}
	// cierra el archivo GTD
    if(fclose(fOut)!=0)
    {
		printf("Error escribiendo el archivo de salidas llamando a fclose\n");
		return;
    }
}

/******************************************* genListaOuts() ********************************************/
void genListaOuts(float *listaOuts, int numPares, int usarSigned, int numEntradas, int numSalidas)
// genera los valores de entrenamiento y los coloca en listaOuts
{
	int i;
	int j=1;
	int k=0;
	// genera aleatoriamente dos entradas 0 o 1 y los coloca como float en el arreglo de entradas
	for(i=0;i<numPares;i++)
	{
		if ((((float)rand())/RAND_MAX)<=0.5)
			j = 0;
		else
			j = 1;
		if ((((float)rand())/RAND_MAX)<=0.5)
			k = 0;
		else
			k = 1;
		// llena el arreglo listaOuts con los valores de j, k y su XOR
		if ((((numEntradas+numSalidas)*i)+2)<100000)
		{
			listaOuts[((numEntradas+numSalidas)*i)] = j==1? 1 : usarSigned==1? -1 : 0;
			listaOuts[((numEntradas+numSalidas)*i)+1] = k==1? 1 : usarSigned==1? -1 : 0;
			listaOuts[((numEntradas+numSalidas)*i)+2] = (k^j)==1? 1 : usarSigned==1? -1 : 0;
		}
	}
}

/******************************************* main() ********************************************/
int main()
{
	long int numPares=1000; //número de pares aleatorios a genera
	int usarSigned=1; //1= generar valores aleatorios con signo, 0 = sin signo.
	int numEntradas=2;
	int numSalidas=1;
	int tamRegistros=4;
	FILE *fOut=NULL;
	char *filenameGTDv1;
	float *listaOuts; // para almacenar los valores en float de las entradas y su resultado.
	filenameGTDv1 = "d:\\h5phet\\trdata.gtd\0";
	// reserva memoria para los  arreglo listaOuts
	if ((listaOuts=(float *)malloc((numSalidas+numEntradas)*numPares*sizeof(float)))==NULL)
	{
		printf("\nError llamando a malloc para listaIns\n");
		return(0);
	}
	// inicializa random seed.
	inicRandomSeed();
	// llena el buffer de datos de entrenamiento listaOuts
	printf("Generado datos de entrenamiento.\n");
	genListaOuts(listaOuts, numPares, usarSigned, numEntradas, numSalidas);
	// abre el archivo de salida GTD
	printf("Abriendo archivo GTD.\n");
	//if ((fOut=fopen("d:/h5phet/test_in.txt","w+"))==NULL){
    if ((fOut=fopen(filenameGTDv1,"wb"))==NULL)
    {
		printf("\nError abriendo el archivo %s para escritura\n",filenameGTDv1);
		free(listaOuts);
		return(0);
    }
	// escribe el encabezado de GTD: tipoArch (3 unsigned char GTD), char version, unsigned char usarSigned, unsigned char tamRegistros, unsigned numEntradas, unsignedNumSalidas.
	printf("Escribiendo encabezado GTD.\n");
    if (escribirEncabezadoGTDv1(fOut, usarSigned, tamRegistros, numEntradas, numSalidas)==0) //fileHandle(debe estar abierto para ewscritura),usarSigned,tamRegistros, numEntradas, numSalidas
	{
		printf("Error escribiendo el encabezado en el archivo GTD\n");
		free(listaOuts);
		return(0);
	}
	// escribe datos en archivo GTD
	printf("Escribiendo datos\n");
	if (fwrite(listaOuts, tamRegistros, (numEntradas+numSalidas)*numPares, fOut) != ((numEntradas+numSalidas)*numPares))
	{
		printf("Error escribiendo archivo GTD llamando a fwrite()\n");
		free(listaOuts);
		return(0);
	}
	printf("Cerrando archivo.\n");
	// cierra el archivo de salida
    if(fclose(fOut)!=0)
    {
		printf("Error escribiendo el archivo GTD llamando a fclose()\n");
		free(listaOuts);
		return(0);
    }
    printf("Probando funcionamiento\n");
    // Lee 20 datos desde el archivo generado y los imprime en pantalla.
    probarGTD(filenameGTDv1);
	// libera memoria de las listas usadas.
    free(listaOuts);
    return 0;
}
