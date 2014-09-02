//ANN
//Evalua o entrena redes de hasta 65536 neuronas (requeriría 34.3GB de memoria, 10k neruronas requieren 800MB)
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#define MILLISEC 1000

//QUITAR: Dimensión de la matriz cuadrada


//Variables globales (globales debido a que se usan variables de más de 100MB y así se reduce uso de memoria).
 float version=0.013;
 float **matrizPesos ; //Matriz cuadrada de pesos float.(FALTA: mejor hacer punteros a punteros para evitar hacer multiplicación durante indexación)
 double **dMatrizPesos;//Matriz cuadrada de pesos double.(FALTA: mejor hacer punteros a punteros para evitar hacer multiplicación durante indexación)
 float *vectorValores ; //Vector de valores de tipo float.
 double *dVectorValores;//Vector de valores de tipo double.
 int useFloat=0; //por defecto usa matriz de pesos (y entradas) double 
 int trainingMode=0; //por defecto NO usa modo de entrenamiento (modo de evaluación de ANN en entradas)
 int tSigma=5; //Por defecto usa elliot unitario
 float fSigmaD=1; //Por defecto usa un D=1 para la función sigma Normal float(0)
 double dSigmaD=1; //Por defecto usa un D=1 para la función sigma Normal double(0)
 int sizeMatrizPesos=0; //tamaño de la matriz cuadrada de pesos (solo número de filas o número de neuronas)
 char* fileName; //nombre de archivo de pesos bin32/64 de entrada
 char* inputFileName; //nombre de archivo de entradas para evaluación/entrenamiento de ANN
 char* genomaFileName; // nombre de archivo de genoma
 int numSalidas=4; //número de neuronas en capa de salida.
 char* outputFileName; //nombre de archivo de salida
 char* outputMatrixFileName; //nombre de archivo de matriz de pesos entrenada de salida (solo en modo entrenamiento).
 long fileSize=0; //Tamaño del archivo de matriz de pesos de entrada en bytes.
 int neuronas=0; //para red FF (param -w<neuronas>) especifica el número total de neuronas.
 int capasOcultas=1; //para red FF (param -w<neuronas>)especifica el número de capas ocultas por defecto 1.
 int numEntradas=2; //número de neuronas de entrada, se colocan con el parametro -i<num>, pueden generarse desde genoma.
 FILE *f_in;//puntero a archivo de entradas.
 float threshold=0; //número usado como umbral de la función de activación con parámetro -r<num>.
 float dThreshold=0; //número usado como umbral de la función de activación con parámetro -r<num>.
 /******************************************/
/*       función delay										*/
/******************************************/

//Función delay, espera el número de ms especificado.
void delay(clock_t lMillisec){
	clock_t finish;
	finish = clock() + (lMillisec / MILLISEC * CLOCKS_PER_SEC);
	while (finish > clock());
}

/******************************************/
/*       función fSigma										*/
/******************************************/

//fSigna, retorna un float corespondiente a la función de activación seleccionada con param para una entrada X
float fSigma(float fX, int param, float fD){
	float y=0;
	//Param =	0 = sigma(0,1),		y = 1 / (1 + exp (- D * x))
	//		1 = sigma aprox (0,1)	y = 0.5 + x * (1 - abs(x) / 2), y=0 si x<=-1, y=1 si x>=1
	//		2 = tanh(-1,1),		y = 2 / (1 + exp(-2 * x)) - 1
	//		3 = gauss(-1,1),	y = exp(- x * x)
	//		4 = elliot,(-1,1)	y = x / (1 + |x|)
	//		5 = elliot (0,1 )	y = (x / 2) / (1 + |x|) + 0.5 
	if (param==0)
		return (1 / (1 + exp (- fD * fX)));
	if (param==1){
		if ((fX>-1)&&(fX<1))
			y = 0.5 + fX * (1 - (abs(fX) / 2));
		if (fX<=-1) y=0;
		if (fX>=1) y=1;
		return (y);
	}
	if (param==2)
		return (2 / (1 + exp(-2 * fX)) - 1);
	if (param==3)
		return (exp(- fX * fX));
	if (param==4)
		return (fX / (1 + abs(fX)));
	if (param==5)
		return ((fX / 2) / (1 + abs(fX)) + 0.5);
	return(0);
}

/******************************************/
/*       función dSigma										*/
/******************************************/

//dSigma, retorna un double correspondiente a la función de activación seleccionada para una entrada X
double dSigma(double dX, int param, double dD){
	double y=0;
	//Param =	0 = sigma(0,1),		y = 1 / (1 + exp (- D * x))
	//		1 = sigma aprox (0,1)	y = 0.5 + x * (1 - abs(x) / 2), y=0 si x<=-1, y=1 si x>=1
	//		2 = tanh(-1,1),		y = 2 / (1 + exp(-2 * x)) - 1
	//		3 = gauss(-1,1),	y = exp(- x * x)
	//		4 = elliot,(-1,1)	y = x / (1 + |x|)
	//		5 = elliot (0,1 )	y = (x / 2) / (1 + |x|) + 0.5 
	if (param==0)
		return (1 / (1 + exp (- dD * dX)));
	if (param==1){
		if ((dX>-1)&&(dX<1))
			y = 0.5 + dX * (1 - (abs(dX) / 2));
		if (dX<=-1) y=0;
		if (dX>=1) y=1;
		return (y);
	}
	if (param==2)
		return (2 / (1 + exp(-2 * dX)) - 1);
	if (param==3)
		return (exp(- dX * dX));
	if (param==4)
		return (dX / (1 + abs(dX)));
	if (param==5)
		return ((dX / 2) / (1 + abs(dX)) + 0.5);
	if (param>5)
		return (0);
	return(0);
}

/*  tipos de funciones de activación sacados de: http://www.dontveter.com/bpr/activate.html */

/**************************************/
/*       Función help									*/
/**************************************/

//Función help, imprime ayuda.
void help(float version){
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

/****************************************/
/*     función procParameters						*/
/****************************************/

//procesa los command-line parameters retorna 0 si hay error 1 si ok.
int procParameters(int argc, char *argv[]){
 
  // Si no hay parámetros muestra error.
	if (argc==1){
	printf("\nANN V%3.3f\nSe necesitan parámetros. \nUse: ann -h para obtener ayuda.\n",version);
		return(0);
	}
	//Verifica cada parámetro.
	while ((argc > 1) && (argv[1][0] == '-')) {
		switch (argv[1][1]) {
		//-f usa float
		case 'f':
			useFloat = 1; 
			break;
		//-s<numero> Selecciona el tipo de función de activación (mirar función fSigma o dSigma)
		case 's':
			tSigma = atoi(&argv[1][2]);
		break;
		//-d<numero> valor de D solo para  tSigma=0
		case 'd':
			fSigmaD = atof(&argv[1][2]);
			break;
		//-b<filename> nombre de archivo bin32/64 (matriz cuadrada de pesos)
		case 'b':
			fileName = &argv[1][2];
			break;
		//-t<iteraciones> modo de entrenamiento , especifica el número de iteraciones, necesita un archivo de entradas (opción -e) y un archivo de genoma (opción -g) produce un archivo de pesos de salida (opción -m) 
		case 't':
			trainingMode = atoi(&argv[1][2]);
		break;
		//-e<filename> archivo de entradas, debe contener vectores de 32/64bits (ver opción -f) del mismo tamaño que la matriz de pesos de entrada o el especificado en el genoma.
		case 'e':
			inputFileName = &argv[1][2];
			break;
/************  FALTA: FUNCION DE INICIALIZACION DE FEEDFORWARD   *****************************/
		//-w<neuronas> fuerza creación de matriz para feedforward de <neuronas> neuronas.
		case 'w':
			neuronas = atoi(&argv[1][2]);
			break;
		//-c<capas> especifica número de capas ocultas default=1.
		case 'c':
			capasOcultas = atoi(&argv[1][2]);
			break;
/************  FALTA: DEFINIR FORMATO DE GENOMA   *****************************/
		//-g<filename> archivo de genoma (formato por definir)
	 	case 'g':
			genomaFileName = &argv[1][2];
			break;
/************  FALTA:  DEFINIR FORMATO DE SALIDA  *****************************/
		//-o<filename> archivo de salida (en formato por definir.)
		case 'o':
			outputFileName = &argv[1][2];
			break;
		//-m<filename> archivo de salida con matriz de pesos. (en formato por definir.)
		case 'm':
			outputMatrixFileName = &argv[1][2];
			break;
		//-h imprime ayuda.
		case 'h':
			help(version); //función de ayuda
			return(0);
			break;
		//-n<numero> número de neuronas en la capa de salida. \n
		case 'n':
			numSalidas = atoi(&argv[1][2]);
			break;
		//-i<numero> número de neuronas en la capa de entrada. \n
		case 'i':
			numEntradas = atoi(&argv[1][2]);
			break;
		//-r<numero> umbral de activación. \n
		case 'r':
			threshold = atof(&argv[1][2]);
			dThreshold = atof(&argv[1][2]);
			break;

		default:
			printf("\nANN V%3.3f\nOpción inválida: %s\nUse: ann -h para obtener ayuda.\n", version,argv[1]);
			break;
		}
		//fin  switch
		++argv;
		--argc;
	}//Fin while

	//Comprueba si hay error en tSigmaD
	if ((tSigma!=0)&&(fSigmaD!=1)){
		printf("\nANN V%3.3f\nOpción inválida: -d , solo puede usarse con -s0\nUse: ann -h para obtener ayuda.\n", version);
		return(0);
	}
	//comprueba si se especificó el filename del archivo de matriz de pesos de entrada sin la opción -t
	if ((trainingMode==0)&&(fileName==NULL)){
		printf("\nANN V%3.3f\nFaltan parámetros, archivo de matriz de pesos de entrada con: -f<filename> \ncuando NO se usa modo de entrenamiento con -t<iteraciones>\nUse: ann -h para obtener ayuda.\n", version);
		return(0);
	}
	return(1);
}

/******************************************************/
/*          Función leerMatrix								        */
/******************************************************/
//leerMatrix: inicializa los valores de la matriz de pesos desde el archivo fileName en trainigMode=0
//retorna 0 si hay error, retorna 1 si ok
//usada por la función inicMatrix
/**********************FALTA COMPROBAR QUE LEE ARCHIVO DE PRUEBA EN MEMORIA CORRECTAMENTE*************/
int leerMatrix(void){
	FILE *f;	
	int i9=0;
	unsigned int readed=0;
	//abre el archivo 
	f = fopen(fileName, "rb");
	if (f){
		if (useFloat==1){
			for (i9=0;i9<sizeMatrizPesos;i9++)
				readed=fread(matrizPesos[i9], sizeMatrizPesos*4, 1, f);
			fclose(f);
		}
		else{
			for (i9=0;i9<sizeMatrizPesos;i9++)
				readed=fread(dMatrizPesos[i9], sizeMatrizPesos*8, 1, f);
			fclose(f);
		}
	}
	else{
		printf("\nError abriendo el archivo: %s\n",fileName);
		return(0);
	}
	return(1);
}

/******************************************************/
/*          Función inicializarFF								        */
/******************************************************/
//Inicializa en 1 los pesos de las conexiones para formar una red FF teniendo en cuenta el
//número de neuronas, el número de neuronas de entrada y el de salida.
void inicializarFF(){
	int i3=0;
	int j3=0;
	int k3=0;
	int in=numEntradas;
	int o=numSalidas;
	int n=sizeMatrizPesos;
	int c=capasOcultas;
	int pc=div(n-in-o,c).quot;
	int r=div(n-in-o,c).rem;
	//Coloca todos los elementos de la matriz de pasos  en 0.
	//obtener numCapa_1 y numCapa_n
	//se incializa contador i en numEntradas-1
	for (j3=0;j3<n;j3++){
		for (i3=0;i3<n;i3++){
					if (useFloat==1){
						matrizPesos[j3*sizeMatrizPesos+i3]=0;
					}
					else{
						dMatrizPesos[j3*sizeMatrizPesos+i3]=0;
					}
		//FALTA
		}
	}
	/*
	*primeras conexiones a capa grande
para j desde 0 hasta in-1{
	para i desde in hasta in+(pc+r)-1{
		matr[j][i]=1
	}
}
*/
for (j3=0;j3<=in-1;j3++){
	for (i3=in;i3<=(in+(pc+r)-1);i3++){
		if (useFloat==1){
			matrizPesos[j3][i3]=1;
		}
		else{
			dMatrizPesos[j3][i3]=1;
		}
	}
}


/*

*conexiones para capa grande
para j=in hasta in+pc+r-1{
	para i=in+(pc+r) hasta in+(pc+r)+pc-1{
		matr[j][i]=1
	}
}
*/
for (j3=in;j3<=in+pc+r-1;j3++){
	for (i3=in+pc+r;i3<=(in+(pc+r)+pc-1);i3++){
		if (useFloat==1){
			matrizPesos[j3][i3]=1;
		}
		else{
			dMatrizPesos[j3][i3]=1;
		}
	}
}

/*
*conexiones de capas iguales:
para k desde 0 hasta q-3{
	para j desde in+(pc+r)+(k*pc) hasta in+(pc+r)+((k+1)*pc)-1{
		para i desde in+pc+r+pc+(k*pc) hasta in+pc+r+pc+((k+1)*pc)-1{
			matr[j][i]=1
		}
	}
}
*/
for (k3=0;k3<=c-3;k3++){
	for (j3=in+(pc+r)+(k3*pc);j3<=(in+(pc+r)+((k3+1)*pc)-1);j3++){
		for (i3=in+pc+r+pc+(k3*pc);i3<=(in+pc+r+pc+((k3+1)*pc)-1);i3++){
			if (useFloat==1){
				matrizPesos[j3][i3]=1;
			}
			else{
				dMatrizPesos[j3][i3]=1;
			}
		}
	}
}

/*

*conexiones para capa de salida:
k=q-2
para j=in+(pc+r)+(k*pc) hasta in+(pc+r)+((k+1)*pc)-1{
	para i=in+pc+r+pc+((k)*pc) hasta in+pc+r+pc+((k)*pc)+o-1{
		matr[j][i]=1
	}
}
	*/
for (j3=in+(pc+r)+(k3*pc);j3<=(in+(pc+r)+((k3+1)*pc)-1);j3++){
	for (i3=in+pc+r+pc+(k3*pc);i3<=(in+pc+r+pc+((k3)*pc)+o-1);i3++){
		if (useFloat==1){
			matrizPesos[j3][i3]=1;
		}
		else{
			dMatrizPesos[j3][i3]=1;
		}
	}
}	

}


/******************************************************/
/*          Función inicMatriz								        */
/******************************************************/
//verifica tamaño de archivo correcto e inicializa valores de matriz de pesos
//dependiendo de si ANN está en modo train o no (desde archivo).
//retorna 0 si hay error, 1 si ok.
int inicMatriz(void){
 int i7=0;
 FILE *f;
/*************FALTA seleccionar si generar matriz o leer desde archivo o modo train (debería ser una función.)****/
 if (trainingMode==0){
	 f = fopen(fileName, "rb");
	 if (f){
		fseek(f,0,SEEK_END);
		fileSize = ftell(f);
		fclose(f);
		if (useFloat==1){
			sizeMatrizPesos=(int)(sqrt(ldiv(fileSize,4).quot)+0.5); //para float tam=round(sqrt(filesize/4))
			/****Verifica que el tamaño del archivo bin32 sea exacto.*****/   		
	   		if ((sizeMatrizPesos*sizeMatrizPesos*4)!=fileSize){
					printf("\nNúmero incorrecto de elementos de 4bytes en archivo bin32 (float) en: %s\nEl archivo debe contener valores para una matriz cuadrada\n",fileName);
					return(0);
	   		}
	   	}
	   	else{
	   		sizeMatrizPesos=(int)(sqrt(ldiv(fileSize,8).quot)+0.5); //para double tam=round(sqrt(filesize/8))
	   		//Verifica que el tamaño del archivo bin64 sea exacto.
	   		if ((sizeMatrizPesos*sizeMatrizPesos*8)!=fileSize){
					printf("\nNúmero incorrecto de elementos de 8bytes en archivo bin64 (double) en: %s\nEl archivo debe contener valores para una matriz cuadrada\n",fileName);
					return(0);
	   		}
	  	}
	 }//fin if(f)
	 else{
		printf("\nError abriendo el archivo: %s\n",fileName);
		return(0);
	 }//fin else ->if(f)
 }//fin if(trainingMode==0)
 else{
 		//si se está en modo de FF (-w<neuronas>) se hace el tamaño de la matriz = <neuronas>

 		if (neuronas>0)
 	  	sizeMatrizPesos=neuronas;
 }//fin else->if(trainingMode==0)	
	if (useFloat==1){
		//  Creación de matriz de pesos float 
		matrizPesos=(float**) malloc (sizeMatrizPesos*sizeof(float *)); //Ubica memoria para elementos
		for (i7=0; i7<sizeMatrizPesos; i7++)
			matrizPesos[i7] = (float *)malloc (sizeMatrizPesos*sizeof(float));
	}//fin If(useFloat==1)
	else{
		/**********  Creación de matriz de pesos double  ***********/
		dSigmaD=(double)fSigmaD;
	 	dMatrizPesos=(double**) malloc (sizeMatrizPesos*sizeof(double *)); //Ubica memoria para elementos
		for (i7=0; i7<sizeMatrizPesos; i7++)
			dMatrizPesos[i7] = (double *) malloc (sizeMatrizPesos*sizeof(double));
	}//fin else ->if(useFloat==1)

  //Si trainingMode=0, lee desde archivo de matriz de pesos de entrada.
	//lee los valores de matrix desde el archivo fileName si hay error retorna 0.

	if (neuronas>0){
		inicializarFF();

	}
	if (trainingMode==0){
		if (!leerMatrix())
				return(0);
		}
		else{
			// si se usó FF (-w<neuronas>) se llama a inicializarFF
			/********** FALTA: FUNCION INICIALIZAR FF **********/
			inicializarFF();
		}
	return(1);
}

/******************************************************/
/*          Función freeMatrix								        */
/******************************************************/
//Libera la memoria usada para almacenar la matriz de pesos.
void freeMatrix(void){
	int i8=0;
	if (useFloat==1){
		for (i8=0; i8<sizeMatrizPesos; i8++)
   			free (matrizPesos[i8]);
		free(matrizPesos);
	}   
	else{
		for (i8=0; i8<sizeMatrizPesos; i8++)
   			free (dMatrizPesos[i8]);
		free(dMatrizPesos);
	}	
}

/******************************************************/
/*          Función freeVectorValores					        */
/******************************************************/
//Libera la memoria usada para almacenar la matriz de pesos.
void freeVectorValores(void){
	if (useFloat==1){
		free(vectorValores);
	}   
	else{
		free(dVectorValores);
	}	
}

/********************************************************/
/*          Función inicVectorValores	        	*/
/********************************************************/
//Crea el vector de valores del mismo tamaño que el num de neuronas (sizeMatrizPesos )
void inicVectorValores(void){
	int i4=0;
	if (useFloat==1){
		//  Creación de matriz de pesos float 
		vectorValores=(float*) malloc (sizeMatrizPesos*sizeof(float)); //Ubica memoria para elementos
		//coloca valores iniciales en 0
		for (i4=0;i4<sizeMatrizPesos;i4++)
			vectorValores[i4]=0;
	}//fin If(useFloat==1)
	else{
		/**********  Creación de matriz de pesos double  ***********/
		dSigmaD=(double)fSigmaD;
	 	dVectorValores=(double*) malloc (sizeMatrizPesos*sizeof(double)); //Ubica memoria para elementos
	 	//coloca valores iniciales en 0
	 	for (i4=0;i4<sizeMatrizPesos;i4++)
	 		dVectorValores[i4]=0;
	}//fin else ->if(useFloat==1)
}


/********************************************************/
/*          Función leerEntrada		        	*/
/********************************************************/
//Lee los valores de entrada correspondientes al  número especificado por el parámetro -i<num> o por el Genoma (el genoma tiene prioridad FALTA)
//retorna 1 si la lectura fué correcta, 2 si se llegó al EOF y 0 si el archivo no existe o no se puede abrir.
//coloca los valores leidos en los primeros elemantos del vector de valores.
int leerEntrada(){

	
	unsigned int readed_in=0;
	int i5=0;
	//abre el archivo 
	
	f_in = fopen(inputFileName, "rb");
	if (f_in){
		if (useFloat==1){
			for (i5=0;i5<numEntradas;i5++){
				readed_in=fread(&vectorValores[i5], sizeof(float), 1, f_in);
				if (readed_in==0)
					return(2);
			}
		}
		else{
			for (i5=0;i5<numEntradas;i5++){
				readed_in=fread(&dVectorValores[i5], sizeof(double), 1, f_in);
				if (readed_in==0)
					return(2);
			}
		}
	}
	else{
		printf("\nError abriendo el archivo: %s\n",inputFileName);
		return(0);
	}
	return(1);
}

//FALTA:OPTIMIZAR MATRIZ DE PESOS CON DOBLE PUNTERO PARA EVITAR MULTIPLICACIÓN en CADA REFERENCIA
//Calcula el resto de los elementos del vector de valores (Sip, incluyendo salidas :)  )
void calcularElResto(void){
	int i6=0;
	int j6=0;
	float facum=0;
	double dacum=0;
	if (useFloat==1){
		for (i6=numEntradas;i6<sizeMatrizPesos;i6++){
			facum=0;
			for (j6=0;j6<sizeMatrizPesos;j6++){
				facum=facum+(matrizPesos[j6][i6]*(vectorValores[j6]));
			}
			//fSigma(float fX, int param, float fD)
			vectorValores[i6]=fSigma((-threshold+facum),tSigma,fSigmaD);
		}
	}
	else{
		for (i6=numEntradas;i6<sizeMatrizPesos;i6++){
			dacum=0;
			for (j6=0;j6<sizeMatrizPesos;j6++){
				dacum=dacum+(dMatrizPesos[j6][i6]*(dVectorValores[j6]));
			}
			//fSigma(float fX, int param, float fD)
			dVectorValores[i6]=dSigma((-dThreshold+dacum),tSigma,dSigmaD);
		}
	}
}

/*****************************************************************************************************************************/

//main. ANN

/*****************************************************************************************************************************/

int main(int argc, char *argv[]){
	unsigned long int i,j;
	//verifica y asigna los parámetros de linea de comandos, si hay algún error, sale.
	if (!procParameters(argc,argv))
	return(0);
	//Verifica tamaño de matriz de pesos y inicializa sus elementos si trainingMode=1
	if (!inicMatriz())
		return(0);
	//Crea a inicializa el vector de valores.	
	inicVectorValores();
	
	//Lee un número de valores de entrada = numEntradas y retorna 2 si llegó a EOF, 1 so OK y 0 si hay error abriendo el archivo.
	//No cierra el archivo.
	if (!leerEntrada())
		return(0);
	
	//Calcula el resto de los elementos del vector de valores
	calcularElResto();
	
	printf("\nfileName = %s, fileSize=%li bytes, matrixSize=%i, useFloat=%d, tSigma=%d, tSigmaD=%3.3f \n\n",fileName,fileSize,sizeMatrizPesos,useFloat,tSigma,fSigmaD);
	if (useFloat==1){
		for (i=0;i<sizeMatrizPesos;i++){
			for (j=0;j<sizeMatrizPesos;j++)
				printf ("%1.1f, ",matrizPesos[i][j]);
				printf("-> %1.1f\n",vectorValores[i]);
		}
	}
	else{
		for (i=0;i<sizeMatrizPesos;i++){
			for (j=0;j<sizeMatrizPesos;j++)
				printf ("%1.1f, ",dMatrizPesos[i][j]);
				printf("  -> %1.1f\n",dVectorValores[i]);
		}
	}
	//Libera la memoria usada por el vector de valores y cierra el archivo de entradas.
	fclose(f_in);
	freeVectorValores();
	freeMatrix();	
	return(0);
} //Fin main  :)

//Teléfono Jairo: 321 6429465

