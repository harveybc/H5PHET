/** Funciones para menejo de especies - C file
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
#ifndef GENOMA_H_INCLUDED
#include "genoma.h"
#define GENOMA_H_INCLUDED
#endif
#include "especie.h"


void imprimirListasInnov(TConfig* conf)
{
    unsigned i,j;
    //Imprime innovaciones de nodo formato: In
    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nImprime la lista de innovaciones de nodos y conex\nFormato: **i<in>(<contInnov>),o<out>=(<innovnum>,<nodoout>)... ");
    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nLista de innov Nodos: contInnovNodo=%u<br>",conf->contInnovNodo);
    for (i=0; i<conf->contInnovNodo; i++)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"**i%u(%u),",i,conf->listaInnovNodo[i].numOut);
        for (j=0; j<conf->listaInnovNodo[i].numOut; j++)
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"o%u=(%u,%u) ",j,conf->listaInnovNodo[i].nodoOut[j].innovNum,conf->listaInnovNodo[i].nodoOut[j].nodoOut);
    }
    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nLista de innov Conexiones: contInnovConex=%u<br>",conf->contInnovCon);
    for (i=0; i<conf->contInnovNodo; i++)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"**i%u(%u),",i,conf->listaInnovCon[i].numOut);
        for (j=0; j<conf->listaInnovCon[i].numOut; j++){
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"o%u=(%u,%u) ",j,conf->listaInnovCon[i].nodoOut[j].innovNum,conf->listaInnovCon[i].nodoOut[j].nodoOut);
        }
    }

}

unsigned verificarListasInnov(TConfig* conf)
{
    //retorna 0 si hay valores inválidos en la lista de innovación y imprime el maxInnovNum
    //TODO: quitar cuando no haya errores de innovnums
    unsigned i,j;
//    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>miN=%i,miC=%i,",conf->contInnovNodo,conf->contInnovCon);
    for (i=0; i<conf->contInnovNodo; i++)  //para la lista de innovaciones de nodos
    {
        for (j=0; j<conf->listaInnovNodo[i].numOut; j++)  //busca límites inferiores para innovnum y nodoOut para cada innovación
        {
            if (conf->listaInnovNodo[i].nodoOut[j].innovNum>=conf->contInnovNodo)  //verifica maximo innovnumNodo
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 51 en verificarListasInnov(), liNodo[%u].no[%u]=%u > ciN=%u",i,j,conf->listaInnovNodo[i].nodoOut[j].innovNum,conf->contInnovNodo);
                return(0);
            }
            if (conf->listaInnovNodo[i].nodoOut[j].innovNum<0)  //verifica que el innovnumnodo sea >=0
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 52 en verificarListasInnov(), liNodo[%u].no[%u].iN=%u < 0",i,j,conf->listaInnovNodo[i].nodoOut[j].innovNum);
                return(0);
            }
            if (conf->listaInnovNodo[i].nodoOut[j].nodoOut>=conf->contInnovNodo)  //verifica que el innovnumnodo sea >=0
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 53 en verificarListasInnov(), liNodo[%u].no[%u].nO=%u > ciN=%u",i,j,conf->listaInnovNodo[i].nodoOut[j].nodoOut,conf->contInnovNodo);
                return(0);
            }
            if (conf->listaInnovNodo[i].nodoOut[j].nodoOut<0)  //verifica que el innovnumnodo sea >=0
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 54 en verificarListasInnov(), liNodo[%u].no[%u].nO=%u < 0",i,j,conf->listaInnovNodo[i].nodoOut[j].nodoOut);
                return(0);
            }
        }
    }
    for (i=0; i<conf->contInnovNodo; i++)  //para la lista de innovaciones de conex
    {
        for (j=0; j<conf->listaInnovCon[i].numOut; j++)  //busca límites inferiores para innovnum y nodoOut para cada innovación
        {
            if (conf->listaInnovCon[i].nodoOut[j].innovNum>=conf->contInnovCon)  //verifica maximo innovnumNodo
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 55 en verificarListasInnov(), liCon[%u].no[%u]=%u > ciC=%u",i,j,conf->listaInnovCon[i].nodoOut[j].innovNum,conf->contInnovCon);
                return(0);
            }
            if (conf->listaInnovCon[i].nodoOut[j].innovNum<0)  //verifica que el innovnumnodo sea >=0
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 56 en verificarListasInnov(), liCon[%u].no[%u].iN=%u < 0",i,j,conf->listaInnovCon[i].nodoOut[j].innovNum);
                return(0);
            }
            if (conf->listaInnovCon[i].nodoOut[j].nodoOut>=conf->contInnovNodo)  //verifica que el innovnumnodo sea >=0
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 57 en verificarListasInnov(), liCon[%u].no[%u].nO=%u > ciN=%u",i,j,conf->listaInnovCon[i].nodoOut[j].nodoOut,conf->contInnovNodo);
                return(0);
            }
            if (conf->listaInnovCon[i].nodoOut[j].nodoOut<0)  //verifica que el innovnumnodo sea >=0
            {
                fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 58 en verificarListasInnov(), liCon[%u].no[%u].nO=%u < 0",i,j,conf->listaInnovCon[i].nodoOut[j].nodoOut);
                return(0);
            }
        }
    }
    return(1);
}

float calcularDist(unsigned indexpob1,unsigned indexpob2,float c1, float c2, float c3,unsigned eG_t, TConfig* conf)  //OPTIMIZADA
{
//Retorna la distancia entre dos genomas. uusando la formula D=((((C1)*E)/n)+(((C2)*D)/n)+C3*W)
//			Donde: 	C1,C2,C3 = constantes de proporcionalidad (mirar ejemplos para valores iniciales.)
//					n = Número de genes en el genoma más grande.
//					E = Número de genes de exceso
//					D = Número de genes disjounsigned  int
//					W = Promedio de diferencia de pesos entre genes correspondientes. Puede ser seteado a 1 para valores que no sean excesivamente grandes.
//Parámetros: 	indexpob1 , indexpob2, = index de los genomas en la conf->pob.
//				c1 = constante de proporcionalidad para Excess genes
//				c2 = constante de proporcionalidad para Disjounsigned genes
//				c3 = constante de proporcionalidad para el promedio de diferencias de pesos.
//				eG_t = threshold para considerar número de genes "excesivamente" grandes (y hacer n=1)
    ////TODO: Quitar inicializaciones innecesarias en todas las funciones.
    unsigned i=0;
    float disjointC=0;
    float excessC=0;
    float wAver=0;
    unsigned tmpInnov1=0;
    unsigned tmpInnov2=0;
    unsigned indexMenorInnovNumConex=indexpob1;
    unsigned indexMayorInnovNumConex=indexpob2;
    unsigned cont=0;
    unsigned n=0;
    unsigned tmpMaxInnovNumConex; //usado para acelerar calculos
    //determinar el mayor y el menor (maxInnovNum)
    if (conf->pob[indexMenorInnovNumConex].maxInnovNumConex>conf->pob[indexMayorInnovNumConex].maxInnovNumConex)
    {
        indexMayorInnovNumConex=indexpob1;
        indexMenorInnovNumConex=indexpob2;
    }
    // excess = buscar todos los innovNums (desde el menor conf->maxInnovNumConex hasta el mayor) en el genoma que tiene el mayor.
    if (conf->pob[indexMenorInnovNumConex].maxInnovNumConex<conf->pob[indexMayorInnovNumConex].maxInnovNumConex)
    {
        tmpMaxInnovNumConex=conf->pob[indexMayorInnovNumConex].maxInnovNumConex; // para acelerar evaluación en for.
        for (i=conf->pob[indexMenorInnovNumConex].maxInnovNumConex+1; i<=conf->pob[indexMayorInnovNumConex].maxInnovNumConex; i++)
            if (buscarInnovConex(indexMayorInnovNumConex,i,conf)!=UINT_MAX) excessC++;
    }
    // disjounsigned  int= contar disjounsigned  ints hasta el menor maxInnovNum entre los dos genes de conexiones
    tmpMaxInnovNumConex=conf->pob[indexMenorInnovNumConex].maxInnovNumConex; // para acelerar evaluación en for.
    for (i=0; i<=tmpMaxInnovNumConex; i++) //MAL
    {
        tmpInnov1=buscarInnovConex(indexpob1,i,conf);
        tmpInnov2=buscarInnovConex(indexpob2,i,conf);
        if(tmpInnov1 == UINT_MAX)
            if (tmpInnov2 != UINT_MAX ) disjointC++;
        if(tmpInnov2 == UINT_MAX)
            if (tmpInnov1 != UINT_MAX ) disjointC++;
        if(tmpInnov1!=UINT_MAX) //obtiene promedio de diferencias
            if (tmpInnov2!=UINT_MAX)
            {
                cont++;
                wAver=wAver+fabs(conf->pob[indexpob1].conex[tmpInnov1].peso-conf->pob[indexpob2].conex[tmpInnov2].peso);
            }
    }
    // TODO: faltan optimizar todos los for del programa para que evalúen una variable local en lugar de un valor referenciado.
    // calcula wAver
    if (cont>0)
        //wAver=floor(wAver/(float)cont); CORREGIDO
        wAver=wAver/(float)cont;
    // calcula n=mayor número de genes entre los dos genomas
    if (conf->pob[indexpob1].totalConexiones>conf->pob[indexpob2].totalConexiones) n=conf->pob[indexpob1].totalConexiones;
    else n=conf->pob[indexpob2].totalConexiones;
    // ((((C1)*E)/n)+(((C2)*D)/n)+C3*W)
    if(n>0) return(((((c1)*excessC)+((c2)*disjointC))/(float)n)+c3*wAver);
    // si error o los genomas son iguales retorna 0;
    return(0);
}

unsigned especieMinDist(unsigned indexpob,float c1, float c2, float c3,unsigned eG_t, TConfig* conf)  //OPTIMIZADA
{
// Retorna la especie del genoma al que se tiene la mínima distancia entre los conf->representantes
// Parámetros: indexpob = indice del genoma para el cual se busca la mindist.
//				c1 = constante de proporcionalidad para Excess genes
//				c2 = constante de proporcionalidad para Disjounsigned genes
//				c3 = constante de proporcionalidad para el promedio de diferencias de pesos.
//				eG_t = threshold para considerar número de genes "excesivamente" grandes (y hacer n=1)
    unsigned i=0;
    float minDist=99999999;
    float tmpDist=0;
    unsigned foundIndex=0;
    for (i=0; i<conf->numEspecies; i++)
    {
        tmpDist=calcularDist(indexpob,conf->representantes[i],c1,c2,c3,eG_t,conf);
        if (tmpDist<minDist)
        {
            minDist=tmpDist;
            foundIndex=i;
        }
    }
    return (foundIndex);
}

float distEspecieCercana(unsigned indexpob,float c1, float c2, float c3,unsigned eG_t, TConfig* conf)  //OPTIMIZADA
{
//retorna la distancia a la especie más cercana diferente a la del indexpob
//Parámetros: indexpob = indice del genoma para el cual se busca la mindist.
//				c1 = constante de proporcionalidad para Excess genes
//				c2 = constante de proporcionalidad para Disjounsigned genes
//				c3 = constante de proporcionalidad para el promedio de diferencias de pesos.
//				eG_t = threshold para considerar número de genes "excesivamente" grandes (y hacer n=1)
    unsigned i=0;
    float minDist=99999999; // retorna este valor si hay una sola especie
    float tmpDist=0;
    for (i=0; i<conf->numEspecies; i++)
    {
        if (conf->pob[indexpob].especie!=i){
            tmpDist=calcularDist(conf->representantes[conf->pob[indexpob].especie],conf->representantes[i],c1,c2,c3,eG_t,conf);
            if (tmpDist<minDist)
            {
                minDist=tmpDist;
            }
        }
    }
    return (minDist);
}

unsigned asignarEspecie(unsigned indexpob, float espThreshold,float c1, float c2, float c3,unsigned eG_t, TConfig* conf)  //OPTIMIZADA
{
// Usando especieMinDist obtiene la especie más compatible, y la compara su distancia con el threeshold, si es menor asigna la especie, si es mayor, crea
// una nueva y la asigna al genoma en cuestión
// retorna el número de la especie asignada., -1 si hay error
// parámetros:	indexpob	= indice del genoma al que se asignará la especie
//				threeshold	= límite(inferior) de distancia para pertenecier a una especie (//TODO: Variación dinámica de este parám)
    unsigned especieCercana=0;
    unsigned i=0;
    float minDist=99999999;
    float tmpDist=0;
    // busca especie más cercana entre los representantes (como Funcion especieMinDist)
    for (i=0; i<conf->numEspecies; i++)
    {
        tmpDist=calcularDist(indexpob,conf->representantes[i],c1,c2,c3,eG_t,conf);
        if (tmpDist<minDist)
        {
            minDist=tmpDist;
            especieCercana=i;
        }
    }
    //Si la distancia entre la especie más cercana y el indexpob es menor que el threshold, retorna el número de la especie (i)
    if (minDist<(conf->threshold+(((randL(conf))-0.5)*conf->threshold*conf->maxDesvThEspecies)))
        return(especieCercana);
    else
    {
        //Si es mayor que el threshold, crea una nueva especie
        ////TODO: Eliminar especies que lleven varias generaciones sin mejorar fitness (Hacer nuevo arreglo de numGeneracsSinCambioEnFitnessconf->representantes).
        //Incrementa el número de especies
        conf->numEspecies++;
        // Si genera nueva especie, incrementa el threshold para mantenerlo en el máximo nivel posible.
        conf->threshold*=(1+(2*conf->porcentVarTh));
        //reserva memoria para el nuevo arreglo de especies con tamaño conf->numEspeciessizeof(unsigned  int)
        if((conf->representantes=(unsigned *)realloc(conf->representantes,sizeof(unsigned  int)*conf->numEspecies))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\Error 654 en asignarEspecie() llamando a realloc");
            return(UINT_MAX);
        }
        //reserva memoria para el nuevo arreglo de conf->contGeneracSinMejora
        if((conf->contGeneracSinMejora=(unsigned *)realloc(conf->contGeneracSinMejora,sizeof(unsigned  int)*conf->numEspecies))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\Error 654 en asignarEspecie() llamando a realloc");
            return(UINT_MAX);
        }
        //coloca al indexpob como representante de la nueva especie
        conf->representantes[conf->numEspecies-1]=indexpob;
        //coloca en 0 el número de generaciones sin mejora para la nueva especie
        conf->contGeneracSinMejora[conf->numEspecies-1]=0;
        //POSIBLE PROBLEMA: se generan demasiadas especies.
        return(conf->numEspecies-1);
    }
    return(1);
}
/*
unsigned evaluarEspecie(unsigned inicializar,unsigned especie,unsigned maxBufferSize,char *fileNameGTDv1, TConfig* conf)  // OPTIMIZADA, NMR
{
// evalúa toda la población y deja el valor post fSigma en cada nodo.
// y calcula el fitness basado en el que se va acumulando con cada evaluación de cada genoma.
// Los archivos de entrada y salida deben estar previamente abiertos para lectura binaria br
// retorna 0 si hay error, 1 si ok.
// //TODO: URGENTE PARA FX parámetro repeticiones para pasar los archivos de entrada repetidas veces por las redes neuronales al realizar las evaluaciones.
    unsigned leidosIn=0; // guarda el número de datos leidos del archivo de entrada.
    unsigned mBuffer;
    float* inputs; // arreglo de datos leidos del archivo de entrada
    unsigned menorLeido;
    unsigned tmpSize;
    unsigned i;
    unsigned j;
    unsigned k;
    float* tmpin;
    float* tmpout;
    double pasadas=0;// usado para normalizar el fitness
    hdrGTDv1 header; // usado apra leer el ancabezado de GTDv1
    // Abre archivo GTDv1
    if((conf->fIn=fopen(fileNameGTDv1,"rb"))==NULL)
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 57 en funcion evaluarPob() llamando a fopen(%s,\"br\")",fileNameGTDv1);
        return(0);
    }
    // Lee encabezado GTDv1
	if (fread(&header, sizeof(hdrGTDv1), 1, conf->fIn)!=1)
	{
		fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>Error 58 en función evaluarPob() llamando a fread(%s)\n",fileNameGTDv1);
		return(0);
	}
    // calcula cuantos grupos de float debe leer multiplos de nEntradas menores que maxBuffersize
    mBuffer = (unsigned  int)floor((float)maxBufferSize/(float)(header.numEntradas+header.numSalidas));
    // compara los multiplicadores y selecciona el menor como mBuffer
	// selecciona por defecto el de entradas debido a que generalmente se usan más entradas que salidas
    // reserva memoria para los arreglos de datos leidos  inputs y outpues
    tmpSize=sizeof(float)*mBuffer*(header.numEntradas+header.numSalidas);
    if (tmpSize<32000*maxBufferSize*4)// verifica que no haya más de 32k entradas y salidas especificadas en el archivo de entrada
    {
		if((inputs=malloc(tmpSize))==NULL)
		{
			fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 57 en funcion evaluarEspecie() llamando a malloc(%u)",(unsigned  int)(sizeof(float)*mBuffer*header.numEntradas));
			return(0);
		}
    }
	else
	{
		fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 58 en funcion evaluarEspecie() numEntradas + numSalidas >32k leyendo desde archivo GTD");
		return(0);
	}
	//7am 7328064 avianca pasto maletas# AV986484 y 82 y 83 reporte a nombre de Susana Bastidas
    // coloca los valores de cada neurona en 0 para comenzar la evaluación. si incializar=1;
    if (inicializar == 1)  // no debe ser siempre 1? igual que en evaluarpob?
    {
        for (i=0; i<conf->sizePob; i++)
        {
            if (conf->pob[i].especie == especie )
            {
                // Actualiza los punteros a conecHijo de los nodos del genoma.
                if (actualizarPNodos(i,conf)==0)
                {
                    fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 1150 en fnción evaluarPob() llamando a actualizarPNodos()\n");
                    free(inputs);
                    return(0);
                }
                for (j=0; j<conf->pob[i].totalNodos; j++)
                {
                    if (conf->pob[i].nodo[j].nodeFunction!=3)
                    {
                        conf->pob[i].nodo[j].valor=0;
                    }
                }
                conf->pob[i].fitness=0; // Inicializa el acumulador de fitness
            }
        }
    }
    // hace repeticiones deben ser 0 o más
    for(k=0; k<=conf->repTrain; k++)
    {
        // evalúa toos los genomas en todas las entradas y salidas
        while(!feof(conf->fIn)) // lee hasta que termine cualquiera de los dos archivos
        {
            if (header.tamRegistros==4)
            {
            	tmpSize=4*(header.numEntradas+header.numSalidas);
            	if (tmpSize<4*32000) leidosIn = (unsigned int)fread(inputs,tmpSize,(size_t)mBuffer,conf->fIn);
				menorLeido=leidosIn;
				for (i=0; i<menorLeido; i++)
				{
				    tmpSize= i*(header.numEntradas+header.numSalidas);
					if (tmpSize<conf->sizePob*32000) tmpin=inputs+(tmpSize);
                    tmpSize=(i*(header.numEntradas+header.numSalidas))+header.numEntradas;
					if (tmpSize<conf->sizePob*32000) tmpout=inputs+(tmpSize);
					pasadas=pasadas+1;
					// evalúa entradas y salidas para todos los genomas HERE!!!!ERROR
					for (j=0; j<conf->sizePob; j++)
					{
						if (conf->pob[j].especie == especie)
						{

							// fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>Ins = %3.3f,%3.3f Out = %3.3f\n",tmpin[0],tmpin[1],tmpout[0]);
							if(evaluarGenoma(j,&(conf->pob[j]),0, tmpin,tmpout,conf)==0) ////TODO: HACER PARAMETRO GLOBAL conf->PRIMERO PARA DIFERENTES APLICACIONES
							{
								fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 59 en funcion evaluarPob() llamando a evaluarGenoma(genoma=%u,iteracion=%u )\n",j,i);
								free(inputs);
								return(0);
							}
						}
                    }
                }
            }
            else //FALTA: implementar manejo para cuendo se usa dato double.
            {
            	fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 60 en función evaluarEspecie() manejo de pesos double aún no se ha implementado.\n");
            	free(inputs);
            	return(0);
            }
        }
        // envía al inicio el prt del archivo de entradas y salidas
        rewind(conf->fIn);
    }
    // Calcula los fitness a partir del acumulado de cada genoma
    if (pasadas>0)
	{
        for (i=0; i<conf->sizePob; i++)
        {
            if (conf->pob[i].especie == especie )
            {
                conf->pob[i].fitness/=pasadas;
                conf->pob[i].fitness=1.0-conf->pob[i].fitness;
            }
        }
	}
    // cierra archivos de entrada y salida
    fclose(conf->fIn);
    // libera los punteros usados por los punteros a las conexhijos usados durante eval (deben ser rearmados nuevamente con actualizarPNodos).
    for (i=0; i<conf->sizePob; i++)
    {
        if (conf->pob[i].especie == especie )
        {
            for(j=0; j<conf->pob[i].totalNodos; j++)
            {
                if ((conf->pob[i].nodo[j].conexHijo!=NULL)&&(conf->pob[i].nodo[j].contHijos>0)&&(conf->pob[i].nodo[j].nodeFunction!=3))
                {
                    free(conf->pob[i].nodo[j].conexHijo);
                }
            }
        }
    }
    // libera memoria de las listas inputs y outputs
    if (inputs!=NULL) free(inputs);
    // TODO: hacer Funcion bufferize para leer las entradas y salidas en arreglos globales para que se lean una sola vez del disco durante todo el ciclo ppal
    return(1);
}
*/
unsigned nuevaInnovNodo(unsigned nodoIn, unsigned nodoOut, TConfig* conf)  //OPTIMIZADA
{
//Retorna el número innovación(número de nodo) buscandolo en la lista de innovaciones de Nodo, si no lo encuentr, lo crea.
//Parámetros: nodoIn, nodoOut = innovnum de nodos de conexión unsigned  interrumpida para agregar el nuevo nodo.
//RETORNA -1 si hay error, numero de innovación si ok.
//OPTIMIZACIÓN: Puede univicarse nuevaInnovNodo y nuevaInnovCon en una sola pasando el puntero a la lista correspondiente como parámetro.
    unsigned i=0;
    TListaInnov* innovNodoIn=NULL;
    TNodoOut* ra_result=NULL;
    ////TODO: VERIFICAR tamaño de conf->listaInnovNodo antes de referenciar nidoIn como index
    if ((nodoIn>=conf->contInnovNodo)||(nodoOut>=conf->contInnovNodo))  //Verifica que los nodos In y Out sean menores que el contInnovIn
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\Error 1 en nuevaInnovNodo(%u,%u) : nIn o nOut >maxN (%u)",nodoIn,nodoOut,conf->contInnovNodo);
        return(UINT_MAX);
    }
    if ((innovNodoIn=&(conf->listaInnovNodo[nodoIn]))->numOut==0) //Si es la primeravez que se agrega una innovación para este nodo, la adiciona.
    {
//TODO: probando quitando el free ya que no hay necesidad porque es null
        //if (!=NULL) free((void *)innovNodoIn->nodoOut);
        if((innovNodoIn->nodoOut=(TNodoOut *)calloc(1,(unsigned  int)sizeof(TNodoOut)))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 2 en funcion nuevaInnovNodo(%u,%u) llamando a calloc(1,%u)\n",nodoIn,nodoOut,(unsigned  int)sizeof(TNodoOut));
            return(UINT_MAX);
        }
        innovNodoIn->nodoOut[0].nodoOut=nodoOut;
        innovNodoIn->nodoOut[0].innovNum=conf->contInnovNodo;
        innovNodoIn->numOut++;
        // reserva memoria para el nuevo nodo como entrada de la lista de innovcon y nodo.(se debe hacer en las 2)
        if((conf->listaInnovNodo=(TListaInnov*) realloc(conf->listaInnovNodo,(sizeof(TListaInnov)*(conf->contInnovNodo+1))))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 2.1 en funcion nuevaInnovNodo(%u,%u) llamando a realloc(%s,%u)\n",nodoIn,nodoOut,innovNodoIn->nodoOut==NULL?"NULL":"iNI->nodoOut",(innovNodoIn->numOut+1)*(unsigned int)sizeof(TNodoOut));
            return(UINT_MAX);
        }
        if((conf->listaInnovCon=(TListaInnov*) realloc(conf->listaInnovCon,(sizeof(TListaInnov)*(conf->contInnovNodo+1))))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 2.2 en funcion nuevaInnovNodo(%u,%u) llamando a realloc(%s,%u)\n",nodoIn,nodoOut,innovNodoIn->nodoOut==NULL?"NULL":"iNI->nodoOut",(innovNodoIn->numOut+1)*(unsigned int)sizeof(TNodoOut));
            return(UINT_MAX);
        }
        //inicializa los nuevos elementos de las listas (el index del arreglo para las dos es el mismo)
        conf->listaInnovNodo[conf->contInnovNodo].numOut=0;
        conf->listaInnovNodo[conf->contInnovNodo].nodoOut=NULL;
        conf->listaInnovCon[conf->contInnovNodo].numOut=0;
        conf->listaInnovCon[conf->contInnovNodo].nodoOut=NULL;
        conf->contInnovNodo++;
        return(conf->contInnovNodo-1);
    }
    else //Si no es la primera vez, busca LA innovación en la lista y la retorna si la encuentra
    {
        for (i=0; i<innovNodoIn->numOut; i++)
        {
            if(innovNodoIn->nodoOut[i].nodoOut==nodoOut)
            {
                return(innovNodoIn->nodoOut[i].innovNum);
            }
        }//Si no la encuentra, en la lista, la crea con realloc para obntener más memoria para el arreglo de salidas.
        if((ra_result=(TNodoOut *) realloc(innovNodoIn->nodoOut,(sizeof(TNodoOut)*(innovNodoIn->numOut+1))))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 3 en funcion nuevaInnovNodo(%u,%u) llamando a realloc(%s,%u)\n",nodoIn,nodoOut,innovNodoIn->nodoOut==NULL?"NULL":"iNI->nodoOut",(innovNodoIn->numOut+1)*(unsigned int)sizeof(TNodoOut));
            return(UINT_MAX);
        }
        innovNodoIn->nodoOut=ra_result;
        innovNodoIn->nodoOut[innovNodoIn->numOut].nodoOut=nodoOut;
        innovNodoIn->nodoOut[innovNodoIn->numOut].innovNum=conf->contInnovNodo;
        innovNodoIn->numOut++;
        // reserva memoria para el nuevo nodo como entrada de la lista de innovcon y nodo.(se debe hacer en las 2)
        if((conf->listaInnovNodo=(TListaInnov*) realloc(conf->listaInnovNodo,(sizeof(TListaInnov)*(conf->contInnovNodo+1))))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 2.1 en funcion nuevaInnovNodo(%u,%u) llamando a realloc(%s,%u)\n",nodoIn,nodoOut,innovNodoIn->nodoOut==NULL?"NULL":"iNI->nodoOut",(innovNodoIn->numOut+1)*(unsigned int)sizeof(TNodoOut));
            return(UINT_MAX);
        }
        if((conf->listaInnovCon=(TListaInnov*) realloc(conf->listaInnovCon,(sizeof(TListaInnov)*(conf->contInnovNodo+1))))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 2.2 en funcion nuevaInnovNodo(%u,%u) llamando a realloc(%s,%u)\n",nodoIn,nodoOut,innovNodoIn->nodoOut==NULL?"NULL":"iNI->nodoOut",(innovNodoIn->numOut+1)*(unsigned int)sizeof(TNodoOut));
            return(UINT_MAX);
        }
        //inicializa los nuevos elementos de las listas (el index del arreglo para las dos es el mismo)
        conf->listaInnovNodo[conf->contInnovNodo].numOut=0;
        conf->listaInnovNodo[conf->contInnovNodo].nodoOut=NULL;
        conf->listaInnovCon[conf->contInnovNodo].numOut=0;
        conf->listaInnovCon[conf->contInnovNodo].nodoOut=NULL;
        conf->contInnovNodo++;
        return(conf->contInnovNodo-1);
    }

}

unsigned nuevaInnovCon(unsigned nodoIn, unsigned nodoOut, TConfig* conf)  // OPTIMIZADA
{
//Retorna el número de innovación para una conexión buscandolo en la lista de innovaciones de conex. si no lo encuentra lo crea.
//retrona -1 si hay error
//Parámetros: nodoIn,nodoOut = nodos (innovNums ) de origen y destino de la conex.
//OPTIMIZACIÓN: Puede univicarse nuevaInnovNodo y nuevaInnovCon en una sola pasando el puntero a la lista correspondiente como parámetro.
    unsigned i=0;
    TListaInnov* innovNodoIn=NULL;
    TNodoOut* ra_result=NULL;
    if ((nodoIn>=conf->contInnovNodo)||(nodoOut>=conf->contInnovNodo))  //Verifica que los nodos In y Out sean menores que el contInnovIn
    {
        fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\Error 4 en nuevaInnovCon(%u,%u) : nIn o nOut >maxN (%u)",nodoIn,nodoOut,conf->contInnovNodo);
        return(UINT_MAX);
    }
////TODO: cada vex que se haga malloc o realloc hacer inicialización de los elementos especialmente listas y sus elementos y sublistas.
//y probablemente haciendo realloc de la conf->listaInnovNodo actual a una nueva de tamaño i si i >conf->maxInnovNumConex
    innovNodoIn=&(conf->listaInnovCon[nodoIn]);
    if (innovNodoIn->numOut>0)
    {
        if ((innovNodoIn->nodoOut[innovNodoIn->numOut-1].innovNum>=conf->contInnovCon)||(innovNodoIn->nodoOut[innovNodoIn->numOut-1].innovNum<0))
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 59.45 en nuevaInnovCon(%u,%u) innovNum de lastConexIN=%u, maxIC=%u",nodoIn,nodoOut,innovNodoIn->nodoOut[innovNodoIn->numOut-1].innovNum,conf->contInnovCon);
            return(UINT_MAX);
        }
    }
    if (innovNodoIn->numOut==0) //Si es la primeravez que se agrega una innovación para este nodo, la adiciona.
    {
//TODO: probando quitando el free ya que no hay necesidad porque es null
        //if (!=NULL) free((void *)innovNodoIn->nodoOut);
        if((innovNodoIn->nodoOut=(TNodoOut *)calloc(1,(unsigned  int)sizeof(TNodoOut)))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 5 en funcion nuevaInnovCon(%u,%u) llamando a calloc(1,%u)",nodoIn,nodoOut,(unsigned  int)sizeof(TNodoOut));
            return(UINT_MAX);
        }
        innovNodoIn->nodoOut[0].nodoOut=nodoOut;
        innovNodoIn->nodoOut[0].innovNum=conf->contInnovCon;
        innovNodoIn->numOut=1;
        conf->contInnovCon++;
        return(conf->contInnovCon-1);
    }
    else //Si no es la primera vez, busca LA innovación en la lista y la retorna si la encuentra
    {

        for (i=0; i<innovNodoIn->numOut; i++)
        {
            //fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nHello9.5 i=%u, numout=%u,nodoIn=%u,nodoOut=%u",i,innovNodoIn->numOut,nodoIn,nodoOut); //Problema en arreglo nodoOut para index 0 parece que es NULL
            if(innovNodoIn->nodoOut[i].nodoOut==nodoOut)
            {
                return(innovNodoIn->nodoOut[i].innovNum);
            }
        }//Si no la encuentra, en la lista, la crea con realloc para obntener más memoria para el arreglode salidas.
        if((ra_result=(TNodoOut*)malloc(sizeof(TNodoOut)*(innovNodoIn->numOut+1)))==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 6 en funcion nuevaInnovCon(%u,%u) llamando a realloc(%s,%u)\n",nodoIn,nodoOut,ra_result==NULL? "NULL":"iNI->nodoOut",(innovNodoIn->numOut+1)*(unsigned int)sizeof(TNodoOut));
            return(UINT_MAX);
        }
        ra_result=(TNodoOut*)memcpy(ra_result,innovNodoIn->nodoOut,sizeof(TNodoOut)*innovNodoIn->numOut);
        if (ra_result!=NULL) free(innovNodoIn->nodoOut);
        innovNodoIn->nodoOut=ra_result;
        innovNodoIn->nodoOut[innovNodoIn->numOut].nodoOut=nodoOut;
        innovNodoIn->nodoOut[innovNodoIn->numOut].innovNum=conf->contInnovCon;
        innovNodoIn->numOut++;
        conf->contInnovCon++;
        return(conf->contInnovCon-1);
    }
}
//Mi nombre es exterminio.



