/** Funciones de manejo de genes - C file
	usan/modifican genes únicamente.
*/

#ifndef PARAMS_H_INCLUDED
#include "params.h"
#define PARAMS_H_INCLUDED
#endif
#ifndef AUXILIARES_H_INCLUDED
#include "auxiliares.h"
#define AUXILIARES_H_INCLUDED
#endif
#include "gen.h"


void calcularLimThOneNode(unsigned indexPob, unsigned nodeIndex, TConfig* conf)  //OPTIMIZADA
{
//Calcula maxth y minTh para un solo nodo de un genoma
    unsigned k;
    //coloca en 0 maxth y minTh
    conf->pob[indexPob].nodo[nodeIndex].maxTh=0;
    conf->pob[indexPob].nodo[nodeIndex].minTh=0;
    //acumula
    for (k=0; k<conf->pob[indexPob].totalConexiones; k++)
    {
        if (conf->pob[indexPob].conex[k].indexOut==nodeIndex) // //TODO: NODOOUT ES EL IINOVNUM O EL INDEX DEL NODO?????
            if (conf->pob[indexPob].conex[k].enabled==1)
            {
                if (conf->pob[indexPob].conex[k].peso>0)
                {
                    conf->pob[indexPob].nodo[conf->pob[indexPob].conex[k].indexOut].maxTh+=conf->pob[indexPob].conex[k].peso;
                }
                else
                {
                    conf->pob[indexPob].nodo[conf->pob[indexPob].conex[k].indexOut].minTh+=conf->pob[indexPob].conex[k].peso;
                }
            }
    }
    //en caso de tsigma=-1,1 se ajusta cada max y min th
    if (conf->tSigma>1000)
    {
        conf->pob[indexPob].nodo[nodeIndex].maxTh+=fabs(conf->pob[indexPob].nodo[nodeIndex].minTh);
        conf->pob[indexPob].nodo[nodeIndex].minTh= -conf->pob[indexPob].nodo[nodeIndex].maxTh;
    }
    //controla que th se encuentre entre los límites.
    /*	if (conf->pob[indexPob].nodo[nodeIndex].thNodo<=conf->pob[indexPob].nodo[nodeIndex].maxTh){
    		if (conf->pob[indexPob].nodo[nodeIndex].thNodo<conf->pob[indexPob].nodo[nodeIndex].minTh){
    			conf->pob[indexPob].nodo[nodeIndex].thNodo=conf->pob[indexPob].nodo[nodeIndex].minTh;//corrige si min
    		}
    	}
    	else{
    			conf->pob[indexPob].nodo[nodeIndex].thNodo=conf->pob[indexPob].nodo[nodeIndex].maxTh; //corrige si max
    	}
    */
}
////TODO:  MODIFICAR EVALUARGENOMA PARA USAR INDEXOUT y TODO LO QUE USE NODOOUT, también debe usar indexout.
//Luego terminar de optimizar funciones  de GEN
void calcularLimThOneGenome(unsigned indexPob, TConfig* conf)  //OPTIMIZADA
{
//Calcula maxTh y minTh para un solo genoma indexPob.
    unsigned j;
    //coloca en 0 maxth y minTh
    for (j=0; j<conf->pob[indexPob].totalNodos; j++)
    {
        conf->pob[indexPob].nodo[j].maxTh=0;
        conf->pob[indexPob].nodo[j].minTh=0;
    }
    //acumula
    for (j=0; j<conf->pob[indexPob].totalConexiones; j++)
    {
        if (conf->pob[indexPob].conex[j].enabled==1)
        {
            if (conf->pob[indexPob].conex[j].peso>0)
            {
                conf->pob[indexPob].nodo[conf->pob[indexPob].conex[j].indexOut].maxTh+=conf->pob[indexPob].conex[j].peso;
            }
            else
            {
                conf->pob[indexPob].nodo[conf->pob[indexPob].conex[j].indexOut].minTh+=conf->pob[indexPob].conex[j].peso;
            }
        }
    }
    //en caso de tsigma=-1,1 se ajusta cada max y min th
    if (conf->tSigma>1000)
    {
        for (j=0; j<conf->pob[indexPob].totalNodos; j++)
        {
            conf->pob[indexPob].nodo[j].maxTh+=fabs(conf->pob[indexPob].nodo[j].minTh);
            conf->pob[indexPob].nodo[j].minTh= -conf->pob[indexPob].nodo[j].maxTh;
        }
    }
    //controla que th se encuentre entre los límites.
    /*	for (j=0;j<conf->pob[indexPob].totalNodos;j++){
    		if (conf->pob[indexPob].nodo[j].thNodo<=conf->pob[indexPob].nodo[j].maxTh){
    			if (conf->pob[indexPob].nodo[j].thNodo<conf->pob[indexPob].nodo[j].minTh){
    				conf->pob[indexPob].nodo[j].thNodo=conf->pob[indexPob].nodo[j].minTh;//corrige si min
    			}
    		}
    		else{
    			conf->pob[indexPob].nodo[j].thNodo=conf->pob[indexPob].nodo[j].maxTh; //corrige si max
    		}
    	}*/
}

void calcularLimTh(TConfig* conf)  //OPTIMIZADA
{
// Calcula los límites maxTh y minTh para cada nodo de cada genoma de la pob, también altera Th si se sale de los límites..
// se debe llamar DESPUES DE evaluarpob, es decir, una sola vez al final del ciclo ppal. y una vez después de primeragen.
    //Para cada indexpob, pera cada conex, acumule max y minTh para el nodo de destino de la conex.
    //si tsigma es -1,1 se debe colocar el max=abs(sum de negativos) + sumpositivos y el min=-max
    //tsigma=0 = sigma(0,1),		y = 1 / (1 + exp (- D * x))
    //		1 = sigma aprox (0,1)	y = 0.5 + x * (1 - abs(x) / 2), y=0 si x<=-1, y=1 si x>=1
    //		2 = elliot (0,1 )	y = (x / 2) / (1 + |x|) + 0.5
    //		3 = gauss(-1,1),	y = exp(- x * x)
    //		4 = elliot,(-1,1)	y = x / (1 + |x|)
    //		5 = tanh(-1,1),		y = 2 / (1 + exp(-2 * x)) - 1

    unsigned i;
    unsigned j;
    //para cada indexpob

    //acumula en max y minTh
    for (i=0; i<conf->sizePob; i++ )
    {
        //coloca en 0 maxth y minTh
        for (j=0; j<conf->pob[i].totalNodos; j++)
        {
            conf->pob[i].nodo[j].maxTh=0;
            conf->pob[i].nodo[j].minTh=0;
        }
        //acumula
        for (j=0; j<conf->pob[i].totalConexiones; j++)
        {
            if (conf->pob[i].conex[j].enabled==1)
            {
                if (conf->pob[i].conex[j].peso>0)
                {
                    conf->pob[i].nodo[conf->pob[i].conex[j].indexOut].maxTh+=conf->pob[i].conex[j].peso;
                }
                else
                {
                    conf->pob[i].nodo[conf->pob[i].conex[j].indexOut].minTh+=conf->pob[i].conex[j].peso;
                }
            }
        }
        //en caso de tsigma=-1,1 se ajusta cada max y min th
        if (conf->tSigma>1000)
        {
            for (j=0; j<conf->pob[i].totalNodos; j++)
            {
                conf->pob[i].nodo[j].maxTh+=fabs(conf->pob[i].nodo[j].minTh);
                conf->pob[i].nodo[j].minTh= -conf->pob[i].nodo[j].maxTh;
            }
        }
        //controla que th se encuentre entre los límites.
        /*	for (j=0;j<conf->pob[i].totalNodos;j++){
        		if (conf->pob[i].nodo[j].thNodo<=conf->pob[i].nodo[j].maxTh){
        			if (conf->pob[i].nodo[j].thNodo<conf->pob[i].nodo[j].minTh){
        				conf->pob[i].nodo[j].thNodo=conf->pob[i].nodo[j].minTh;//corrige si min
        			}
        		}
        		else{
        			conf->pob[i].nodo[j].thNodo=conf->pob[i].nodo[j].maxTh; //corrige si max
        		}
        	}*/
    }
}

void perturbarPeso(unsigned index, float porcent, float porcentConex, unsigned usarQuadPorcent, unsigned quadBase, unsigned tipoProb, float minProb, TConfig* conf)  //OPTIMIZADA: //TODO deberia ir a genoma.c y  otros tipos de pert, quadPorcent y hacer Funcion perturbarTh si es necesaria
{
// Perturba uniformemente en un valor aleatorio máximo de un porcentaje del rango de pesos, un porcentaje del total de conexiones conexiones con
// probabilidad uniforme, lineal o cuadrática dependiente la antiguedad de la conexión. También se cuenta con la opción de decrementar el porcentaje
// de conexiones a perturbqar si se supera cierto número de conexiones en total para evitar grandes cambios de funcionamiento en genomas grandes.
// Parámetros: 	index = indexpob del genoma al que se le aplicará la perturbación.
//				porcent = (0,1) porcentaje del rango inicial de pesos máximo con el que se perturbarán (dependiendo de tipoProb) las conexiones seleccionadas
//				porcentConex = (0,1) porcentaje máximo de conexiones que se perturbarán uniformemente, el resto se coloca aleatorio en rango inicial  (def=0.9)
//				usarQuadPorcent = 	si es 1 usa el param quadBase y porcentConex para calcular el porcentaje de conexiones a perturbar, puede
//									servir paqra genomas grandes, donde se deben modificar menos conexiones a medida que el genoma incrementa el número de conexiones.
//				quadBase = número de conexiones del genoma a partir del cual se usa la regla de reducción cuadrática del porcentaje de conex perturbadas
//				tipoProb = 	0 = uniforme : (default) utiliza un mismo valor de perturbación para cada genoma
//							1 = lineal : selecciona aleatoriamente en una distribución lineal que incrementa desde minProb en las conex más antiguas hasta 100% en las conex más nuevas.
//							2 = cuadrática : selecciona aleatoriamente en una distrib cuadratica que se incrementa desde minProb para la conex más antigua hasta 100% en la más nueva.
//							3 o más = exponencial : selecciona aleatoriamente en una distribuación exponencial que se incrementa desde minProb hasta 100% para la conex más nueva.
//				minProb = mínimo de probabilidad para la conexión más antigua que se e3escoge para las distribuciones lineal, cuadrática o exponencial.
    unsigned i;
    unsigned totConex2;
    float pertU=0;
    float pertN=0;

    //Perturba pesos uniformemente
    if (tipoProb==0)
    {
        //para todas las conexiones:
        for (i=0; i<conf->pob[index].totalConexiones; i++)
        {
            if ((randL(conf))<porcentConex)
            {
                pertU=(((2*randL(conf))-1)*(conf->pesoMaxInicial-conf->pesoMinInicial))*porcent; //magnitud de la pert uniforme
                conf->pob[index].conex[i].peso+=pertU;
            }
            else  //para el porcentaje restante, calcula perturbación aleatoria para la conex con el 100% del rango inicial
            {
                pertN=(((2*randL(conf))-1)*(conf->pesoMaxInicial-conf->pesoMinInicial)); //magnitud de la pert
                conf->pob[index].conex[i].peso+=pertN;
            }
        }
    }
    if (tipoProb==1)  //perturbación de pesos con magnitud de variación lineal
    {
        //para todas las conexiones:
        for (i=0; i<conf->pob[index].totalConexiones; i++)
        {
            if (((float)randL(conf))<porcentConex)
            {
                pertU=(((2*randL(conf))-1)*(conf->pesoMaxInicial-conf->pesoMinInicial))*porcent*((float)(i+1)/conf->pob[index].totalConexiones); //magnitud de la pert uniforme
                conf->pob[index].conex[i].peso+=pertU;
            }
            else  //para el porcentaje restante, calcula perturbación aleatoria para la conex con el 100% del rango inicial
            {
                pertN=(((2*randL(conf))-1)*(conf->pesoMaxInicial-conf->pesoMinInicial)); //magnitud de la pert
                conf->pob[index].conex[i].peso+=pertN;
            }
        }
    }
    if (tipoProb==2)  //perturbación de pesos con magnitud de variación cuadrática
    {
        totConex2=(conf->pob[index].totalConexiones*conf->pob[index].totalConexiones);
        //para todas las conexiones:
        for (i=0; i<conf->pob[index].totalConexiones; i++)
        {
            if (((float)randL(conf))<porcentConex)
            {
                pertU=(((2*randL(conf))-1)*(conf->pesoMaxInicial-conf->pesoMinInicial))*porcent*((float)(i+1)*(i+1)/totConex2); //magnitud de la pert uniforme
                conf->pob[index].conex[i].peso+=pertU;
            }
            else  //para el porcentaje restante, calcula perturbación aleatoria para la conex con el 100% del rango inicial
            {
                pertN=(((2*randL(conf))-1)*(conf->pesoMaxInicial-conf->pesoMinInicial)); //magnitud de la pert
                conf->pob[index].conex[i].peso+=pertN;
            }
        }
    }
}

unsigned randomizarPesos(unsigned indexpob, float minPeso, float maxPeso, TConfig* conf)  //OPTIMIZADA
{
//Randomiza los pesos de cada conexión del genoma y los thresholds en cada nodo aleatoriamente en un rango determinado por los pesos de las conexiones entrantes
//Retorna 0 si hay error (NULL pounsigned  inter), 1 si ok.
//Parámetros: 	indexpob =indexpob a randomizar
//				float minPeso = límite inferior del rango de números aleatorios
//				float maxPeso = límite superior del rango de números aleatorios
    unsigned i;
    float rango=maxPeso-minPeso;
    //perturba los pesos de todas las conexiones del genoma
    for (i=0; i<conf->pob[indexpob].totalConexiones; i++) //para cada conexión
    {
        if (conf->pob[indexpob].conex==NULL)
        {
            fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 27 en Funcion randomizarPesos(%u,%1.1f,%1.1f) conf->pob[%u].conex=NULL\n",indexpob,minPeso,maxPeso,indexpob);
            return(0);
        }
        else  // escoge ealeatoriamente  el peso de la conex entre minPeso y maxPeso
        {
            conf->pob[indexpob].conex[i].peso=(((float)randL(conf))*rango)+minPeso;
        }
    }
    //calcula los max y min th de este genoma a partir de los pesos
    calcularLimThOneGenome(indexpob,conf);
    /*	for (i=0;i<conf->pob[indexpob].totalNodos;i++){//para cada nodo
    		if (conf->pob[indexpob].nodo==NULL){
    			fclose(conf->logFile);conf->logFile=fopen(conf->fileNameLog,"a+");fprintf(conf->logFile,"<br>\nError 27h en Funcion randomizarPesos(%u,%1.1f,%1.1f) conf->pob[%u].conex=NULL\n",indexpob,minPeso,maxPeso,indexpob);
    			return(0);
    		}
    		else{  // escoge aleatoriamente el th entre max y minTh
    			conf->pob[indexpob].nodo[i].thNodo=(((float)randL(conf))*(conf->pob[indexpob].nodo[i].maxTh-conf->pob[indexpob].nodo[i].minTh))+conf->pob[indexpob].nodo[i].minTh;
    		}
    	} */
    //TODO: verificar si es necesaria la ubicación aleatoria de thresholds
    return(1);
}
