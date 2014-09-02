// H5PHET - FX_GTDGen
// Genera archivo GTD de entrenamiento con datos de forex
// Se deben configurar manualmente los parámetros del GTD de salida.


#property indicator_separate_window
#property indicator_buffers 1  // actual y variación por cada uno de los tres TF: 1M, 15M y 4H
                               // si se va a hacer la salida, debe colocarse en 1 y el buffer
                               // de la salida debe tener index 0 en inicialización
#property indicator_color1 Magenta
//---- indicator parameters
//Agregado por Harvey
int numCandles15M_Max=667; // ANTES ERA 10k.
extern int SkipInicio=1500; // 6336 offset para el inicio.
extern int SkipFinal=1500; // offset para el final.
extern double A_Sigma=2.435; // factor usado en sigma out=exp(-A(In*In))   de las fdt de las neuronas.
extern int pipsTPSL=400; // nuero de pips sobre los que se calcula el masterTrainer2 (default 400)
extern int maxDias=10; // periodos en el futuro para calcular el masterTrainer (4 dias)

//---- indicator buffers
double SNN_Indicator[];//Para el indicador de red neuronal
double MasterIndicator[];  //indicador compuesto 0 = máxima ganancia en venta, 1= máx en compra
double CompraE_S1[];  //es 1 si el valor es el mínimo posible para cada zzorder y 0 si es el max.
double VentaE_S1[];  //para venta
double ZigzagBuffer[];
double HighMapBuffer[];
double LowMapBuffer[];



// Para cargar ANN desde SNN máximo 10000 conexiones
int numIns,numBias,numOuts,numHiddens,numConex;
int conexIn[10000];
int conexOut[10000];
double peso[10000];
double tmpPeso;
double tmpVolume;
int inicStat;
int error;

double valor[9000]; //usado para guardar los valores de salida de las neuronas
double normDMax[1000]; // usado para la normalización dinámica de las entradas max 1k.
double normDMin[1000]; // usado para la normalización dinámica de las entradas max 1k.
double decNormD=0.001; // decremento pir shift de la norma dinámica.
int valorCalculado[9000];  //usado para "marcar" cuando el valor de cada neurona ya se procesó(1) o está pendiente (0).
int totalNodos=0; //numero total de nodos.

// - Variables globales
int level=3; // recounting's depth
bool downloadhistory=false;
int PrimeraVez; // para garantizar que solo se corra y guarde una vez.


int handle; // para manejo del archivos de salida.
//int counted_bars = IndicatorCounted();
int counted_bars = 0;  //Hace 0 el counted bars para que siempre se calculen todos los valores.
int limit,counterZ,whatlookfor;
int shift,back,lasthighpos,lastlowpos;
double val,res;
double curlow,curhigh,lasthigh,lastlow;
double acumulado=0;
double A=0;
int ultimo=0;
int ultimoAnt=0;
int zPoints=0; //número total de picos +valles en el intérvalo
int CandleStart=0; //index del primer candlestick a dibujar
int cuentaEscritos=0; //para obtener el número de floats escritos en el archivo.
double maximo1M, minimo1M, maximo15M, minimo15M, maximo4H, minimo4H; // usados para los límites de normalización
double maximo1Mo=1, minimo1Mo=1, maximo15Mo=1, minimo15Mo=1; //usados para buscar los límites de normalización de hobv.
double maximoTam1M, maximoTam15M, maximoTam4H; //usados apra normalizar tamaño y color de candlesticks
string fileName="";
double tmpHigh,tmpLow;
int idle15M=1, idle4H=1; //usados para regular el periodo de calculo de valores de normalización para hobv
int caso1,caso2,caso3_1,caso3_2;
int maxMinutos;

// Encabezado de GTDv1
int version; //(GTD[1])
int usarSigned;
int tamRegistros;
int numEntradas;
int numSalidas;

//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
int init()
{
    IndicatorBuffers(7);
//---- drawing settings
    SetIndexStyle(0,DRAW_SECTION); 
//---- indicator buffers mapping
    SetIndexBuffer(0,MasterIndicator); 
    SetIndexBuffer(1,SNN_Indicator);
    SetIndexBuffer(2,ZigzagBuffer);
    SetIndexBuffer(3,HighMapBuffer);
    SetIndexBuffer(4,LowMapBuffer);
    SetIndexBuffer(5,CompraE_S1);
    SetIndexBuffer(6,VentaE_S1);

    SetIndexEmptyValue(0,0.0);
    SetIndexEmptyValue(1,0.0);
    SetIndexEmptyValue(2,0.0);    
    SetIndexEmptyValue(5,0.0);
    SetIndexEmptyValue(6,0.0);    
//---- indicator short name
    IndicatorShortName("H5PHET - Master Trainer 2 ("+pipsTPSL+","+maxDias+","+numCandles15M_Max+")");
//---- initialization done
   PrimeraVez=0; // Para que se ejecute solo una vez después de la inicialización.
    

// agregado por H
    // new label object
    ObjectDelete("label_object2");
    if(!ObjectCreate("label_object2", OBJ_LABEL, 0, 0, 0))
    {
        Print("error: cant create label_object! code #",GetLastError());
        return(0);
    }
    ObjectSet("label_object2", OBJPROP_XDISTANCE, 250);
    ObjectSet("label_object2", OBJPROP_YDISTANCE, 0);
    ObjectSetText("label_object2", "H5PHET - Master Trainer 2", 10, "Times New Roman", Green);
    caso1=0;
    caso2=0;
    caso3_1=0;
    caso3_2=0;
    maxMinutos=maxDias*1440;
    return(0);
}

      // esta recorre el arreglo conexOut buscando los que tienen como conexOut a i.
        // y retorna FSIGMA de la sumatoria de sus (calcularValor(j)*peso[j]).
        
    
/***********************************   FUNCIONES    ********************************************************/
/******************************************************************************************/    

/***********************************  escribirGTDhdr()  *************************/
void escribirGTDhdr()
{
   // Escribe el encabezado GTD en el archivo de salida (handler).
   // estructura del encabezado de GTDv1 (Generic Training Data).
// formato GDTv1: [encabezado][d1_entrada_0]..[d1_entrada_n][d1_salida_0]..[d1_salida_n].._entrada_0]..[dn_entrada_0]..[dn_entrada_n][dn_salida_0]..[dn_salida_n]
 	int fileIDv1 = 21255239; // Debe ser siempre GTD<#1> para identificar el tipo de archivo version 1.
   int k;
   int numIns2;
   usarSigned = 1;
   tamRegistros = 4; //4=float, 8=double
   numIns2 = inicStat+24;
   numSalidas = 1;
   k=FileWriteInteger(handle, fileIDv1, LONG_VALUE);
   if (k<0)
   {
      Print("Error escribiendo dato 1   Desc:",GetLastError());   
   }
   k=FileWriteInteger(handle, usarSigned, LONG_VALUE);
   if (k<0)
   {
      Print("Error escribiendo dato 2   Desc:",GetLastError());   
   }
   k=FileWriteInteger(handle, tamRegistros, LONG_VALUE);
   if (k<0)
   {
      Print("Error escribiendo dato 3   Desc:",GetLastError());   
   }
   k=FileWriteInteger(handle, numIns2, LONG_VALUE);
   if (k<0)
   {
      Print("Error escribiendo dato 4   Desc:",GetLastError());   
   }
   k=FileWriteInteger(handle, numSalidas, LONG_VALUE);
   if (k<0)
   {
      Print("Error escribiendo dato 5   Desc:",GetLastError());   
   }
    
}

/***************************************  CalcularEntradas             ****************************************************************/

void calcularEntradas(int i)
// i =shift.
{
    int j,tmi,tmj,tmpj,k,i15,i240;
    double tmp,tmp2;
        // se busca el maximo y el mínimo close entre i y i+Tnormalizac, en los ultimos 96 15M
        // entrada 0: CloseN 1M
        if (MathMod(tmi,100)==0) Print("Progreso=",i);
        maximo1M=-9999;
        minimo1M=9999;
        maximoTam1M=-999;
        tmj=(MathFloor(i/15)+96);
        
        
        for (j=MathFloor(i/15); j<tmj; j++)
        {
            tmp=iHigh(NULL,PERIOD_M15,j);
            if(tmp>maximo1M)
            {
                maximo1M=tmp;
            }
            tmp2=iLow(NULL,PERIOD_M15,j);
            if(tmp2<minimo1M)
            {
                minimo1M=tmp2;
            }
            if ((tmp-tmp2)>maximoTam1M)
            {
                maximoTam1M=(tmp-tmp2);
            }
        }
        // se calcula el valor de salida i con formula en MasterIndicator
        if (maximo1M>minimo1M) // para evitar división por 0
        {
            valor[inicStat+0]=(iClose(NULL,PERIOD_M1,i)-minimo1M)/(maximo1M-minimo1M);  //entre -1 y 1
        }
        else
        {
            valor[inicStat+0]=0;
        }
        // se recortan salidas fuera de los límites
        if (valor[inicStat+0]<-1)
        {
            valor[inicStat+0] = -1;
        }
        if (valor[inicStat+0]>1)
        {
            valor[inicStat+0] = 1;
        }

        // entrada 1: CoseN 1M diff
        valor[inicStat+2]=(((iClose(NULL,PERIOD_M1,i)-iClose(NULL,PERIOD_M1,i+1))/(maximo1M-minimo1M)));


        // entrada 2: CloseN 15M
        // se busca el maximo y el mínimo close entre i y i+Tnormalizac, en los ultimos 100 4H (15 Dias)
        maximo15M=-9999;
        minimo15M=9999;
        maximoTam15M=-999;        
        tmj=(MathFloor(i/240)+100); //4horas = 240 minutos
        for (j=MathFloor(i/240); j<tmj; j++)
        {
            tmp=iHigh(NULL,PERIOD_H4,j);
            if(tmp>maximo15M)
            {
                maximo15M=tmp;
            }
            tmp2=iLow(NULL,PERIOD_H4,j);
            if(tmp2<minimo15M)
            {
                minimo15M=tmp2;
            }
            if ((tmp-tmp2)>maximoTam15M)
            {
                maximoTam15M=(tmp-tmp2);
            }
        }
        // se calcula el valor de salida i con formula
        if (maximo15M>minimo15M) // para evitar división por 0
        {
            valor[inicStat+2]=(iClose(NULL,PERIOD_M1,i)-minimo15M)/(maximo15M-minimo15M);  //entre -1 y 1
        }
        else
        {
            valor[inicStat+2]=0;
        }
        // se recortan salidas fuera de los límites
        if (valor[inicStat+2]<-1)
        {
            valor[inicStat+2] = -1;
        }
        if (valor[inicStat+2]>1)
        {
            valor[inicStat+2] = 1;
        }

        // entrada 3: CloseN 15M diff
        valor[inicStat+4]=(((iClose(NULL,PERIOD_M1,i)-iClose(NULL,PERIOD_M1,i+15))/(maximo15M-minimo15M)));

        // entrada 4: CloseN 4H
        // se busca el maximo y el mínimo close entre i y i+Tnormalizac, en los ultimos 41 semanas (225 dias).
        maximo4H=-9999;
        minimo4H=9999;
        maximoTam4H=-999;
        tmj=(MathFloor(i/7920)+41); //1 semana = 7920 minutos 41 semanas = 225d = 15*15d
        for (j=MathFloor(i/7920); j<tmj; j++)
        {
            tmp=iHigh(NULL,PERIOD_W1,j);
            if(tmp>maximo4H)
            {
                maximo4H=tmp;
            }
            tmp2=iLow(NULL,PERIOD_W1,j);
            if(tmp2<minimo4H)
            {
                minimo4H=tmp2;
            }
            if ((tmp-tmp2)>maximoTam4H)
            {
                maximoTam4H=(tmp-tmp2);
            }
        }
        // se calcula el valor de salida i con formula en MasterIndicator
        valor[inicStat+4]=(iClose(NULL,PERIOD_M1,i)-minimo4H)/(maximo4H-minimo4H);  //entre -1 y 1
        // se recortan salidas fuera de los límites
        if (valor[inicStat+4]<-1)
        {
            valor[inicStat+4] = -1;
        }
        if (valor[inicStat+4]>1)
        {
            valor[inicStat+4] = 1;
        }
        // entrada 5: CloseN 4H diff
        valor[inicStat+6]=(((iClose(NULL,PERIOD_M1,i)-iClose(NULL,PERIOD_M1,i+240))/(maximo4H-minimo4H)));
        // entrada 6: Tamaño de candlestick 1M
        valor[inicStat+6]=(iHigh(NULL,PERIOD_M1,i)-iLow(NULL,PERIOD_M1,i))/(maximoTam1M);  //entre -1 y 1
        // se recortan salidas fuera de los límites (porque se normaliza con aproximación por ahora TODO: normalizar exáctamente)
        if (valor[inicStat+6]<-1)
        {
            valor[inicStat+6] = -1;
        }
        if (valor[inicStat+6]>1)
        {
            valor[inicStat+6] = 1;
        }
        //entrada 7: Tamaño del candlestick 1M diff
        valor[inicStat+8]=(((iHigh(NULL,PERIOD_M1,i)-iLow(NULL,PERIOD_M1,i)-iHigh(NULL,PERIOD_M1,i+1)+iLow(NULL,PERIOD_M1,i+1))/(maximoTam1M)));
        //entrada 8: Tamaño del candlestick 15M
        valor[inicStat+8]=(iHigh(NULL,PERIOD_M1,iHighest(NULL,PERIOD_M1,MODE_HIGH,15,i))-iLow(NULL,PERIOD_M1,iLowest(NULL,PERIOD_M1,MODE_LOW,15,i)) )/(maximoTam15M);  //entre -1 y 1
        // se recortan salidas fuera de los límites
        if (valor[inicStat+8]<-1)
        {
            valor[inicStat+8] = -1;
        }
        if (valor[inicStat+8]>1)
        {
            valor[inicStat+8] = 1;
        }
        // entrada 9: Tamaño de candlestick 15M diff
        valor[inicStat+10]=(((iHigh(NULL,PERIOD_M1,iHighest(NULL,PERIOD_M1,MODE_HIGH,15,i))-iLow(NULL,PERIOD_M1,iLowest(NULL,PERIOD_M1,MODE_LOW,15,i))-iHigh(NULL,PERIOD_M1,iHighest(NULL,PERIOD_M1,MODE_HIGH,15,i+15))+iLow(NULL,PERIOD_M1,iLowest(NULL,PERIOD_M1,MODE_LOW,15,i+15)))/(maximoTam15M)));
        // entrada 10: Tamaño de candlestick 4H
        valor[inicStat+10]=(iHigh(NULL,PERIOD_M1,iHighest(NULL,PERIOD_M1,MODE_HIGH,15,i))-iLow(NULL,PERIOD_M1,iLowest(NULL,PERIOD_M1,MODE_LOW,15,i)))/(maximoTam4H);  //entre -1 y 1
        // se recortan salidas fuera de los límites
        if (valor[inicStat+10]<-1)
        {
            valor[inicStat+10]= -1;
        }
        if (valor[inicStat+10]>1)
        {
            valor[inicStat+10]= 1;
        }
        // entrada 11: Tamaño de candlestick 4H diff
        valor[inicStat+12]=(((iHigh(NULL,PERIOD_M1,iHighest(NULL,PERIOD_M1,MODE_HIGH,240,i))-iLow(NULL,PERIOD_M1,iLowest(NULL,PERIOD_M1,MODE_LOW,240,i))-iHigh(NULL,PERIOD_M1,iHighest(NULL,PERIOD_M1,MODE_HIGH,240,i+240))+iLow(NULL,PERIOD_M1,iLowest(NULL,PERIOD_M1,MODE_LOW,240,i+240)))/(maximoTam4H)));
        // entrada 12: Color de candlestick 1M
        valor[inicStat+12]=(((iClose(NULL,PERIOD_M1,i)-iOpen(NULL,PERIOD_M1,i))/(maximoTam1M)));  //entre -1 y 1
        // se recortan salidas fuera de los límites
        if (valor[inicStat+12]<-1)
        {
            valor[inicStat+12] = -1;
        }
        if (valor[inicStat+12]>1)
        {
            valor[inicStat+12] = 1;
        }
        //guarda color de 1M en ZigzagBuffer para ser usardo por hobv
        ZigzagBuffer[i]=valor[inicStat+12];
        //entrada 13: Color de caldlestick 1M diff
        valor[inicStat+14]=(((iClose(NULL,PERIOD_M1,i)-iOpen(NULL,PERIOD_M1,i)-iClose(NULL,PERIOD_M1,i+1)+iOpen(NULL,PERIOD_M1,i+1))/(maximoTam1M)));
        //entrada 14: Color de candlestick 15M
        valor[inicStat+14]=(((iClose(NULL,PERIOD_M1,i)-iOpen(NULL,PERIOD_M1,i+15))/(maximoTam15M)));  //entre -1 y 1
        // se recortan salidas fuera de los límites
        if (valor[inicStat+14]<-1)
        {
            valor[inicStat+14]= -1;
        }
        if (valor[inicStat+14]>1)
        {
            valor[inicStat+14]= 1;
        }
        //guarda color de 15M en compraE_S para ser usado por hobv
        CompraE_S1[i]=valor[inicStat+14];
        // entrada 15: Color de candlestick 15M diff
        valor[inicStat+16]=(((iClose(NULL,PERIOD_M1,i)-iOpen(NULL,PERIOD_M1,i+14)-iClose(NULL,PERIOD_M1,i+15)+iOpen(NULL,PERIOD_M1,i+29))/(maximoTam15M)));
        // entrada 16: Color de candlestick 4H
        valor[inicStat+16]=(((iClose(NULL,PERIOD_M1,i)-iOpen(NULL,PERIOD_M1,i+240))/(maximoTam4H)));  //entre -1 y 1
        // se recortan salidas fuera de los límites
        if (valor[inicStat+16]<-1)
        {
            valor[inicStat+16] = -1;
        }
        if (valor[inicStat+16]>1)
        {
            valor[inicStat+16] = 1;
        }
        //guarda color de 4H en ventaE_S para ser usado por hobv        
        VentaE_S1[i]=valor[inicStat+16];
        // entrada 17: Color de candlestick 4H diff
        valor[inicStat+18]=(((iClose(NULL,PERIOD_M1,i)-iOpen(NULL,PERIOD_M1,i+239)-iClose(NULL,PERIOD_M1,i+240)+iOpen(NULL,PERIOD_M1,i+479))/(maximoTam4H)));
        // entrada 18: hobv (Harvey's On-Balance-Volume) 1M
        // con el color en ZizagBuffer, se calcula el hobv actual
        ZigzagBuffer[i]=ZigzagBuffer[i+1]+(ZigzagBuffer[i]*iVolume(NULL,PERIOD_M1,i));
        // busca entre los últimos 1440 el más alto y el más bajo y los almacena en maximo y minimo
        // debido a lo intensivo de los calculos se calculan los límites de normalización solo cada 15 minutos
        idle15M--;
        if (idle15M<=0)
        {
            //normaliza respecto a 1dia 1440m
            minimo1Mo=2000000000;
            maximo1Mo=-2000000000;
            for (j=0; j<1440; j++)
            {
                if (ZigzagBuffer[i+j]>maximo1Mo) maximo1Mo=ZigzagBuffer[i+j];
                if (ZigzagBuffer[i+j]<minimo1Mo) minimo1Mo=ZigzagBuffer[i+j];
            }
            idle15M=15;
        }
        // Asigna el valor de la entrada 18 (hobv 1M):
        //se hace 1 cada 15 minutos
        valor[inicStat+18]=(ZigzagBuffer[i]-minimo1Mo)/(maximo1Mo-minimo1Mo); //0,1
        // entrada 19: hobv 1M diff
        valor[inicStat+20]=(((ZigzagBuffer[i]-ZigzagBuffer[i+1])/(maximo1Mo-minimo1Mo)));
        // entrada 20: hobv 15M
        // normaliza respecto a 15dias
        tmpVolume=0;
        for (j=0; j<15; j++)
        {
            tmpVolume+=iVolume(NULL,PERIOD_M1,k+j);
        }
        CompraE_S1[i]=CompraE_S1[i+1]+(CompraE_S1[i]*tmpVolume); //calcula el hobv atualsin nomralizar en CompraE_S1
        idle4H--;        // debido a los calculos intensivos necesarios, solo se buscan los periodos de normalizacion cada 4H (240m)
        if (idle4H<=0)
        {
            minimo15Mo=2000000000;
            maximo15Mo=-2000000000;
            for (j=0; j<21600; j++)        // busca entre los últimos 21600 el más alto y el más bajo y los almacena en maximo y minimo
            {
                if (CompraE_S1[i+j]>maximo15Mo) maximo15Mo=CompraE_S1[i+j];
                if (CompraE_S1[i+j]<minimo15Mo) minimo15Mo=CompraE_S1[i+j];
            }
            idle4H=240;
        }
        valor[inicStat+20]=(CompraE_S1[i]-minimo15Mo)/(maximo15Mo-minimo15Mo); // -1,1
        // verifica límites
        if (valor[inicStat+20]>1) valor[inicStat+20]=1;
        if (valor[inicStat+20]<-1) valor[inicStat+20]=-1;
        // entrada 21: hobv 15M diff simplemente es la variacion en los últimos 15 min APROXIMACION
        valor[inicStat+21]=(((valor[inicStat+20]-((CompraE_S1[i+15]-minimo4H)/(maximo15Mo-minimo15Mo)))/2));
        if (valor[inicStat+21]>1) valor[inicStat+20]=1;
        if (valor[inicStat+22]<-1) valor[inicStat+20]=-1;
        // entrada 22: hobv 4H
        // se calcula el color para 15dias (21600m) necesario para calcular el indicador hobv. (lo calcula en MasterIndicator)
         // con el obv sin normalizar en CompraE_S1 (calculado al inicio) se normaliza.
        // normaliza respecto a 15dias m
        tmpVolume=0;
        minimo1M=2000000000;
        maximo1M=-2000000000;
        //busca entre los últimos 1440 el más alto y el más bajo y los almacena en maximo y minimo
        for (j=0;j<41;j++) // 1w = 5.5d = 7920m
        {
            tmpVolume+=iVolume(NULL,PERIOD_W1,MathFloor(i/7920)+j);
        }
        valor[inicStat+22]=((VentaE_S1[i]/tmpVolume));
        // verifica límites
        if (valor[inicStat+22]>1) valor[inicStat+22]=1;
        if (valor[inicStat+22]<-1) valor[inicStat+22]=-1;
        // entrada 23: hobv 4H diff simplemente es la variacion en los últimos 15 min APROXIMACION
        valor[inicStat+23]=(((valor[inicStat+22]-(VentaE_S1[i+15]/tmpVolume))/2));
        if (valor[inicStat+23]>1) valor[inicStat+20]=1;
        if (valor[inicStat+24]<-1) valor[inicStat+20]=-1;        
        // FIN DE calculo de entradas estándar..
        // Para cargar ANN desde SNN máximo 10000 conexiones y 1000 entradas
        i15=i/15;
        i240=i/240;
        for (j=0;j<numEntradas;j++)
        {
            // llama la función calcularIndicador (const Moneda, const TF, shift, double param1, param2, param3,param4, param5)
            //par alos 3 timeframes y 3 diferenciales:
            tmpj=j*6;            
            // tf 1M
            valor[tmpj]=calcularIndicador("EURUSD",PERIOD_M1,j,i);
            // diferencial
            valor[tmpj+1]=valor[tmpj]-calcularIndicador("EURUSD",PERIOD_M1,j, i+1);
            //ft 15M
            valor[tmpj+2]=calcularIndicador("EURUSD",PERIOD_M15,j,i15);
            //diferencial
            valor[tmpj+3]=valor[tmpj+2]-calcularIndicador("EURUSD",PERIOD_M15,j,i15+1);
            //tf 4H
            valor[tmpj+4]=calcularIndicador("EURUSD",PERIOD_H4,j,i240);
            //diferencial
            valor[tmpj+5]=valor[tmpj+4]-calcularIndicador("EURUSD",PERIOD_H4,j,i240+1);
            // para las 6
            for (k=0;k<6;k++)
            {
                // varía (acerca) normas dinámicas
                normDMax[tmpj+k]+=(decNormD*(normDMax[tmpj+k]-normDMin[tmpj+k]));
                normDMin[tmpj+k]-=(decNormD*(normDMax[tmpj+k]-normDMin[tmpj+k]));
                // compara normas con cada valor
                if (valor[tmpj+k]>normDMax[tmpj+k])
                    normDMax[tmpj+k]=valor[tmpj+k];
                if (valor[tmpj+k]<normDMin[tmpj+k])
                    normDMin[tmpj+k]=valor[tmpj+k];
                // verifica que las normas tengan datos válidos
                if (normDMax[tmpj+k]<normDMin[tmpj+k])
                {   
                    // los invierte
                    tmp2=normDMax[tmpj+k];
                    normDMax[tmpj+k]=normDMin[tmpj+k];
                    normDMin[tmpj+k]=tmp2;
                }
                // normaliza cada valor (Máximo=1, mínimo=-1).
                if (normDMax[tmpj+k]!=normDMin[tmpj+k])
                {
                    valor[tmpj+k]=(valor[tmpj+k]-normDMin[tmpj+k])/(normDMax[tmpj+k]-normDMin[tmpj+k]);  //para 0,1
                    valor[tmpj+k]=((2*valor[tmpj+k])-1); //para -1,1

                }
            }
            
/*
con 8
  X
 XX
XXXX
x=enabled.
*/            
        //    Print("EH",tmpj,"=",valor[tmpj]);
        }
        return;
    
}


double calcularIndicador(string symbol, int timeframe, int indicNum, int shift)
{
    // calcula el indicador indicnum con sus params y retorna su valor.
    if (indicNum==0)
        return(iAC(symbol, timeframe, shift));
    else if (indicNum==1)
        return(iAD(symbol, timeframe, shift));
    else if (indicNum==2)
        return(iAlligator(symbol, timeframe, 13, 8, 8, 5, 5, 3, MODE_EMA, PRICE_MEDIAN, MODE_GATORJAW, shift));
    else if (indicNum==3)
        return(iADX(symbol,timeframe,14,PRICE_HIGH,MODE_PLUSDI,shift));
    else if (indicNum==4)
        return(iATR(symbol,timeframe, 12, shift));
    else if (indicNum==5)
        return(iAO(symbol, timeframe, shift));
    else if (indicNum==6)
        return(iBearsPower(symbol, timeframe, 13,PRICE_CLOSE,shift));
    else if (indicNum==7)
        return(iBands(symbol,timeframe,20,2,0,PRICE_LOW,MODE_LOWER,shift));
    else if (indicNum==8)
        return(iBullsPower(symbol,timeframe, 13,PRICE_CLOSE,shift));
    else if (indicNum==9)
        return(iCCI(symbol,timeframe,20,PRICE_TYPICAL,shift));
    else if (indicNum==10)
        return(iDeMarker(symbol,timeframe, 13, shift));
    else if (indicNum==11)
        return(iEnvelopes(symbol,timeframe, 13,MODE_EMA,10,PRICE_CLOSE,0.2,MODE_UPPER,shift));             
    else if (indicNum==12)
        return(iForce(symbol, timeframe, 13,MODE_EMA,PRICE_CLOSE,shift));
    else if (indicNum==13)
        return(iFractals(symbol,timeframe, MODE_UPPER, shift));
    else if (indicNum==14)
        return(iGator(symbol,timeframe, 13, 8, 8, 5, 5, 3, MODE_EMA, PRICE_MEDIAN, MODE_UPPER, shift));
    else if (indicNum==15)
        return(iIchimoku(symbol,timeframe, 9, 26, 52, MODE_TENKANSEN, shift));
    else if (indicNum==16)
        return(iBWMFI(symbol,timeframe,shift));
    else if (indicNum==17)
        return(iMomentum(symbol,timeframe,20,PRICE_CLOSE,shift));
    else if (indicNum==18)
        return(iMFI(symbol,timeframe,14,shift));
    else if (indicNum==19)
        return(iMA(symbol,timeframe,13,8,MODE_EMA,PRICE_MEDIAN,shift));
    else if (indicNum==20)
        return(iOsMA(symbol,timeframe,12,26,9,PRICE_MEDIAN,shift));
    else if (indicNum==21)
        return(iMACD(symbol,timeframe,12,26,9,PRICE_CLOSE,MODE_MAIN,shift));
    else if (indicNum==22)
        return(iOBV(symbol, timeframe, PRICE_CLOSE, shift));
    else if (indicNum==23)
        return(iSAR(symbol,timeframe,0.02,0.2,shift));
    else if (indicNum==24)
        return(iRSI(symbol,timeframe,14,PRICE_CLOSE,shift));
    else if (indicNum==25)
        return(iRVI(symbol, timeframe, 10,MODE_MAIN,shift));
    else if (indicNum==26)
        return(iStdDev(symbol,timeframe,10,0,MODE_EMA,PRICE_CLOSE,shift));
    else if (indicNum==27)
        return(iStochastic(symbol,timeframe,5,3,3,MODE_EMA,0,MODE_SIGNAL,shift));
    else if (indicNum==28)
        return(iWPR(symbol,timeframe,14,shift));
    else
    {
        Print("Error en calcularIndicador() indice no encontrado");
        return(0);
    }
    
}



/**************************************   fSigma        ********************************************************


double fSigma(double fX)
{
    return(MathExp(-A_Sigma*(fX * fX)));
}
*/

/*************************************    calcularValor          ***************************************************************


double calcularValor(int index)
// Calcula el valor de un nodo de la red neuronal simple.   

{
    double acum=0;
    int im;
    // verifica si el valor ya ha sido calculado
    if (valorCalculado[index]==1) return(valor[index]);
    valorCalculado[index]=1;
    // calcula la sumatoria de los valor*peso de las conex que la tengan como destino.
    for (im=0;im<numConex;im++)
    {
        if (conexOut[im]==index) // si alguien la tiene como destino
        {   
            if (valorCalculado[conexIn[im]]==1) // por optimización de velocidad se deja el else.
            {
                acum+=(peso[im]*valor[conexIn[im]]);
            }
            else
            {
                acum+=(peso[im]*calcularValor(conexIn[im]));
            }
        }
    }
    // calcula el valor del nodo
    valor[index]=fSigma(1-acum);
    return(valor[index]);
}
*/

/*************************  Funcion calcularMI ************************-***********************************************************/

double calcularMI(int i)
//  busca en el periodo minutosMax si se cumpe el pipsTPSL para buy o sell y le resta el maximo encontrado hasta ese momento en el sentido contrario
{
    int j;
    double tmpMax;
    double tmpMin;
    double tmpHigh;//variabes usadas para hacer cache de funcs y acelerar calculo
    double tmpLow;
    double tmpClose;
    tmpMax = -200000000;
    tmpMin = 200000000;
    tmpClose=iClose(NULL,PERIOD_M1,i);
    for (j=i; j>(i-maxMinutos); j--)
    {
        tmpHigh=iHigh(NULL,PERIOD_M1,j);
        tmpLow=iLow(NULL,PERIOD_M1,j);
        if (tmpHigh>tmpMax)
        {
            tmpMax=tmpHigh;
            // caso 1: se encontro un TP alto (para sell)
            if ((tmpClose+(Point*pipsTPSL))<tmpMax)
            {
                caso1++;
                return(((Point*pipsTPSL)-(tmpClose-tmpMin))/(2*(Point*pipsTPSL)));
            }
        }
        if (tmpLow<tmpMin)
        {
            tmpMin=tmpLow;
            // caso 2: se encontro un TP bajo (para buy)      
            if ((tmpClose-(Point*pipsTPSL))>tmpMin)
            {
                caso2++;
                return((-(Point*pipsTPSL)-(tmpClose-tmpMax))/(2*(Point*pipsTPSL)));
            }
            
        }
    }
    // caso 3: no se encontró en el periodo de búsqueda ningún TP, se imprime la mayor ganancia en el mejor sentido
    if ((tmpMax-tmpClose)>(tmpClose-tmpMin))
    {
        caso3_1++;
        return(((tmpMax-tmpClose)-(tmpClose-tmpMin))/(2*(Point*pipsTPSL)));
    }
    else
    {
        caso3_2++;
        return(((tmpClose-tmpMin)-(tmpMax-tmpClose))/(2*(Point*pipsTPSL)));    
    }
}

/*************************  Funcion AbrirArchivos ************************-***********************************************************/

void abrirArchivos()
//abre todos los archivos de entrada y salida
{
    int i;
    // para las entradas
    handle=FileOpen("h5phet\\trdata.gtd", FILE_BIN|FILE_WRITE);
    // revisa si hubo error
    if(handle<1)
    {
       Print("No se puede abrir el archivo -",GetLastError());
       error=1;
       return;
    }    
}

/*************************  Funcion cerrarArchivos ************************-***********************************************************/

void cerrarArchivos()
//cierra todos los archivos.
{
        FileClose(handle);
}


/*************************  Funcion escribirTodo ************************-***********************************************************/

void escribirTodo()
// escribe todos los valores en los archivos correspondiente, asume que los archivos están abiertos.
{
    int j,k;
    for (j=0;j<(inicStat+24+numSalidas);j++)
    {
        // TODO: escribir pero en un solo handle. FALTA función escribirencabezado.
        if (tamRegistros==4)
        {
           k=FileWriteDouble(handle, valor[j], FLOAT_VALUE);
        }
        if (tamRegistros==8)
        {
           k=FileWriteDouble(handle, valor[j], DOUBLE_VALUE);
        }
        if (k<0)
        {
            Print("Error escribiendo dato ",j,"   Desc:",GetLastError());   
        }
    }
}


/*************************  Funcion calcularTodo ************************-***********************************************************/

   
// calcula las entradas en valor[]=(0,1) y la señal de entrenamiento en MasterIndicator[]
void calcularTodo()
{

    int err;
    int i,j,k,tmi,tmj;
    datetime tmpFecha;
    /*// Para cargar ANN desde SNN máximo 10000 conexiones
    int numIns,numBias,numOuts,numHiddens,numConex;
    int conexIn[10000];
    int conexOut[10000];
    double peso[10000];
    double tmpPeso;
    */
// variable de valores usada para calcular la salida de cada nodo. las neuronas de entrada y bias no se evaluan, sus valores
// de salida se cargan directamente de las entradas de la red neuronal (4 Entradas de MasterTrainer).
// maximo 9000 neuronas.
    // inicializa los valores en 0
    ArrayInitialize(valor,0.0);
    ArrayInitialize(valorCalculado,0.0);
    ArrayInitialize(normDMax,-1000.0);
    ArrayInitialize(normDMin,1000.0);
    // INICIO de barrido por candlesticks
    tmi=(counted_bars*15); // corregido, antes se sumaban 1500? pero no recuerdo por que.
    CompraE_S1[tmi+1]=0; // NECESARIOs PARA CALCULO DE OBV
    VentaE_S1[tmi+1]=0;
    VentaE_S1[tmi+1]=0;
    if (tmi>0)
    Print("Inicia calculo de entradas");   
    for (i=tmi; i>=0; i--)
    {
        calcularEntradas(i); //Calcula las 24 entradas TODO: optimizar calculo de factores de normalización para velocidad.
        //imprime las entradas y la fecha de cada cstick
        tmpFecha=iTime(NULL,PERIOD_M1,i);
        Print("Fecha: ", TimeYear(tmpFecha),"-",TimeMonth(tmpFecha),"-",TimeDay(tmpFecha)," ",TimeHour(tmpFecha),":",TimeMinute(tmpFecha),"  Candlestick ",i);
/*        Print("E0=",valor[inicStat+1]," E1=",valor[inicStat+1]," E2=",valor[inicStat+3]," E3=",valor[inicStat+3]," E4=",valor[inicStat+5]);
        Print("E5=",valor[inicStat+5]," E6=",valor[inicStat+7]," E7=",valor[inicStat+7]," E8=",valor[inicStat+9]," E9=",valor[inicStat+9]);
        Print("E10=",valor[inicStat+11]," E11=",valor[inicStat+11]," E12=",valor[inicStat+13]," E13=",valor[inicStat+13]," E14=",valor[inicStat+15]);
        Print("E15=",valor[inicStat+15]," E16=",valor[inicStat+17]," E17=",valor[inicStat+17]," E18=",valor[inicStat+19]," E19=",valor[inicStat+19]);
        Print("E20=",valor[inicStat+21]," E21=",valor[inicStat+21]," E22=",valor[inicStat+23]," E23=",valor[inicStat+23]);
*/
        // calcula el entrenador en el buffer MasterIndicator y lo guarda todas las entradas en sus respectivos archivos.
        // usa valores normalizados entre 1 y 0
        if (i>maxMinutos){
            // calcula para cada candlestick en MasterIndicator2 
            MasterIndicator[i] = -calcularMI(i); //EL INDICADOR Quedó al revés
            valor[inicStat+25]=MasterIndicator[i];
            escribirTodo(); //guarda todos los valores en los archivos de salida.
            Print("MI = ",MasterIndicator[i]);
        }
        else
        {
            Print("MI no puede ser calculado pq i=",i," y max=",maxMinutos);
        }
        Print("caso1=",caso1,"  caso2=",caso2, "  caso3_1=",caso3_1,"  caso3_2=",caso3_2);
    } // FIN FOR RECORRIDO DE Candlesticks
    
// Fin de Indicador SNN parte(final) de H5PHET-FX.

//********************************************************************************************************************/
//
// Fin de plataforma de neuroevolución para Forex H5PHET-FX.
// Autor: Harvey Demian Bastidas Caicedo.
// Fecha de Inicio: 02 de noviembre de 2009
// Fecha de Finalización: 08 de febrero de 2011
// Dedicatoria:
//       A mi Dios a mi Madre María del Carmen Caicedo B. a mi padre Duque Manuel H. Bastidas S. a mis Hermanas
//       y a mi novia Juliana Andrade Gómez.
//       "Gracias por darme el mejor mundo para vivir y soñar que pude haber imaginado."
//
//********************************************************************************************************************/
}


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
int start()
{
    numEntradas=29;
    counted_bars = 0;  //Hace 0 el counted bars para que siempre se calculen todos los valores.
    acumulado=0;
    A=0;
    ultimo=0;
    ultimoAnt=0;
    zPoints=0; //número total de picos +valles en el intérvalo
    CandleStart=0; //index del primer candlestick a dibujar
    cuentaEscritos=0; //para obtener el número de floats escritos en el archivo.
    error=0;
    if (PrimeraVez>0) return;     // garantiza que solo se ejecute una vez después de inicializacoón y inicializa buffers
    PrimeraVez++;
            // Inicio de variables estáticas
    inicStat=numEntradas*6;
    Print("Iniciado.");
    if (counted_bars==0 && downloadhistory) // history was downloaded
    {
        ArrayInitialize(MasterIndicator,0.0);
        ArrayInitialize(VentaE_S1,0.0);
        ArrayInitialize(CompraE_S1,0.0);
        ArrayInitialize(ZigzagBuffer,0.0);
        ArrayInitialize(HighMapBuffer,0.0);
        ArrayInitialize(LowMapBuffer,0.0);
    }
    counted_bars=iBars(NULL,PERIOD_M15);      // calcula el máximo posible de candlesticks a procesar.
    if (counted_bars>(iBars(NULL,PERIOD_M1)/15)) // si los de 1m/15 son menos los usa.
    {
        counted_bars=iBars(NULL,PERIOD_M1)/15;
    }
    if (numCandles15M_Max<counted_bars)//    // limita el número de candlesticks de 15m a procesar.
    {
        counted_bars=numCandles15M_Max+960;
    }
    // abre todos los archivos para escritura
    abrirArchivos();
    // escribe encabezado GTDv1
    if (error==0)
    {
          escribirGTDhdr();
    }
    else
    { 
      Print("Error 12: handle no asignado");
      return;
    }
    // calcula y guarda todas las entradas y el indicador de entrenamiento.
    calcularTodo(); 
    // cierra todos los archivos
    cerrarArchivos();
    Print("Terminado.");
    return(0);
}

