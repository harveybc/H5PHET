// H5PHET - FX_GTDGen
// Genera archivo GTD (Generic Training Data)de entrenamiento con datos de forex
// Se deben configurar manualmente los parámetros del GTD de salida.

// estructura del encabezado de SNNv1 (Simple Neural Network).
// formato SNNv1: [encabezado][unsigned conexIn[numConex]][unsigned conexOut[numconex]][char enabled[numconex]][float/double peso[numConex]]

#property indicator_separate_window
#property indicator_buffers 1  // actual y variación por cada uno de los tres TF: 1M, 15M y 4H
                               // si se va a hacer la salida, debe colocarse en 1 y el buffer
                               // de la salida debe tener index 0 en inicialización
#property indicator_color1 Magenta
//---- indicator parameters
//Agregado por Harvey
extern int numDias=105; // Días a calcular
extern int skipInicio=7; //  Dias para saltar desde el tiempo actual
extern int skipFinal=7; // Dias que deben contener datos posteriores a los datos leidos. 
extern double A_Sigma=2.435; // factor usado en sigma out=exp(-A(In*In))   de las fdt de las neuronas.
extern int pipsTPSL=1000; // nuero de pips sobre los que se calcula el masterTrainer2 (default 400)
extern int maxDias=6; // periodos en el futuro para calcular el masterTrainer (4 dias)

//---- indicator buffers
double SNN_Indicator[];//Para el indicador de red neuronal
double MasterIndicator[];  //indicador compuesto 0 = máxima ganancia en venta, 1= máx en compra
double CompraE_S1[];  //es 1 si el valor es el mínimo posible para cada zzorder y 0 si es el max.
double VentaE_S1[];  //para venta
double hOBV_buffer[];
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
double valorAnt[1000];
double normDMax[1000]; // usado para la normalización dinámica de las entradas max 1k.
double normDMin[1000]; // usado para la normalización dinámica de las entradas max 1k.
extern double decNormD=0.00000001; //(0.000001)decremento pir shift de la norma dinámica.0.001=4.1770,0.0001=9.7558,0.00001=10.8345,0.000001=11.0309,0.0000001=11.0901
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
    SetIndexBuffer(2,hOBV_buffer);
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
    IndicatorShortName("H5PHET - Master Trainer 2 ("+pipsTPSL+","+maxDias+","+numDias+")");
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
   numIns2 = inicStat;
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
        // FIN DE calculo de entradas estándar..
        // Para cargar ANN desde SNN máximo 10000 conexiones y 1000 entradas
        i15=i/15;
        i240=i/240;
        for (j=0;j<numEntradas;j++)
        {
            // llama la función calcularIndicador (const Moneda, const TF, shift, double param1, param2, param3,param4, param5)
            // par alos 3 timeframes y 3 diferenciales:
            tmpj=j*6;            
            // tf 1M
            valor[tmpj]=calcularIndicador("EURUSD",PERIOD_M1,j,i);
            // diferencial
            valor[tmpj+1]=valor[tmpj]-calcularIndicador("EURUSD",PERIOD_M1,j, i+1);
            //ft 15M es el anterior + diferencial/15
            valor[tmpj+2]=valorAnt[tmpj+2]+(valor[tmpj+1]/15);
            //diferencial igual a anterior+(diferencial1m/15)-diferencial1m(i+15)
            valor[tmpj+3]=valorAnt[tmpj+3]+(valor[tmpj+1]/15)-(calcularIndicador("EURUSD",PERIOD_M1,j, i+15)-calcularIndicador("EURUSD",PERIOD_M1,j, i+16));
            //tf 4H es el anterior + diferencial/240
            valor[tmpj+4]=valorAnt[tmpj+4]+(valor[tmpj+1]/240);
            //diferencial igual a anterior+(dif1m/240)-dif1m(i+240)
            valor[tmpj+5]=valorAnt[tmpj+5]+(valor[tmpj+1]/240)-(calcularIndicador("EURUSD",PERIOD_M1,j, i+240)-calcularIndicador("EURUSD",PERIOD_M1,j, i+241));
            // para las 6, guarda en anterior y normaliza
            // para las 6 normaliza
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
                // guarda el valor anterior
                valorAnt[tmpj+k]=valor[tmpj+k];
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
    else if (indicNum==29) // cloae
        return(iClose(symbol,timeframe,shift));
    else if (indicNum==30) // Tamaño de candlestick
        return(iHigh(symbol,timeframe,shift)-iLow(symbol,timeframe,shift));
    else if (indicNum==31) // Color y tamaño de cuerpo
        return(iClose(symbol,timeframe,shift)-iOpen(symbol,timeframe,shift));
    else if (indicNum==32)//  Entrada HOBV(i) = HOBV(i+1)+(color[i]*iVolume(i)) solo funciona en  tf=1m
    {
        hOBV_buffer[shift]=hOBV_buffer[shift+1]+((iClose(symbol,timeframe,shift)-iOpen(symbol,timeframe,shift))*iVolume(symbol,timeframe,shift));
        return(hOBV_buffer[shift]);
    }
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
    if (i<maxMinutos)
    {
        Print("Error, i<maxMinutos");
        return(0);
    }
    for (j=i; j>(i-maxMinutos) ; j--)
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
    for (j=0;j<(inicStat+numSalidas);j++)
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
    double tmp2;
    datetime tmpFecha;
    // Para cargar ANN desde SNN máximo 10000 conexiones
    int numIns,numBias,numOuts,numHiddens,numConex;
    int conexIn[10000];
    int conexOut[10000];
    double peso[10000];
    double tmpPeso;
   
// variable de valores usada para calcular la salida de cada nodo. las neuronas de entrada y bias no se evaluan, sus valores
// de salida se cargan directamente de las entradas de la red neuronal (4 Entradas de MasterTrainer).
// maximo 9000 neuronas.
    // inicializa los valores en 0
    ArrayInitialize(valor,0.0);
    ArrayInitialize(valorAnt,0.0);
    ArrayInitialize(valorCalculado,0.0);
    ArrayInitialize(normDMax,-100000.0);
    ArrayInitialize(normDMin,100000.0);
    ArrayInitialize(hOBV_buffer,0);    
    // INICIO de barrido por candlesticks
    tmi=iBars("EURUSD",PERIOD_M1); //el número total de candlesticks
    // verifica límites
    if (tmi<(1+(numDias*1440)+(skipInicio*1440)+(skipInicio*1440)))
    {
      Print("Error calculando las hOBV_buffer");   
      return;
    }
    //establece los límites de tiempo para exportar tomando en cuenta los skip en dias
    for (i=((numDias*1440)+(skipInicio*1440)); i>=(skipInicio*1440); i--)
    {
        // calcula todas las entradas par ael candlestick i
        calcularEntradas(i);
        // calcula el MI.
        MasterIndicator[i] = -calcularMI(i); //EL INDICADOR Quedó al revés    
  /* FALLA CUANDO SE NORMALIZA.             
                // Normaliza dinámicamente el Master Indicator
                normDMax[inicStat]+=(decNormD*(normDMax[inicStat]-normDMin[inicStat]));
                normDMin[inicStat]-=(decNormD*(normDMax[inicStat]-normDMin[inicStat]));
                // compara normas con cada valor
                if (MasterIndicator[i]>normDMax[inicStat])
                    normDMax[inicStat]=MasterIndicator[i];
                if (MasterIndicator[i]<normDMin[inicStat])
                    normDMin[inicStat]=MasterIndicator[i];
                // verifica que las normas tengan datos válidos
                if (normDMax[inicStat]<normDMin[inicStat])
                {   
                    // los invierte
                    tmp2=normDMax[inicStat];
                    normDMax[inicStat]=normDMin[inicStat];
                    normDMin[inicStat]=tmp2;
                }
                // normaliza cada valor (Máximo=1, mínimo=-1).
                if (normDMax[inicStat]!=normDMin[inicStat])
                {
                    MasterIndicator[i]=(MasterIndicator[i]-normDMin[inicStat])/(normDMax[inicStat]-normDMin[inicStat]);  //para 0,1
                    MasterIndicator[i]=((2*MasterIndicator[i])-1); //para -1,1

                }            
*/

        // guarda MI en posición de salida.
        valor[inicStat]=MasterIndicator[i];
        // escribe los valores en el archivo GTD.
        escribirTodo(); //guarda todos los valores en los archivos de salida.
        //imprime la info de cada periodo.
        tmpFecha=iTime(NULL,PERIOD_M1,i);
        Print("Fecha: ", TimeYear(tmpFecha),"-",TimeMonth(tmpFecha),"-",TimeDay(tmpFecha)," ",TimeHour(tmpFecha),":",TimeMinute(tmpFecha),"  Candlestick ",i);
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
    numEntradas=33;
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
    ArrayInitialize(MasterIndicator,0.0);
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

