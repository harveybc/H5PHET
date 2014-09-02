//#include <stdio.h>
#include <cuda.h>
//#include <stdio.h>
//#include <stdlib.h>
#define BLOCK_SIZE 16
#define TAMVECTOR 4000
// kernel de CUDA para sumar dos vectores. (buffsize=tam copiado de global a shared)
__global__ void sumaVec (float* vec1_d, float* vec2_d, float* vec3_d, int numDatos, int buffSize)
{
    int i,j;
    int a,b,c, acum;
    // calcula el id del thread
    int idX=blockDim.x*blockIdx.x+threadIdx.x;
    // calcula el número de repeticiones para floor de todos los datos/buffSize
    int numReps=numDatos/TAMVECTOR;
    // declara el vector shared para vec1,2 y 3 de tamaño buffsize/(3*sizeof(float))
    int tamVector=TAMVECTOR;
    __shared__ float vec1_s[TAMVECTOR];
    __shared__ float vec2_s[TAMVECTOR];
    __shared__ float vec3_s[TAMVECTOR];
    // para i=0;i<numReps;i++ , para cada repetición
    acum=0;
    for (i=0;i<numReps;i++)
    {
        // si el id==0 (para que solo lo haga un core por block)   NECESARIO?
        if (idX==0)
        {
            // copia de global vec1[i*buffsize], buffsize datos a shared.4=sizeof(float)
            // FALTA: colocar los datos consecutivos en bancos diferentes para optimizar lectura
            cudaMemcpy(vec1_s,vec1_d+(i*tamVector),tamVector*4,cudaMemcpyDeviceToDevice);
            cudaMemcpy(vec2_s,vec2_d+(i*tamVector),tamVector*4,cudaMemcpyDeviceToDevice));
        }
        //sincroniza : FALTA sincronixar lecturas?
        __syncthreads();
        // para j=0;j<tamVector;j++
        for (j=0;j<tamVector;j++)
        {
            // coloca el dato j en los registros necesarios
            a=vec1_s[j];
            b=vec2_s[j];
            // realiza los calculos con los registros.
            c=a+b;
            // acumula valor de salida de entrenam y calculada para posterior calculo de error
            acum+=c;
            // escribe en shared[j] los resultados
            vec3_s[j]=c;
            //sincroniza los threads para leer al mismo tiempo de la memoria shared.
            __syncthreads();
        }
        // si el id==0 (para que solo lo haga un core por block)   NECESARIO?
        if (idX==0)
        {
            // copia el vector de resultados de shared a global.
            // FALTA: colocar los datos consecutivos en bancos diferentes para optimizar escritura?
            cudaMemcpy(vec3_d,vec3_s+(i*tamVector),tamVector*4,,cudaMemcpyDeviceToDevice));
        }
        __syncthreads();
    }
    // calcula el fitness
    // para i=0;i<numReps;i++
        // si el id==0 (para que solo lo haga un core por block)   NECESARIO?
            // copia de global vec1[i*buffsize], buffsize datos a shared.
        // para j=0;j<buffsize;j++
        // calcula las 3 sumatorias del coeficiente de correlación
}


int main()
// prueba de CUDA: suma dos vectores componente a componente
{
    int i,j;
    float* vec1_h;
    float* vec2_h;
    float* vec3_h;
    float* vec1_d;
    float* vec2_d;
    float* vec3_d;
    int tamVectores=numDatos*sizeof(float);
    int numDatos=1000000;
    int buffSize=12000;
    // reserva memoria para vector 1,2 y 3 en host
    vec1_h=(float*)malloc(tamVectores);
    vec2_h=(float*)malloc(tamVectores);
    vec3_h=(float*)malloc(tamVectores);
    // llena con números aleatorios los vectores 1 y 2 en host
    printf("Generando vectores\n");
    for (i=0;i<numDatos;i++)
    {
        vec1_h=(float)rand()/RAND_MAX;
        vec2_h=(float)rand()/RAND_MAX;
    }
    printf("Reservando memoria en device\n");
    // reserva memoria para vector 1,2 y 3 en device
    cudaMalloc(&vec1_d, tamVectores); //TODO: POSIBLE PROBLEMA: SE MACE MALLOC DE &vec1_d en lugar de vec1_d (es puntero)
    cudaMalloc(&vec2_d, tamVectores); //TODO: POSIBLE PROBLEMA: SE MACE MALLOC DE &vec1_d en lugar de vec1_d (es puntero)
    cudaMalloc(&vec3_d, tamVectores); //TODO: POSIBLE PROBLEMA: SE MACE MALLOC DE &vec1_d en lugar de vec1_d (es puntero)
    printf("Copiando vectores de host a device\n");
    // copia en la memoria global de device los vec 1 y 2 de host
    cudaMemcpy(vec1_d, vec1_h, tamVectores, cudaMemcpyHostToDevice);
    cudaMemcpy(vec2_d, vec2_h, tamVectores, cudaMemcpyHostToDevice);
    // ejecuta kernel de CUDA con los 3 vectores device, su tamaño y ShBuffSize como parámetros
    printf("Ejecutando kernel de CUDA para todos los threads\n");
    int threadsPorBlock=256;
    int blocksPorGrid=(N+threadsPorBlock-1)/threadsPorBlock;
    sumaVec<<<blocksPorGrid,threadsPorBlock>>>(vec1_d,vec2_d,vec3_d,numDatos,buffSize);
    // copia el vector 3 de la memoria global de device a host.
    cudaMemcpy(vec3_h, vec1_d, tamVectores, cudaMemcpyDeviceToHost);
    // imprime los resultados
    j=0;
    if (numDatos>1000)
        j=numDatos-1000;
    for (i=j;i<numDatos;i++)
    {
        printf("[%i]  %3.3f + %3.3f = %3.3f\n",i,vec1_h,vec2_h,vec3_h);
    }
    return 0;
}

