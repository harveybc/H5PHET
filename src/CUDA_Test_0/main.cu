#include <stdio.h>
#include <cuda.h>
//#include <stdio.h>
//#include <stdlib.h>
#define BLOCK_SIZE 16
// Es 1000 porque son 3 vectores de float =12kbytes, supuestamente tienen 16kbytes/multiproc (bloque)
#define THREADSPB 96
// kernel de CUDA para sumar dos vectores. (buffsize=tam copiado de global a shared)
__global__ void sumaVec (float* vec1_d, float* vec2_d, float* vec3_d, int numDatos)
{
    int i,j;
    float v1,v2,v3,acum; //registros para cada entrada y salida, acumulador par calculo de media
    // calcula el id del thread
    int idx=blockDim.x*blockIdx.x+threadIdx.x;
    // calcula el número de repeticiones para floor de todos los datos/buffSize
    int numReps=numDatos/THREADSPB;

	//THREADSPB también es el número de datos completos leidos de global a shared.
    // declara el vector shared para vec1,2 y 3 de tamaño buffsize/(3*sizeof(float))
    __shared__ float vec1_s[THREADSPB];
    __shared__ float vec2_s[THREADSPB];
    __shared__ float vec3_s[THREADSPB];
    // para i=0;i<numReps;i++ , para cada repetición
    acum=0;
    for (i=0;i<numReps;i++)
    {
		// coloca en memoria shared los valores de global, los lee en paralelo
		vec1_s[idx]=vec1_d[i*THREADSPB+idx];
		vec2_s[idx]=vec2_d[i*THREADSPB+idx];
		//sincroniza
		__syncthreads();
		// copia de shared a registros
        // para j=0;j<THREADSPB;j++
		for(j=0;j<THREADSPB;j++)
		{
            //sincroniza los threads para leer al mismo tiempo de la memoria shared.
            __syncthreads();
			// coloca el dato j en los registros necesarios
            v1=vec1_s[j];
            v2=vec2_s[j];
            // realiza los calculos con los registros.
            v3=v1*v2;
            // acumula valor de salida de entrenam y calculada para posterior calculo de error
            acum+=v3;
            // escribe en shared[j] los resultados
            vec3_s[j]=v3;
        }
        __syncthreads();
		// copia v3 parcial calculado de shared a global
		vec3_d[i*THREADSPB+idx]=vec3_s[idx];
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
    int numDatos=1000000;
	int tamVectores=numDatos*sizeof(float);
    // reserva memoria para vector 1,2 y 3 en host
    vec1_h=(float*)malloc(tamVectores);
    vec2_h=(float*)malloc(tamVectores);
    vec3_h=(float*)malloc(tamVectores);
    // llena con números aleatorios los vectores 1 y 2 en host
    printf("Generando vectores\n");
    for (i=0;i<numDatos;i++)
    {
        vec1_h[i]=(float)rand()/RAND_MAX;
        vec2_h[i]=(float)rand()/RAND_MAX;
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
    int threadsPorBlock=THREADSPB;
    int blocksPorGrid=(1+threadsPorBlock-1)/threadsPorBlock;
    sumaVec<<<blocksPorGrid,threadsPorBlock>>>(vec1_d,vec2_d,vec3_d,numDatos);
    // copia el vector 3 de la memoria global de device a host.
    cudaMemcpy(vec3_h, vec3_d, tamVectores, cudaMemcpyDeviceToHost);
    // imprime los resultados
    j=0;
    if (numDatos>1000)
        j=numDatos-1000;
    for (i=j;i<numDatos;i++)
    {
        printf("[%i]  %3.3f + %3.3f = %3.3f\n",i,vec1_h[i],vec2_h[i],vec3_h[i]);
    }
	while (1)
	{
	    if ('n' == getchar())
		break;
	}
    return 0;
}

