cls
copy main.c main.cu
echo LISTO
nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin" -o test.exe main.cu
echo FUCCKYOU


