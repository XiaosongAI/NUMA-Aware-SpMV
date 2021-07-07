# NUMA-Aware-SpMV
稀疏矩阵-向量乘的NUMA感知优化（c/pthread）
基本环境：
  1. 需安装超图分割工具metis。
  2. 需安装libnuma和numactl。
编译：
  gcc -o pthread_metis_numa pthread_metis_numa.c -lnuma -lpthread -lm -O3
运行：
  ./pthread_metis_numa 矩阵名 线程数（例：./pthread_metis_numa M6.mtx 16）
