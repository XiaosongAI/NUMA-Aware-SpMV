#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <numa.h>
#include <pthread.h>
#include <sched.h>
#include <sys/time.h>
#include "sys/sysinfo.h"
#include "mmio_highlevel.h"
#define VALUE_TYPE float

#ifndef NTIMES
#define NTIMES 100
#endif

int i, j, k;

typedef struct thread
{
    int coreidx, *rowptr, *colidx, start, end, tasknum, m, nthreads, threadidx;
    VALUE_TYPE *value, *Y, *X;
} thread;

typedef struct sparse_mtx
{
    int rownum;
    int colnum;
    int nnznum;
    int *rowptr;
    int *colidx;
    VALUE_TYPE *val;
}sparse_mtx;

int readmtx(char *filename, sparse_mtx *sourcematrix)
{
    int isSymmetric;
    //printf ("filename = %s\n", filename);
    mmio_info(&sourcematrix->rownum, &sourcematrix->colnum, &sourcematrix->nnznum, &isSymmetric, filename);
    sourcematrix->rowptr = (int *)malloc((sourcematrix->rownum+1) * sizeof(int));
    sourcematrix->colidx = (int *)malloc(sourcematrix->nnznum * sizeof(int));
    sourcematrix->val = (VALUE_TYPE *)malloc(sourcematrix->nnznum * sizeof(VALUE_TYPE));
    mmio_data(sourcematrix->rowptr, sourcematrix->colidx, sourcematrix->val, filename);
    //printf("Matrix A is %i by %i, #nonzeros = %i\n", sourcematrix->rownum, sourcematrix->colnum, sourcematrix->nnznum);
    return 1;
}

void spmv_serial(sparse_mtx *mtx, VALUE_TYPE *vec, VALUE_TYPE *res)
{
    for (i = 0; i < mtx->rownum; i++)
    {
        res[i] = 0;
        for (j = mtx->rowptr[i]; j < mtx->rowptr[i+1]; j++)
        {
            res[i] += mtx->val[j] * vec[mtx->colidx[j]];
        }
    }
}

/*int *mtx_partition(sparse_mtx *matrix, int part)
{
    idx_t rownum = matrix->rownum;
    idx_t *rowptr = (idx_t *)malloc(sizeof(idx_t)*(matrix->rownum+1));
    for (i = 0; i <= rownum; i++) rowptr[i] = matrix->rowptr[i];
    idx_t *colidx = (idx_t *)malloc(sizeof(idx_t)*matrix->nnznum);
    for (i = 0; i < matrix->nnznum; i++) colidx[i] = matrix->colidx[i];
    idx_t nweights = 1;
    idx_t *reorder = (idx_t *)malloc(sizeof(idx_t)*matrix->rownum);
    idx_t objval;
    idx_t parts = part;
    METIS_PartGraphKway(&rownum, &nweights, rowptr, colidx, NULL, NULL, NULL, &parts, NULL, NULL, NULL, &objval, reorder);
    int *order = (int *)malloc(sizeof(int) * matrix->rownum);
    for (i = 0; i < matrix->rownum; i++) order[i] = reorder[i];

    return order;

    free(rowptr);
    free(colidx);
    free(reorder);
    free(order);
}*/

int mtx_reorder(sparse_mtx *matrix, int *reorder_res, int num_parts)
{
        sparse_mtx reordered_matrix;
        reordered_matrix.rownum = matrix->rownum;
        reordered_matrix.colnum = matrix->colnum;
        reordered_matrix.nnznum = matrix->nnznum;
        reordered_matrix.rowptr = (int *)malloc(sizeof(int) * (reordered_matrix.rownum+1));
        int * temp_rowptr = (int *)malloc(sizeof(int) * (reordered_matrix.rownum+1));
        reordered_matrix.colidx = (int *)malloc(sizeof(int) * reordered_matrix.nnznum);
        for (i = 0; i <= reordered_matrix.rownum; i++)
        {
            reordered_matrix.rowptr[i] = 0;
            temp_rowptr[i] = 0;
        }
        for (i = 0; i < reordered_matrix.nnznum; i++) reordered_matrix.colidx[i] = 0;
        reordered_matrix.val = (VALUE_TYPE *)malloc(sizeof(VALUE_TYPE) * reordered_matrix.nnznum);
        int *cnt = (int *)malloc(sizeof(int)*(num_parts+1));
        for (i = 0; i <= num_parts; i++) cnt[i] = 0;
        int *index = (int *)malloc(sizeof(int)*num_parts);
        for (i = 0; i < num_parts; i++) index[i] = 0;
        int *temp_reorder = (int *)malloc(sizeof(int)*matrix->rownum);
        for (i = 0; i < matrix->rownum; i++)
        {
            cnt[reorder_res[i]]++;
            temp_reorder[i] = 0;
        }
        exclusive_scan(cnt, num_parts+1);
        for (i = 0; i < num_parts; i++) index[i] = 0;
        for (i = 0; i < matrix->rownum; i++)
        {
            temp_reorder[i] = cnt[reorder_res[i]] + index[reorder_res[i]];
            index[reorder_res[i]] = index[reorder_res[i]] + 1;
        }
        
        for (i = 0; i < matrix->rownum; i++)
            reordered_matrix.rowptr[temp_reorder[i]] = matrix->rowptr[i+1] - matrix->rowptr[i];
        exclusive_scan(reordered_matrix.rowptr, (matrix->rownum+1));
        for (i = 0; i < matrix->rownum; i++)
        {
            int nnz_cnt = 0;
            for (j = matrix->rowptr[i]; j < matrix->rowptr[i+1]; j++)
            {
                reordered_matrix.colidx[reordered_matrix.rowptr[temp_reorder[i]]+nnz_cnt] = matrix->colidx[j];
                nnz_cnt += 1;
            }
        }

        for (i = 0; i < reordered_matrix.nnznum; i++)
            reordered_matrix.colidx[i] = temp_reorder[reordered_matrix.colidx[i]];
        for (i = 0; i < reordered_matrix.nnznum; i++) reordered_matrix.val[i] = 1.0;

        *matrix = reordered_matrix;

        free(cnt);
        free(index);
        free(temp_reorder);
        free(temp_rowptr);

        return 1;
}

int vec_reorder(VALUE_TYPE *vec, int *reorder, int num_parts, int length)
{
    int *cnt = (int *)malloc(sizeof(int)*(num_parts+1));
    for (i = 0; i <= num_parts; i++) cnt[i] = 0;
    int *index = (int *)malloc(sizeof(int)*num_parts);
    for (i = 0; i < num_parts; i++) index[i] = 0;
    int *temp_reorder = (int *)malloc(sizeof(int)*length);
    VALUE_TYPE *temp_vec = (VALUE_TYPE *)malloc(sizeof(VALUE_TYPE )*length);
    for (i = 0; i < length; i++)
    {
        cnt[reorder[i]]++;
        temp_reorder[i] = 0;
        temp_vec[i] = 0;
    }
    exclusive_scan(cnt, num_parts+1);
    for (i = 0; i < num_parts; i++) index[i] = 0;
    for (i = 0; i < length; i++)
    {
        temp_reorder[i] = cnt[reorder[i]] + index[reorder[i]];
        index[reorder[i]] = index[reorder[i]] + 1;
    }
    for (i = 0; i < length; i++)
        temp_vec[temp_reorder[i]] = vec[i];
    for (i = 0; i < length; i++)
      vec[i] = temp_vec[i];
}

void writeresults(char *filename_res, char *filename, int m, int n, int nnzR, double time_metis, double GFlops_metis, double time_metis_numa, double GFlops_metis_numa, int nthreads)
{
    FILE *fres = fopen(filename_res, "a");
    if (fres == NULL) printf("Writing results fails.\n");
    fprintf(fres, "%s, %i, %i, %d, %lf ms, %lf GFlops, %lf ms, %lf GFlops, %i\n", filename, m, n, nnzR, time_metis, GFlops_metis, time_metis_numa, GFlops_metis_numa, nthreads);
    fclose(fres);
}

void pthread_process (void *arg)
{
    thread *pn = (thread*)arg;
    for (int i = pn->start; i < pn->end; i++)
    {
        VALUE_TYPE sum = 0;
        for (int j = pn->rowptr[i]; j < pn->rowptr[i+1]; j++)
        {
            sum += pn->value[j] * pn->X[pn->colidx[j]];
        }
        pn->Y[i] = sum;
    }
}

int **subrowptr, **subcolidx;
VALUE_TYPE **subval, **X, **Y;

typedef struct numaspmv
{
    int nthreads, numanodes, m, coreidx, alloc, *subX_ex, *subX;
    VALUE_TYPE *value, **Y, **X;
} numaspmv;

void *spmv(void *arg){
    numaspmv *pn = (numaspmv*)arg;
    int me = pn->alloc;
    numa_run_on_node(me);
    int m = pn->m;
    int nthreads = pn->nthreads;
    int numanodes = pn->numanodes;
    int coreidx = pn->coreidx;
    int eachnumathreads = nthreads/numanodes;
    int task = ceil((double)m/(double)eachnumathreads);
    int start = coreidx*task;
    int end = (coreidx+1)*task>m?m:(coreidx+1)*task;
    VALUE_TYPE *val = subval[me];
    VALUE_TYPE *y = Y[me];
    int *rpt = subrowptr[me];
    int *col = subcolidx[me];
    for (int u = start; u < end; u++)
    {
        VALUE_TYPE sum = 0;
        for (int j = rpt[u]; j < rpt[u+1]; j++)
        {
            int Xpos = col[j]/pn->subX[0];
            int remainder = col[j]-pn->subX_ex[Xpos];
            sum += val[j] * X[Xpos][remainder];
        }
        y[u] = sum;
    }
}

int main(int argc, char ** argv)
{
    struct timeval t1, t2;
    char *filename=argv[1];
    sparse_mtx mtx;
    int numanodes = numa_max_node()+1;
    int PARTS = numanodes;
    //int partition_num = atoi(argv[2]);
    int nthreads = atoi(argv[2]);
    readmtx(filename, &mtx);
    for (i = 0; i < mtx.nnznum; i++)
    {
        mtx.val[i] = 1.0;
    }
    VALUE_TYPE *vector = (VALUE_TYPE *)malloc(sizeof(VALUE_TYPE)*mtx.colnum);
    memset(vector, 0, sizeof(VALUE_TYPE)*mtx.colnum);
    VALUE_TYPE *y_golden = (VALUE_TYPE *)malloc(sizeof(VALUE_TYPE) * mtx.colnum);
    VALUE_TYPE *y_pthread = (VALUE_TYPE *)malloc(sizeof(VALUE_TYPE) * mtx.colnum);

    for (i = 0; i < mtx.colnum; i++)
    {
        vector[i] = rand()%10+1;
    }
    spmv_serial(&mtx, vector, y_golden);

    int task = ceil((double)mtx.rownum/(double)nthreads);
    thread *th = (thread*)malloc(nthreads*sizeof(thread));
    for(i = 0; i < nthreads; i++)
    {
        th[i].start = i * task;
        th[i].end = (i+1)*task>mtx.rownum?mtx.rownum:(i+1)*task;
        th[i].tasknum = i;
        th[i].nthreads = nthreads;
        th[i].value = mtx.val;
        th[i].m = mtx.rownum;
        th[i].X = vector;
        th[i].Y = y_pthread;
        th[i].rowptr = mtx.rowptr;
        th[i].colidx = mtx.colidx;
    }
    pthread_t* threads_metis=(pthread_t *)malloc(nthreads*sizeof(pthread_t));
    pthread_attr_t pthread_custom_attr_metis;
    pthread_attr_init(&pthread_custom_attr_metis);
    gettimeofday(&t1, NULL);
    for(int r = 0; r < NTIMES; r++)
    {
        for (i = 0; i < nthreads; i++)
        {
            pthread_create(&threads_metis[i], &pthread_custom_attr_metis, (void*)pthread_process, (void *)(th+i));
        }
        for (i = 0; i < nthreads; i++)
        {
            pthread_join(threads_metis[i], NULL);
        }
    }
    gettimeofday(&t2, NULL);
    double time_pthread = ((t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0) / NTIMES;
    double GFlops_pthread = 2 * mtx.nnznum / time_pthread / pow(10,6);
    printf("pthread_metis time %.2f ms GFlops %.2f \n", time_pthread, GFlops_pthread);

    /***************************************************************************************/
    /*memset(y_pthread, 0, sizeof(VALUE_TYPE)*mtx.rownum);
    int *reorder = (int *)malloc(sizeof(int)*mtx.rownum);
    reorder = mtx_partition(&mtx, partition_num);
    mtx_reorder(&mtx, reorder, partition_num);
    vec_reorder(vector, reorder, partition_num, mtx.colnum);

    gettimeofday(&t1, NULL);
    for(int r = 0; r < NTIMES; r++)
    {
        for (i = 0; i < nthreads; i++)
        {
            pthread_create(&threads_metis[i], &pthread_custom_attr_metis, (void*)pthread_process, (void *)(th+i));
        }
        for (i = 0; i < nthreads; i++)
        {
            pthread_join(threads_metis[i], NULL);
        }
    }
    gettimeofday(&t2, NULL);
    double time_pthread_metis = ((t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0) / NTIMES;
    double GFlops_pthread_metis = 2 * mtx.nnznum / time_pthread_metis / pow(10,6);
    printf("pthread_metis time %.2f ms GFlops %.2f \n", time_pthread_metis, GFlops_pthread_metis);*/


    /***************************************************************************************/

    int eachnumacores = nthreads/numanodes;
    int m = mtx.rownum;
    int n = mtx.colnum;
    int nnz = mtx.nnznum;
    int *rowptr = mtx.rowptr;
    int *colidx = mtx.colidx;
    VALUE_TYPE* val = mtx.val;
    if (numa_available()<0)
    {
        printf("Your system does not support NUMA API\n");
        return 0;
    }

    int *subrowpos = (int *)malloc(sizeof(int) * (PARTS+1));
    int *subm = (int *)malloc(sizeof(int) * PARTS);
    int *subm_ex = (int *)malloc(sizeof(int) * (PARTS+1));
    int *subX = (int *)malloc(sizeof(int) * PARTS);
    int *subX_ex = (int *)malloc(sizeof(int) * (PARTS+1));
    int *subnnz = (int *)malloc(sizeof(int) * PARTS);
    int *subnnz_ex = (int *)malloc(sizeof(int) * (PARTS+1));
    for (i = 0; i < PARTS; i++)
        subrowpos[i] = (ceil((double)m/(double)PARTS))*i>m?m:((ceil((double)m/(double)PARTS))*i);
    subrowpos[PARTS] = m;
    for (i = 0; i < (PARTS-1); i++)
    {
        subX[i] = ceil((double)n/(double)PARTS);
    }
    subX[PARTS-1] = n - subX[0]*(PARTS-1);
    for (i = 0; i < PARTS; i++)
    {
        subX_ex[i] = subX[i];
    }
    exclusive_scan(subX_ex, PARTS+1);
    for (i = 0; i < PARTS; i++)
    {
        subm[i] = subrowpos[i+1] - subrowpos[i];
        subm_ex[i] = subrowpos[i+1] - subrowpos[i];
    }
    exclusive_scan(subm_ex, PARTS+1);
    for (i = 0; i < PARTS; i++)
    {
        subnnz[i] = rowptr[subrowpos[i]+subm[i]]-rowptr[subrowpos[i]];
        subnnz_ex[i] = subnnz[i];
    }
    exclusive_scan(subnnz_ex, PARTS+1);
    numaspmv *p = (numaspmv*)malloc(nthreads*sizeof(numaspmv));
    pthread_t* threads=(pthread_t *)malloc(nthreads*sizeof(pthread_t));
    pthread_attr_t pthread_custom_attr;
    pthread_attr_init(&pthread_custom_attr);
    for(i = 0; i < nthreads; i++)
    {
        p[i].alloc = i%numanodes;
        p[i].numanodes = numanodes;
        p[i].nthreads = nthreads;
        p[i].subX = subX;
        p[i].subX_ex = subX_ex;
    }
    for(i = 0; i < eachnumacores; i++)
    {
        for(j = 0; j < numanodes; j++)
        {
            p[i*numanodes+j].coreidx = i;
            p[i*numanodes+j].m = subm[j];
        }
    }
    subrowptr = (int **)malloc(sizeof(int *)*nthreads);
    subcolidx = (int **)malloc(sizeof(int *)*nthreads);
    subval = (VALUE_TYPE **)malloc(sizeof(VALUE_TYPE *)*nthreads);
    X = (VALUE_TYPE **)malloc(sizeof(VALUE_TYPE *)*nthreads);
    Y = (VALUE_TYPE **)malloc(sizeof(VALUE_TYPE *)*nthreads);

    for(i = 0; i < nthreads; i++)
    {
        subrowptr[i] = numa_alloc_onnode(sizeof(int)*(subm[p[i].alloc]+1), p[i].alloc);
        subcolidx[i] = numa_alloc_onnode(sizeof(int)*subnnz[p[i].alloc], p[i].alloc);
        subval[i] = numa_alloc_onnode(sizeof(VALUE_TYPE)*subnnz[p[i].alloc], p[i].alloc);
        X[i] = numa_alloc_onnode(sizeof(VALUE_TYPE)*subX[p[i].alloc], p[i].alloc);
        Y[i] = numa_alloc_onnode(sizeof(VALUE_TYPE)*subm[p[i].alloc], p[i].alloc);
    }
    int curr_core = 0;
    for(i = 0; i < numanodes; i++)
    {
        for(j = 0; j < eachnumacores; j++)
        {
            for (k = 0; k <= subm[i]; k++)
            {
            	curr_core = i+j*eachnumacores;
            	if (curr_core < nthreads)
            	{
                	subrowptr[i+j*eachnumacores][k] = rowptr[subrowpos[i]+k];
            	}
            }
        }
        for(j = 0; j < eachnumacores; j++)
        {
            for (k = 0; k < subnnz[i]; k++)
            {
            	curr_core = i+j*eachnumacores;
            	if (curr_core < nthreads)
            	{
	                subcolidx[i+j*eachnumacores][k] = colidx[subnnz_ex[i]+k];
	                subval[i+j*eachnumacores][k] = val[subnnz_ex[i]+k];
	            }
            }
        }
        for(j = 0; j < eachnumacores; j++)
        {
            for (k = 0; k < subX[i]; k++)
            {
            	curr_core = i+j*eachnumacores;
            	if (curr_core < nthreads)
            	{
                    X[i+j*eachnumacores][k] = vector[subX_ex[i]+k];
            	}
            }
        }

    }
    for (i = 0; i < nthreads; i++)
    {
        if (i%numanodes != 0)
        {
            int temprpt = subrowptr[i][0];
            for (j = 0; j <= subm[i%numanodes]; j++)
            {
                subrowptr[i][j] -= temprpt;
            }
        }
    }

    gettimeofday(&t1, NULL);
    for(int r = 0; r < NTIMES; r++)
    {
        for (i = 0; i < nthreads; i++)
        {
            pthread_create(&threads[i], &pthread_custom_attr, spmv, (void *)(p+i));
        }
        for (i = 0; i < nthreads; i++)
        {
            pthread_join(threads[i], NULL);
        }
    }
    gettimeofday(&t2, NULL);
    double time_numa = ((t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0) / NTIMES;
    double GFlops_numaspmv = 2 * nnz / time_numa / pow(10,6);
    double bandwidth = (((m+1)+nnz)*sizeof(int)+(2*nnz+m)*sizeof(VALUE_TYPE))*nthreads/time_numa/pow(10,6);
    printf("numaspmv time %.2f ms  GFlops_numaspmv %.2f \n", time_numa, GFlops_numaspmv);

    VALUE_TYPE *Y_gather = (VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*m);
    memset(Y_gather,0,sizeof(VALUE_TYPE)*m);
    VALUE_TYPE *Y_golden = (VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*m);
    memset(Y_golden,0,sizeof(VALUE_TYPE)*m);
    for (i = 0; i < PARTS; i++)
    {
        for (j = 0; j < subm[i]; j++)
        {
            Y_gather[subm_ex[i]+j] = Y[i][j];
        }
    }

    //vec_reorder(y_golden, reorder, partition_num, mtx.colnum);
    
    int cnt = 0;
    for (i = 0; i < mtx.colnum; i++)
    if (Y_gather[i] != y_golden[i])
        cnt++;
    //printf("error count %d\n", cnt);
    writeresults("pthread_metis_numa.csv", filename, m, n, nnz, time_pthread, GFlops_pthread, time_numa, GFlops_numaspmv, nthreads);

}
