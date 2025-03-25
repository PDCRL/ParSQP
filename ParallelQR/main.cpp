#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <pthread.h>
#include "include/bn2.h"
#include <tbb/concurrent_priority_queue.h>
#include <tbb/tbb.h>
#include <unistd.h>
#include <csignal>
#include <cstdlib>
#include <mutex>

#define NUM_THREADS 26

#define BETA 11
#define ALPHA 11
#define BETA_DIV_ALPHA ((int)BETA / (int)ALPHA)

#define USE_PRIORITY_MAIN_QUEUE 0
typedef struct
{
    int tid;
    int total_task_rows;
    int total_task_cols;
    int m;
    int n;
    double *mat;
} thread_args_ts;

std::vector<std::stringstream> logstreams(NUM_THREADS);

TaskTable task_table;
DependencyTableAtomic dependency_table;

std::vector<double> global_up_array, global_b_array;


struct TaskComparator {
    bool operator()(const Task* a, const Task* b) const {
        return a->priority < b->priority; 
    }
};

tbb::concurrent_queue<Task *> wait_queue;//, taskPQ;

#if USE_PRIORITY_MAIN_QUEUE
tbb::concurrent_priority_queue<Task *, TaskComparator> taskPQ;
#else
tbb::concurrent_queue<Task *> taskPQ;
#endif

void complete_task1(double *&mat, int m, int n, int row_start, int row_end, int col_start, int col_end)
{

    double sm, sm1, cl, clinv, up, b;
    int _row_start = row_start == 1 ? 0 : row_start;

    for (int lpivot = _row_start; lpivot < row_end; lpivot++)
    {
        cl = fabs(mat[lpivot * n + lpivot]);
        sm1 = 0;

        for (int k = lpivot + 1; k < n; k++)
        {
            sm = fabs(mat[lpivot * n + k]);
            sm1 += sm * sm;
            cl = fmax(sm, cl);
        }

        if (cl <= 0.0)
        {
            continue;
        }
        clinv = 1.0 / cl;

        double d__1 = mat[lpivot * n + lpivot] * clinv;
        sm = d__1 * d__1;
        sm += sm1 * clinv * clinv;

        cl *= sqrt(sm);

        if (mat[lpivot * n + lpivot] > 0.0)
        {
            cl = -cl;
        }

        up = mat[lpivot * n + lpivot] - cl;
        mat[lpivot * n + lpivot] = cl;

        if (row_end - lpivot < 0)
        {
            continue;
        }

        b = up * mat[lpivot * n + lpivot];

        if (b >= 0.0)
        {
            continue;
        }

        b = 1.0 / b;

        global_up_array[lpivot] = up;
        global_b_array[lpivot] = b;

        for (int j = lpivot + 1; j < col_end; j++)
        {
            sm = mat[j * n + lpivot] * up;

            for (int i__ = lpivot + 1; i__ < n; i__++)
            {
                sm += mat[j * n + i__] * mat[lpivot * n + i__];
            }

            if (sm == 0.0)
            {
                continue;
            }

            sm *= b;
            mat[j * n + lpivot] += sm * up;

            for (int i__ = lpivot + 1; i__ < n; i__++)
            {
                mat[j * n + i__] += sm * mat[lpivot * n + i__];
            }
        }
    }
}

void complete_task2(double *&mat, int m, int n, int row_start, int row_end, int col_start, int col_end)
{

    int _row_start = row_start == 1 ? 0 : row_start;
    int _col_start = col_start == 1 ? 0 : col_start;

    double up = 0.0, b = 0.0, sm = 0.0;

    for (int lpivot = _row_start; lpivot < row_end; lpivot++)
    {
        up = global_up_array[lpivot];
        b = global_b_array[lpivot];

        for (int j = _col_start; j < col_end; j++)
        {
            sm = mat[j * n + lpivot] * up;

            for (int i__ = lpivot + 1; i__ < n; i__++)
            {
                sm += mat[j * n + i__] * mat[lpivot * n + i__];
            }

            if (sm == 0.0)
            {
                continue;
            }

            sm *= b;
            mat[j * n + lpivot] += sm * up;

            for (int i__ = lpivot + 1; i__ < n; i__++)
            {
                mat[j * n + i__] += sm * mat[lpivot * n + i__];
            }
        }
    }
}

void *thdwork(void *params)
{
    thread_args_ts *thread_args = (thread_args_ts *)params;
    int total_task_rows = thread_args->total_task_rows;
    int total_task_cols = thread_args->total_task_cols;
    double *mat = thread_args->mat;
    int m = thread_args->m;
    int n = thread_args->n;
    
    while (1)
    {
        Task *new_task = nullptr;
        if (taskPQ.try_pop(new_task)) 
        {
            int i = new_task->chunk_idx_i;
            int j = new_task->chunk_idx_j;

            int row_start = new_task->row_start;
            int row_end = new_task->row_end;
            int col_start = new_task->col_start;
            int col_end = new_task->col_end;

            if (new_task->type == 1)
            {
                complete_task1(mat, m, n, row_start, row_end, col_start, col_end);
                dependency_table.setDependency(i, j, true);
                for (int k = i + 1; k < total_task_rows; k++)
                {
                    Task *next_task = task_table.getTask(k, j);

                    if (j == 0 || dependency_table.getDependency(k, j - 1))
                    {
                        taskPQ.push(next_task);
                    }
                    else
                    {
                        wait_queue.push(next_task);
                    }
                }
            }
            else if (new_task->type == 2)
            {
                complete_task2(mat, m, n, row_start, row_end, col_start, col_end);
                dependency_table.setDependency(i, j, true);
                if (new_task->enq_nxt_t1 && (j + 1) <= total_task_cols)
                {
                    taskPQ.push(task_table.getTask((j + 1) / BETA_DIV_ALPHA, j + 1));  
                }
            }
        }
    
        Task *local_task = nullptr;
        if (wait_queue.try_pop(local_task)) 
        {
            int i = local_task->chunk_idx_i;
            int j = local_task->chunk_idx_j;
            if (dependency_table.getDependency(i, j - 1))
            {
                taskPQ.push(local_task);
            }
            else
            {
                wait_queue.push(local_task);
            }
        }

        if (dependency_table.getDependency(total_task_rows - 1, BETA_DIV_ALPHA * (total_task_rows - 1)))
        {
            break;
        }
    }

    return nullptr;
}

void save(double *&mat , int mat_rows, int mat_cols, const std::string &filename)
{
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        throw std::runtime_error("Error opening file for writing: " + filename);
    }

    for (int i = 0; i < mat_rows; ++i)
    {
        for (int j = 0; j < mat_cols; ++j)
        {
            outfile << mat[i * mat_cols + j];
            if (j < mat_cols - 1)
            {
                outfile << " ";
            }
        }
        outfile << "\n";
        if (outfile.fail())
        {
            throw std::runtime_error("Error writing matrix data to file: " + filename);
        }
    }
    outfile.close();
}

void dagqrf(double *mat , int mat_rows, int mat_cols, std::vector<double> &up_array)
{
    int total_task_rows = std::ceil((double)mat_rows / BETA);
    int total_task_cols = std::ceil((double)mat_rows / ALPHA);

    global_up_array.resize(mat_rows, 0.0);
    global_b_array.resize(mat_rows, 0.0);

    dependency_table.init(total_task_rows, total_task_cols);
    task_table.init(total_task_rows, total_task_cols, ALPHA, BETA, mat_rows, mat_cols);

    std::vector<pthread_t> threads(NUM_THREADS);
    std::vector<thread_args_ts> thread_args(NUM_THREADS);

    for (int i = 0; i < NUM_THREADS; i++)
    {
        thread_args[i].tid = i;
        thread_args[i].total_task_rows = total_task_rows;
        thread_args[i].total_task_cols = total_task_cols;
        thread_args[i].m = mat_rows;
        thread_args[i].n = mat_cols;
        thread_args[i].mat = mat;
    }
    taskPQ.push(task_table.getTask(0, 0));
    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_create(&threads[i], NULL, thdwork, &thread_args[i]);
    }

    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
    }
    up_array = global_up_array;
    //save(mat,mat_rows,mat_cols,"test_output.txt");
}
