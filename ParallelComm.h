/* Common declarations for the functions in comm.c */
#ifndef KRIPKE_COMM_H__
#define KRIPKE_COMM_H__

#include<vector>
#include<mpi.h>
#include "Task.h"

class ParallelComm {
  public:
    ParallelComm();
    virtual ~ParallelComm();

    // Adds a task to the work queue
    virtual void addTask(int sdom_id, Task &task);

    // Checks if there are any outstanding task to complete
    // false indicates all work is done, and all sends have completed
    virtual bool workRemaining(void);

    // Returns a vector of ready tasks, and clears them from the ready queue
    virtual std::vector<int> readyTasks(void);

    // Marks task as complete, and performs downwind communication
    virtual void markComplete(int task_id);

    // Returns when the specified task is ready
    virtual void WaitforReady(int task_id);

    // Sets the task size
    virtual void setSizes(int cells, int angles, int groups);
    
    // Sets the All_Tasks vector so we can get access to downstream task
    virtual void SetTask(Task* Task);

  protected:
    int computeTag(int mpi_rank, int task_id);
    static void computeRankTask(int tag, int &mpi_rank, int &task_id);

    int findTaskinTasks(int task_id);
    int findTask(int task_id);
    Task *dequeueTask(int task_id);
    void postRecvs(int task_id, Task &task);
    void postSends(Task *task, double *buffers[3]);
    void testRecieves(void);
    void waitAllSends(void);
    std::vector<int> getReadyList(void);


    // These vectors contian the recieve requests and the index in the request vector
    std::vector<MPI_Request> recv_requests;
    std::vector<int> recv_tasks;
    std::vector<int> request_index;

    // These vectors have the tasks, and the remaining dependencies
    std::vector<int> queue_task_ids;
    std::vector<Task *> queue_tasks;
    std::vector<int> queue_depends;

    // These vectors have the remaining send requests that are incomplete
    std::vector<MPI_Request> send_requests;

    // These are the sizes of the task (cellset, angleset, groupset)
    int num_cellsets;
    int num_anglesets;
    int num_groupsets;
    
    int num_tasks;
    
    vector<Task*> all_tasks;
};


#endif
