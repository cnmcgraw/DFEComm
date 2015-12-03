//#include <Kripke/Grid.h>
#include "ParallelComm.h"
#include "Problem.h"

ParallelComm::ParallelComm(){

}

ParallelComm::~ParallelComm(){

}

void ParallelComm::setSizes(int cells, int angles, int groups){
  num_tasks = cells * angles * groups;
  
  // We also need to resize the various vectors
  recv_requests.resize(3*num_tasks, MPI_REQUEST_NULL);
  send_requests.resize(3*num_tasks, MPI_REQUEST_NULL);
  already_completed.resize(3*num_tasks, true);
  queue_depends.resize(num_tasks);
  all_tasks.resize(num_tasks, NULL); 
}

void ParallelComm::SetTask(Task* task){

    all_tasks[task->task_id] = task; 
}

int ParallelComm::computeTag(int mpi_rank, int task_id){
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int tag;
  tag = mpi_rank * num_tasks + task_id;

  return tag;
}

void ParallelComm::computeRankTask(int tag, int &mpi_rank, int &task_id){
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  mpi_rank = tag % mpi_size;
  task_id = tag / mpi_size;
}


/**
  Adds a task to the work queue.
  Determines if incoming dependencies require communication, and posts appropriate Irecv's.
*/
void ParallelComm::addTask(int task_id, Task &task){
  // Post recieves for incoming dependencies, and add to the queue
  postRecvs(task_id, task);
}

Task *ParallelComm::dequeueTask(int task_id){

  // Get task pointer before removing it from queue
 // Task *task = queue_tasks[task_id];
  Task *task = all_tasks[task_id];

  return task;
}

/**
  Adds a task to the work queue.
  Determines if incoming dependencies require communication, and posts appropirate Irecv's.
  All recieves use the plane_data[] arrays as recieve buffers.
*/
void ParallelComm::postRecvs(int task_id, Task &task){
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // go thru each dimensions incoming neighbors, and add the dependencies
  int num_depends = 0;
  for(int dim = 0;dim < 3;++ dim){
    // If it's a boundary condition, skip it
    if(task.incoming[dim][1] < 0){
      task.Get_buffer_from_bc(dim);
      continue;
    }
    // If it's an on-rank communication (from another Task)
    if(task.incoming[dim][1] == mpi_rank){
      // skip it, but track the dependency
      num_depends ++;
      continue;
    }
    // Change the flag for already completed to false
    already_completed[3*task_id + dim] = false;

    // compute the tag id of THIS task (tags are always based on destination)
     int tag = computeTag(mpi_rank, task_id);
    // Post the recieve
    MPI_Irecv(&task.plane_data[dim][0], task.plane_data[dim].size(), MPI_DOUBLE, task.incoming[dim][1],
      tag, MPI_COMM_WORLD, &recv_requests[3* task_id + dim]);


    // increment number of dependencies
    num_depends ++;
  }

  // add task to queue
  queue_depends[task_id] = num_depends;
}

void ParallelComm::postSends(Task *task, double *src_buffers[3]){
  // post sends for outgoing dependencies
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  for(int dim = 0;dim < 3;++ dim){
    // If it's a boundary condition, fill plane data with BC's
    if(task->outgoing[dim][1] < 0){
 //     task->Get_buffer_from_bc(dim);
      continue;
    }

    // If it's an on-rank communication (to another task)
    if(task->outgoing[dim][1] == mpi_rank){

      queue_depends[task->outgoing[dim][2]] --;

      // copy the boundary condition data into the outgoings plane data
      int task_index = task->outgoing[dim][2];
      Task& task_outgoing = *(all_tasks[task_index]);
      task_outgoing.plane_data[dim].assign(task->plane_data[dim].begin(), task->plane_data[dim].begin() + task->plane_data[dim].size());
      continue;
    }

    // At this point, we know that we have to send an MPI message
    // compute the tag id of TARGET task (tags are always based on destination)
    int tag = computeTag(task->outgoing[dim][1], task->outgoing[dim][2]);

    // Post the send
    MPI_Irsend(src_buffers[dim], task->plane_data[dim].size(), MPI_DOUBLE, task->outgoing[dim][1],
      tag, MPI_COMM_WORLD, &send_requests[task->task_id + dim]);
      
    int complete_flag = 0;
    MPI_Status status;
    MPI_Test(&send_requests[task->task_id + dim], &complete_flag, &status);
  }
}


/**
  Checks for incomming messages, and does relevant bookkeeping.
*/
void ParallelComm::testRecieves(int task_id, bool wait){

  int complete_flag = 0;
  std::vector<MPI_Status> recv_status(3);
  
   if(!wait){
     for(int i = 0; i < 3; i++){
       MPI_Test(&recv_requests[3*task_id + i], &complete_flag, &recv_status[i]);
       if(complete_flag && !already_completed[3*task_id + i]){
         queue_depends[task_id] -= 1;
         already_completed[3*task_id + i] = true;
       }
     }
   }
   else{
    MPI_Waitall(3, &recv_requests[3*task_id], &recv_status[0]);
    queue_depends[task_id] = 0;
   }


}

void ParallelComm::markComplete(int task_id){
  // Get task pointer and remove from work queue
 // Task *task = dequeueTask(task_id);
  Task *task = all_tasks[task_id];

  // Send new outgoing info for sweep
  double *buf[3] = {
    &(task->plane_data[0][0]),
    &(task->plane_data[1][0]),
    &(task->plane_data[2][0])
  };
  postSends(task, buf);
}

void ParallelComm::WaitforReady(int task_id){
  // Check to see if this task is ready
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  bool ready = false;
  if(queue_depends[task_id] == 0)
    ready = true;
  else{
    testRecieves(task_id, false);
    if(queue_depends[task_id] == 0)
      ready = true;
   }

  if(!ready){
    testRecieves(task_id, true);
  }
}

void ParallelComm::CleanUp(){
  queue_depends.clear();
  recv_requests.clear();
  send_requests.clear();
  already_completed.clear();
}
