//#include <Kripke/Grid.h>
#include "ParallelComm.h"
#include "Problem.h"

ParallelComm::ParallelComm(){

}

ParallelComm::~ParallelComm(){

}

void ParallelComm::setSizes(int cells, int angles, int groups){
  num_tasks = cells * angles * groups;
}

void ParallelComm::SetTask(Task* Task){

    all_tasks.push_back(Task); 
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
  Determines if incoming dependencies require communication, and posts appropirate Irecv's.
*/
void ParallelComm::addTask(int task_id, Task &task){
  // Post recieves for incoming dependencies, and add to the queue
  postRecvs(task_id, task);
}

Task *ParallelComm::dequeueTask(int task_id){

  // Get task pointer
  Task *task = queue_tasks[task_id];
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  return task;
}

/**
  Adds a task to the work queue.
  Determines if incoming dependencies require communication, and posts appropirate Irecv's.
  All recieves use the plane_data[] arrays as recieve buffers.
*/
void ParallelComm::postRecvs(int task_id, Task &task){
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // go thru each dimensions incoming neighbors, and add the dependencies
  int num_depends = 0;
  int num_request = 0;
  for(int dim = 0;dim < 3;++ dim){
    // If it's a boundary condition, skip it
    if(task.incoming[dim][1] < 0){
      continue;
    }
    // If it's an on-rank communication (from another Task)
    if(task.incoming[dim][1] == mpi_rank){
      // skip it, but track the dependency
      num_depends ++;
      continue;
    }
    // Add request to pending list
    recv_requests.push_back(MPI_Request());
    recv_tasks.push_back(task_id);
    request_index.push_back(recv_requests.size()-1);
    already_completed.push_back(false);
    // compute the tag id of THIS task (tags are always based on destination)
     int tag = computeTag(mpi_rank, task_id);
    // Post the recieve
    MPI_Irecv(&task.plane_data[dim][0], task.plane_data[dim].size(), MPI_DOUBLE, task.incoming[dim][1],
      tag, MPI_COMM_WORLD, &recv_requests[recv_requests.size()-1]);


    // increment number of dependencies
    num_depends ++;
    num_request ++;
  }

  // add task to queue
  queue_task_ids.push_back(task_id);
  queue_tasks.push_back(&task);
  queue_depends.push_back(num_depends);
  queue_request.push_back(num_request);
  if(num_request == 0)
    queue_start_recv.push_back(-1);
  else
    queue_start_recv.push_back(request_index.back() - num_request + 1);
}

void ParallelComm::postSends(Task *task, double *src_buffers[3]){
  // post sends for outgoing dependencies
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  for(int dim = 0;dim < 3;++ dim){
    // If it's a boundary condition, fill plane data with BC's
    if(task->outgoing[dim][1] < 0){
      task->Get_buffer_from_bc(dim);
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
    // Add request to send queue
    send_requests.push_back(MPI_Request());

    // compute the tag id of TARGET task (tags are always based on destination)
    int tag = computeTag(task->outgoing[dim][1], task->outgoing[dim][2]);

    // Post the send
    MPI_Irsend(src_buffers[dim], task->plane_data[dim].size(), MPI_DOUBLE, task->outgoing[dim][1],
      tag, MPI_COMM_WORLD, &send_requests[send_requests.size()-1]);
  }
}

/**
  Checks for incomming messages, and does relevant bookkeeping.
*/
void ParallelComm::testRecieves(int task_id){

  int start_index = queue_start_recv[task_id];
  int num_request = queue_request[task_id];

  int complete_flag;
  MPI_Status recv_status;
  for(int i = 0; i < num_request; i++){
    MPI_Test(&recv_requests[start_index + i], &complete_flag, &recv_status);
    if(complete_flag && !already_completed[start_index + i]){
      queue_depends[task_id] -= 1;
      already_completed[start_index + i] = true;
    }
  }

   
  // // Check for any recv requests that have completed
  // int num_requests = recv_requests.size();
  // bool done = false;
  // while(!done && num_requests > 0){
    // // Create array of status variables
    // std::vector<MPI_Status> recv_status(num_requests);

    // // Ask if either one or none of the recvs have completed?
    // int index; // this will be the index of request that completed
    // int complete_flag; // this is set to TRUE if something completed
    // MPI_Testany(num_requests, &recv_requests[0], &index, &complete_flag, &recv_status[0]);

    // if(complete_flag != 0){

      // // get task that this completed for
      // int task_id = recv_tasks[index];

      // // remove the request from the list
      // recv_requests.erase(recv_requests.begin()+index);
      // recv_tasks.erase(recv_tasks.begin()+index);
      // num_requests --;

      // // decrement the dependency count for that task
      // queue_depends[task_id] --;

    // }
    // else{
      // done = true;
    // }
  // }
}

void ParallelComm::markComplete(int task_id){
  // Get task pointer and remove from work queue
  Task *task = dequeueTask(task_id);

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
 // int index = findTask(task_id);
  bool ready = false;
  if(queue_depends[task_id] == 0)
    ready = true;
  else{
    testRecieves(task_id);
    if(queue_depends[task_id] == 0)
      ready = true;
   }

  // If it's not, Wait until it is
  while(!ready){
    testRecieves(task_id);

    if(queue_depends[task_id] == 0){
      ready = true; }
  }
}

void ParallelComm::CleanUp(){
  queue_task_ids.clear();
  queue_tasks.clear();
  queue_depends.clear();
  recv_requests.clear();
  recv_tasks.clear();
  request_index.clear();
  send_requests.clear();
}
