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

int ParallelComm::findTaskinTasks(int task_id){
  // find task in all_tasks
  int index;
  for(index = 0;index < all_tasks.size();++ index){
    if(all_tasks[index]->task_id == task_id){
      break;
    }
  }
  return index;

}

int ParallelComm::computeTag(int mpi_rank, int task_id){
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int tag;
/* tag = mpi_rank * (num_cellsets * num_groupsets * num_anglesets)
             + cs * (num_groupsets * num_anglesets) 
             + gs * num_anglesets + as; */
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
  Finds task in the queue by its task id.
*/
int ParallelComm::findTask(int task_id){

  // find task in queue
  int index;
  for(index = 0;index < queue_task_ids.size();++ index){
    if(queue_task_ids[index] == task_id){
      break;
    }
  }
  if(index == queue_task_ids.size()){
    printf("Cannot find task id %d in work queue\n", task_id);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  return index;
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
  int index = findTask(task_id);

  // Get task pointer before removing it from queue
  Task *task = queue_tasks[index];

  // remove task from queue
  queue_task_ids.erase(queue_task_ids.begin()+index);
  queue_tasks.erase(queue_tasks.begin()+index);
  queue_depends.erase(queue_depends.begin()+index);

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

    // compute the tag id of THIS task (tags are always based on destination)
   // int tag = computeTag(mpi_rank, task.cellset_id, task.angleset_id, task.groupset_id);
     int tag = computeTag(mpi_rank, task_id);

    // Post the recieve
    MPI_Irecv(&task.plane_data[dim][0], task.plane_data[dim].size(), MPI_DOUBLE, task.incoming[dim][1],
      tag, MPI_COMM_WORLD, &recv_requests[recv_requests.size()-1]);


    // increment number of dependencies
    num_depends ++;
  }

  // add task to queue
  queue_task_ids.push_back(task_id);
  queue_tasks.push_back(&task);
  queue_depends.push_back(num_depends);
}

void ParallelComm::postSends(Task *task, double *src_buffers[3]){
  // post sends for outgoing dependencies
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  for(int dim = 0;dim < 3;++ dim){
    // If it's a boundary condition, skip it
    if(task->outgoing[dim][1] < 0)
      continue;

    // If it's an on-rank communication (to another task)
    if(task->outgoing[dim][1] == mpi_rank){
      // find the local task in the queue, and decrement the counter
      for(int i = 0;i < queue_task_ids.size();++ i){
        if(queue_task_ids[i] == task->outgoing[dim][2]){
          queue_depends[i] --;
          break;
        }
      }

      // copy the boundary condition data into the outgoings plane data
      int task_index = findTaskinTasks(task->outgoing[dim][2]);
      Task& task_outgoing = *(all_tasks[task_index]);
      task_outgoing.plane_data[dim].assign(task->plane_data[dim].begin(), task->plane_data[dim].begin() + task->plane_data[dim].size());
      continue;
    }

    // At this point, we know that we have to send an MPI message
    // Add request to send queue
    send_requests.push_back(MPI_Request());

    // compute the tag id of TARGET task (tags are always based on destination)
   // int tag = computeTag(task->outgoing[dim][1], task->outgoing[dim][2], task->angleset_id, task->groupset_id);
    int tag = computeTag(task->outgoing[dim][1], task->outgoing[dim][2]);

    // Post the send
    MPI_Irsend(src_buffers[dim], task->plane_data[dim].size(), MPI_DOUBLE, task->outgoing[dim][1],
      tag, MPI_COMM_WORLD, &send_requests[send_requests.size()-1]);
  }
}


// Checks if there are any outstanding task to complete
bool ParallelComm::workRemaining(void){
  // If there are outstanding task to process, return true
  if (recv_requests.size() > 0 || queue_tasks.size() > 0){
    return true;
  }
  
  // No more work, so make sure all of our sends have completed
  // before we continue
  waitAllSends(); 
  
}

/**
  Checks for incomming messages, and returns a list of ready task id's
*/
std::vector<int> ParallelComm::readyTasks(void){
  // check for incomming messages
  testRecieves();

  // build up a list of ready tasks
  return getReadyList();
}


// Blocks until all sends have completed, and flushes the send queues
void ParallelComm::waitAllSends(void){
  // Wait for all remaining sends to complete, then return false
  int num_sends = send_requests.size();
  if(num_sends > 0){
    std::vector<MPI_Status> status(num_sends);
    MPI_Waitall(num_sends, &send_requests[0], &status[0]);
    send_requests.clear();
  }
}

/**
  Checks for incomming messages, and does relevant bookkeeping.
*/
void ParallelComm::testRecieves(void){

  // Check for any recv requests that have completed
  int num_requests = recv_requests.size();
  bool done = false;
  while(!done && num_requests > 0){
    // Create array of status variables
    std::vector<MPI_Status> recv_status(num_requests);

    // Ask if either one or none of the recvs have completed?
    int index; // this will be the index of request that completed
    int complete_flag; // this is set to TRUE if something completed
    MPI_Testany(num_requests, &recv_requests[0], &index, &complete_flag, &recv_status[0]);

    if(complete_flag != 0){

      // get task that this completed for
      int task_id = recv_tasks[index];

      // remove the request from the list
      recv_requests.erase(recv_requests.begin()+index);
      recv_tasks.erase(recv_tasks.begin()+index);
      num_requests --;

      // decrement the dependency count for that task
      for(int i = 0;i < queue_task_ids.size();++ i){
        if(queue_task_ids[i] == task_id){
          queue_depends[i] --;
          break;
        }
      }
    }
    else{
      done = true;
    }
  }
}


std::vector<int> ParallelComm::getReadyList(void){
  // build up a list of ready tasks
  std::vector<int> ready;
  for(int i = 0;i < queue_depends.size();++ i){
    if(queue_depends[i] == 0){
      ready.push_back(queue_task_ids[i]);
    }
  }
  return ready;
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
 
  testRecieves();
  
  // Check to see if this task is ready
  int index = findTask(task_id);
  bool ready = false;
  if(queue_depends[index] == 0)
    ready = true;

  // If it's not, Wait until it is
  while(!ready){
    testRecieves();

    if(queue_depends[index] == 0){

      ready = true; }
  }
}
