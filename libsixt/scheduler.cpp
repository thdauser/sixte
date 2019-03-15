
#include "scheduler.h"

#include "log.h"

#include "threadsafe_queue.h"
#include "tasksSIRENA.h"

std::mutex end_workers_mut;
std::mutex end_eworkers_mut;
std::mutex records_detected_mut;
std::mutex records_energy_mut;

threadsafe_queue<sirena_data*> detection_queue;
threadsafe_queue<sirena_data*> detected_queue;
threadsafe_queue<sirena_data*> energy_queue;
threadsafe_queue<sirena_data*> end_queue;

scheduler* scheduler::instance = 0;

bool end_workers = false;
bool end_eworkers = false;

unsigned int records_detected = 0;
unsigned int records_energy = 0;

/* ****************************************************************************/
/* Workers ********************************************************************/
/* ****************************************************************************/
void detection_worker()
{
  //log_trace("Starting detection worker...");
  while(1){
    sirena_data* data;
    if(detection_queue.wait_and_pop(data)){
      //log_trace("Extracting detection data from queue...");
      th_runDetect(data->rec,
                data->last_record,
                data->all_pulses,
                &(data->rec_init),
                &(data->record_pulses));
      detected_queue.push(data);
      std::unique_lock<std::mutex> lk(records_detected_mut);
      ++records_detected;
      lk.unlock();
    }
    std::unique_lock<std::mutex> lk(end_workers_mut);
    if(end_workers){
      lk.unlock();
      break;
    }
    lk.unlock();
  }
}

void energy_worker()
{
  //log_trace("Starting energy worker...");
  while(1){
    sirena_data* data;
    if(energy_queue.wait_and_pop(data)){
      //log_trace("Extracting energy data from queue...");
      //log_debug("Energy data in record %i",data->n_record);
      th_runEnergy(data->rec, 
                   &(data->rec_init),
                   &(data->record_pulses),
                   &(data->optimal_filter));
      end_queue.push(data);
      std::unique_lock<std::mutex> lk(records_energy_mut);
      ++records_energy;
      lk.unlock();
    }
    std::unique_lock<std::mutex> lk(end_workers_mut);
    if(end_workers){
      lk.unlock();
      break;
    }
    lk.unlock();
  }
}

void energy_worker_v2()
{
  //log_trace("Starting energy worker...");
  while(1){
    sirena_data* data;
    if(detected_queue.wait_and_pop(data)){
      //log_trace("Extracting energy data from queue...");
      //log_debug("Energy data in record %i",data->n_record);
      th_runEnergy(data->rec, 
                   &(data->rec_init),
                   &(data->record_pulses),
                   &(data->optimal_filter));
      end_queue.push(data);
      std::unique_lock<std::mutex> lk(records_energy_mut);
      ++records_energy;
      lk.unlock();
    }
    std::unique_lock<std::mutex> lk(end_eworkers_mut);
    if(end_eworkers){
      lk.unlock();
      break;
    }
    lk.unlock();
  }
}

/* ****************************************************************************/
/* Scheduler ******************************************************************/
/* ****************************************************************************/
void scheduler::push_detection(TesRecord* record, 
                               int nRecord, 
                               int lastRecord, 
                               PulsesCollection *pulsesAll, 
                               ReconstructInitSIRENA** reconstruct_init, 
                               PulsesCollection** pulsesInRecord,
                               OptimalFilterSIRENA** optimal,
                               TesEventList* event_list)
{
  //log_trace("pushing detection data into the queue...");
  sirena_data* input = new sirena_data;
  tesrecord* rec = new tesrecord(record);
  input->rec = rec->get_TesRecord();
  input->n_record = nRecord;
  input->last_record = lastRecord;
  input->all_pulses = new PulsesCollection;
  if (pulsesAll and pulsesAll->ndetpulses > 0){
    *input->all_pulses = *pulsesAll;
  }

  input->record_pulses = *pulsesInRecord;
  input->rec_init = *reconstruct_init;
  input->optimal_filter = new OptimalFilterSIRENA;
  input->event_list = new TesEventList;
  input->event_list->size = event_list->size;

  input->event_list->index = event_list->index;
  input->event_list->size_energy = event_list->size_energy;
  input->event_list->event_indexes = new double[event_list->size];
  input->event_list->energies = new double[event_list->size];
  input->event_list->avgs_4samplesDerivative = new double[event_list->size];
  input->event_list->Es_lowres = new double[event_list->size];
  input->event_list->grades1 = new int[event_list->size];
  input->event_list->grades2 = new int[event_list->size];
  input->event_list->pulse_heights = new double[event_list->size];
  input->event_list->ph_ids = new long[event_list->size];
  input->event_list->grading = new int[event_list->size];
  input->event_list->phis = new double[event_list->size];
  input->event_list->lagsShifts = new int[event_list->size];
  detection_queue.push(input);
  ++num_records;
}

void scheduler::finish_reconstruction(ReconstructInitSIRENA* reconstruct_init,
                                      PulsesCollection** pulsesAll, 
                                      OptimalFilterSIRENA** optimalFilter)
{
  // Waits until all the records are detected
  // this works because this function should only be called
  // after all the records are queue
  //log_trace("Waiting until all the detection workers end");
  while(1){
    std::unique_lock<std::mutex> lk(records_detected_mut);
    if (records_detected == this->num_records){
      lk.unlock();
      break;
    }
    lk.unlock();
  }
  
  // Sorting the arrays by record number
  //log_trace("Sorting the arrays by record number");
  if ((*pulsesAll)->pulses_detected){
    delete [] (*pulsesAll)->pulses_detected;
  }
  
  (*pulsesAll)->size = this->num_records;
  (*pulsesAll)->ndetpulses = 0;
  (*pulsesAll)->pulses_detected = new PulseDetected[this->num_records];

  //log_debug("Number of records %i", this->num_records);
  this->data_array = new sirena_data*[this->num_records];//+1];
  while(!detected_queue.empty()){
    sirena_data* data;
    if(detected_queue.wait_and_pop(data)){
      data_array[data->n_record-1] = data;
    }
  }

  std::unique_lock<std::mutex> lk(end_workers_mut);
  end_workers = true;
  lk.unlock();
  for(unsigned int i = 0; i < this->max_detection_workers; ++i){
    this->detection_workers[i].join();
  }

  //
  // Energy
  //
  //log_trace("Starting energy workers...");
  if(this->is_running_energy){
    std::unique_lock<std::mutex> lk_end(end_workers_mut);
    end_workers = false;
    lk_end.unlock();
    this->energy_workers = new std::thread[this->max_detection_workers];
    for (unsigned int i = 0; i < this->max_detection_workers; ++i){//
      this->energy_workers[i] = std::thread (energy_worker);
    }
    
    //log_trace("Filling energy queue...");
    for (unsigned int i = 0; i < this->num_records; ++i){
      energy_queue.push(data_array[i]);
    }
  
    // Waits until all the energies are calculated
    //log_trace("Waiting until the energy workers end...");
    while(1){
      std::unique_lock<std::mutex> lk_energy(records_energy_mut);
      if (records_energy == this->num_records){
        lk_energy.unlock();
        break;
      }
      lk_energy.unlock();
    }
    
    std::unique_lock<std::mutex> lk_end2(end_workers_mut);
    end_workers = true;
    lk_end2.unlock();
    for(unsigned int i = 0; i < this->max_detection_workers; ++i){
      this->energy_workers[i].join();
    }
  }//end energy

  //
  // Reconstruction of the pulses array
  //
  //log_trace("Reconstruction of the pulses array...");
  for (unsigned int i = 0; i < this->num_records; ++i){
    PulsesCollection* all = data_array[i]->all_pulses;
    PulsesCollection* in_record = data_array[i]->record_pulses;
    PulsesCollection aux;
    if ((*pulsesAll)->size < 
        ((*pulsesAll)->ndetpulses + in_record->ndetpulses)){
      aux = *(*pulsesAll);
      delete *pulsesAll;
      *pulsesAll = new PulsesCollection;
      //**pulsesAll = aux;
      (*pulsesAll)->ndetpulses = aux.ndetpulses;
      (*pulsesAll)->size = (*pulsesAll)->ndetpulses + in_record->ndetpulses;
      (*pulsesAll)->pulses_detected =  new PulseDetected[(*pulsesAll)->size];
      for (unsigned int i = 0; i < aux.ndetpulses; ++i){
        (*pulsesAll)->pulses_detected[i] = aux.pulses_detected[i];
      }
    }
    for (unsigned int i = 0; i < in_record->ndetpulses; ++i){
      (*pulsesAll)->pulses_detected[i+(*pulsesAll)->ndetpulses] = 
        in_record->pulses_detected[i];
    }
    (*pulsesAll)->ndetpulses += in_record->ndetpulses;
  }// End reconstruction of the pulses array

  //log_trace("Filling eventlist...");
  for(unsigned int i = 0; i < this->num_records; ++i){
    
    TesEventList* event_list = data_array[i]->event_list;
    //PulsesCollection* pulses = data_array[i]->pulses_detected;
    PulsesCollection* record_pulses = data_array[i]->record_pulses;
    TesRecord* rec = data_array[i]->rec;
    /*if (event_list->energies != 0) delete [] event_list->energies;
    if (event_list->avgs_4samplesDerivative != 0) delete [] event_list->avgs_4samplesDerivative;
    if (event_list->grades1 != 0) delete [] event_list->grades1;
    if (event_list->grades2 != 0) delete [] event_list->grades2;
    if (event_list->pulse_heights != 0) delete [] event_list->pulse_heights;
    if (event_list->ph_ids != 0) delete [] event_list->ph_ids;*/

    if (strcmp(data_array[i]->rec_init->EnergyMethod,"PCA") != 0){
      event_list->index = record_pulses->ndetpulses;
      /*event_list->energies = new double[event_list->index];
      event_list->avgs_4samplesDerivative = new double[event_list->index];
      event_list->grades1  = new int[event_list->index];
      event_list->grades2  = new int[event_list->index];
      event_list->pulse_heights  = new double[event_list->index];
      event_list->ph_ids   = new long[event_list->index];*/

      for (int ip=0; ip < record_pulses->ndetpulses; ip++) {
        //log_debug("eventlist index %i record %i", ip, i);
        //log_debug("event_list indexes %p", event_list->event_indexes[ip]);
        //continue;
        event_list->event_indexes[ip] = 
          (record_pulses->pulses_detected[ip].Tstart - rec->time)/rec->delta_t;
        
        event_list->energies[ip] = data_array[i]->record_pulses->pulses_detected[ip].energy;
        
        event_list->avgs_4samplesDerivative[ip] = 
          record_pulses->pulses_detected[ip].avg_4samplesDerivative;
        event_list->Es_lowres[ip] = record_pulses->pulses_detected[ip].E_lowres;
        event_list->grading[ip] = record_pulses->pulses_detected[ip].grading;
        event_list->grades1[ip]  = record_pulses->pulses_detected[ip].grade1;
        event_list->grades2[ip]  = record_pulses->pulses_detected[ip].grade2;
        event_list->pulse_heights[ip]  = record_pulses->pulses_detected[ip].pulse_height;
        event_list->ph_ids[ip]   = 0;
        event_list->phis[ip] = record_pulses->pulses_detected[ip].phi;
        event_list->lagsShifts[ip] = record_pulses->pulses_detected[ip].lagsShift;
      }
      if (data_array[i]->last_record == 1) {
        //log_debug("eventlist last record");
        double numLagsUsed_mean;
        double numLagsUsed_sigma;
        gsl_vector *numLagsUsed_vector = gsl_vector_alloc((*pulsesAll)->ndetpulses);
        
        for (int ip = 0; ip < (*pulsesAll)->ndetpulses; ip++) {
          gsl_vector_set(numLagsUsed_vector,ip,(*pulsesAll)->pulses_detected[ip].numLagsUsed);
        }
        if (findMeanSigma (numLagsUsed_vector, &numLagsUsed_mean, &numLagsUsed_sigma)) {
          EP_EXIT_ERROR("Cannot run findMeanSigma routine for calculating numLagsUsed statistics",EPFAIL);
        }
        gsl_vector_free(numLagsUsed_vector);
      }
    }else{
      if (data_array[i]->last_record == 1) {
        // Free & Fill TesEventList structure
        /*event_list->index = (*pulsesAll)->ndetpulses;
        event_list->event_indexes = new double[event_list->index];
        event_list->energies = new double[event_list->index];
        event_list->avgs_4samplesDerivative = new double[event_list->index];
        event_list->grades1  = new int[event_list->index];
        event_list->grades2  = new int[event_list->index];
        event_list->pulse_heights  = new double[event_list->index];
        event_list->ph_ids   = new long[event_list->index];*/
        
        for (int ip = 0; ip<(*pulsesAll)->ndetpulses; ip++) {
          //log_debug("eventlist index %i record %i", ip, i);
          event_list->event_indexes[ip] = 
            ((*pulsesAll)->pulses_detected[ip].Tstart - rec->time)/rec->delta_t;
          
          event_list->energies[ip] = (*pulsesAll)->pulses_detected[ip].energy;
          
          event_list->avgs_4samplesDerivative[ip]  = (*pulsesAll)->pulses_detected[ip].avg_4samplesDerivative;
          event_list->Es_lowres[ip]  = (*pulsesAll)->pulses_detected[ip].E_lowres;
          event_list->grading[ip] = (*pulsesAll)->pulses_detected[ip].grading;
          event_list->grades1[ip]  = (*pulsesAll)->pulses_detected[ip].grade1;
          event_list->grades2[ip]  = (*pulsesAll)->pulses_detected[ip].grade2;
          event_list->pulse_heights[ip]  = (*pulsesAll)->pulses_detected[ip].pulse_height;
          event_list->ph_ids[ip]   = 0;    
          event_list->phis[ip] = record_pulses->pulses_detected[ip].phi;
          event_list->lagsShifts[ip] = record_pulses->pulses_detected[ip].lagsShift;
        }
      }
    }

    //log_debug("Eventlist from record %i", (i + 1) );
    //log_debug("%i, %i, %i",event_list->size, event_list->size_energy, event_list->index);
    for (int j = 0; j < event_list->index; ++j){
      /*log_debug("%f, %f, %f, %f, %i %i, %ld", event_list->event_indexes[j],
                event_list->pulse_heights[j],
                event_list->avgs_4samplesDerivative[j],
                event_list->energies[j],
                event_list->grades1[j],
                event_list->grades2[j],
                event_list->ph_ids[j]);*/
    }
  }// for event_list
}

void scheduler::finish_reconstruction_v2(ReconstructInitSIRENA* reconstruct_init,
                                         PulsesCollection** pulsesAll, 
                                         OptimalFilterSIRENA** optimalFilter)
{
  // Waits until all the records are detected
  // this works because this function should only be called
  // after all the records are queue
  //log_trace("Waiting until all the detection workers end");
  while(1){
    std::unique_lock<std::mutex> lk(records_detected_mut);
    if (records_detected == this->num_records){
      lk.unlock();
      break;
    }
    lk.unlock();
  }

  std::unique_lock<std::mutex> lk(end_workers_mut);
  end_workers = true;
  lk.unlock();
  for(unsigned int i = 0; i < this->max_detection_workers; ++i){
    this->detection_workers[i].join();
  }

  // Waits until all the energies are calculated
  //log_trace("Waiting until the energy workers end...");
  while(1){
    std::unique_lock<std::mutex> lk_energy(records_energy_mut);
    if (records_energy == this->num_records){
      lk_energy.unlock();
      break;
    }
    lk_energy.unlock();
  }
  
  std::unique_lock<std::mutex> lk_end2(end_eworkers_mut);
  end_eworkers = true;
  lk_end2.unlock();
  for(unsigned int i = 0; i < this->max_energy_workers; ++i){
    this->energy_workers[i].join();
  }

  // Sorting the arrays by record number
  //log_trace("Sorting the arrays by record number");
  if ((*pulsesAll)->pulses_detected){
    delete [] (*pulsesAll)->pulses_detected;
  }
  
  (*pulsesAll)->size = this->num_records;
  (*pulsesAll)->ndetpulses = 0;
  (*pulsesAll)->pulses_detected = new PulseDetected[this->num_records];

  //log_debug("Number of records %i", this->num_records);
  this->data_array = new sirena_data*[this->num_records];//+1];
  while(!end_queue.empty()){
    sirena_data* data;
    if(end_queue.wait_and_pop(data)){
      data_array[data->n_record-1] = data;
    }
  }
  //
  // Reconstruction of the pulses array
  //
  //log_trace("Reconstruction of the pulses array...");
  for (unsigned int i = 0; i < this->num_records; ++i){
    PulsesCollection* all = data_array[i]->all_pulses;
    PulsesCollection* in_record = data_array[i]->record_pulses;
    PulsesCollection aux;
    if ((*pulsesAll)->size < 
        ((*pulsesAll)->ndetpulses + in_record->ndetpulses)){
      aux = *(*pulsesAll);
      delete *pulsesAll;
      *pulsesAll = new PulsesCollection;
      (*pulsesAll)->ndetpulses = aux.ndetpulses;
      (*pulsesAll)->size = (*pulsesAll)->ndetpulses + in_record->ndetpulses;
      (*pulsesAll)->pulses_detected =  new PulseDetected[(*pulsesAll)->size];
      for (unsigned int i = 0; i < aux.ndetpulses; ++i){
        (*pulsesAll)->pulses_detected[i] = aux.pulses_detected[i];
      }
    }
    for (unsigned int i = 0; i < in_record->ndetpulses; ++i){
      (*pulsesAll)->pulses_detected[i+(*pulsesAll)->ndetpulses] = 
        in_record->pulses_detected[i];
    }
    (*pulsesAll)->ndetpulses += in_record->ndetpulses;
  }// End reconstruction of the pulses array

  //log_debug("pulsesAll");
  for (int i = 0; i < (*pulsesAll)->ndetpulses; ++i){
    /*log_debug("data - %i, %i, %i, %i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %i",
              (*pulsesAll)->pulses_detected[i].pulse_duration,
              (*pulsesAll)->pulses_detected[i].grade1,
              (*pulsesAll)->pulses_detected[i].grade2,
              (*pulsesAll)->pulses_detected[i].pixid,
              (*pulsesAll)->pulses_detected[i].Tstart,
              (*pulsesAll)->pulses_detected[i].Tend,
              (*pulsesAll)->pulses_detected[i].riseTime,
              (*pulsesAll)->pulses_detected[i].fallTime,
              (*pulsesAll)->pulses_detected[i].pulse_height,
              (*pulsesAll)->pulses_detected[i].maxDER,
              (*pulsesAll)->pulses_detected[i].samp1DER,
              (*pulsesAll)->pulses_detected[i].avg_4samplesDerivative,
              (*pulsesAll)->pulses_detected[i].quality,
              (*pulsesAll)->pulses_detected[i].numLagsUsed);*/
    for (int j = 0; j < (*pulsesAll)->pulses_detected[i].pulse_adc->size; ++j){
      //log_debug("%d, ",*((*pulsesAll)->pulses_detected[i].pulse_adc->data));
    }
  }

  //log_trace("Filling eventlist...");
  for(unsigned int i = 0; i < this->num_records; ++i){
    
    TesEventList* event_list = data_array[i]->event_list;
    //PulsesCollection* pulses = data_array[i]->pulses_detected;
    PulsesCollection* record_pulses = data_array[i]->record_pulses;
    TesRecord* rec = data_array[i]->rec;
    /*if (event_list->energies != 0) delete [] event_list->energies;
    if (event_list->avgs_4samplesDerivative != 0) delete [] event_list->avgs_4samplesDerivative;
    if (event_list->grades1 != 0) delete [] event_list->grades1;
    if (event_list->grades2 != 0) delete [] event_list->grades2;
    if (event_list->pulse_heights != 0) delete [] event_list->pulse_heights;
    if (event_list->ph_ids != 0) delete [] event_list->ph_ids;*/

    if (strcmp(data_array[i]->rec_init->EnergyMethod,"PCA") != 0){
      event_list->index = record_pulses->ndetpulses;
      /*event_list->energies = new double[event_list->index];
      event_list->avgs_4samplesDerivative = new double[event_list->index];
      event_list->grades1  = new int[event_list->index];
      event_list->grades2  = new int[event_list->index];
      event_list->pulse_heights  = new double[event_list->index];
      event_list->ph_ids   = new long[event_list->index];*/

      for (int ip=0; ip < record_pulses->ndetpulses; ip++) {
        //log_debug("eventlist index %i record %i", ip, i);
        //log_debug("event_list indexes %p", event_list->event_indexes[ip]);
        //continue;
        event_list->event_indexes[ip] = 
          (record_pulses->pulses_detected[ip].Tstart - rec->time)/rec->delta_t;
        
        event_list->energies[ip] = data_array[i]->record_pulses->pulses_detected[ip].energy;
        
        event_list->avgs_4samplesDerivative[ip] = 
          record_pulses->pulses_detected[ip].avg_4samplesDerivative;
        event_list->Es_lowres[ip] = record_pulses->pulses_detected[ip].E_lowres;
        event_list->grading[ip] = record_pulses->pulses_detected[ip].grading;
        event_list->grades1[ip]  = record_pulses->pulses_detected[ip].grade1;
        event_list->grades2[ip]  = record_pulses->pulses_detected[ip].grade2;
        event_list->pulse_heights[ip]  = record_pulses->pulses_detected[ip].pulse_height;
        event_list->ph_ids[ip]   = 0;
        event_list->phis[ip] = record_pulses->pulses_detected[ip].phi;
        event_list->lagsShifts[ip] = record_pulses->pulses_detected[ip].lagsShift;
      }
      if (data_array[i]->last_record == 1) {
        //log_debug("eventlist last record");
        double numLagsUsed_mean;
        double numLagsUsed_sigma;
        gsl_vector *numLagsUsed_vector = gsl_vector_alloc((*pulsesAll)->ndetpulses);
        
        for (int ip = 0; ip < (*pulsesAll)->ndetpulses; ip++) {
          gsl_vector_set(numLagsUsed_vector,ip,(*pulsesAll)->pulses_detected[ip].numLagsUsed);
        }
        if (findMeanSigma (numLagsUsed_vector, &numLagsUsed_mean, &numLagsUsed_sigma)) {
          EP_EXIT_ERROR("Cannot run findMeanSigma routine for calculating numLagsUsed statistics",EPFAIL);
        }
        gsl_vector_free(numLagsUsed_vector);
      }
    }else{
      if (data_array[i]->last_record == 1) {
        // Free & Fill TesEventList structure
        /*event_list->index = (*pulsesAll)->ndetpulses;
        event_list->event_indexes = new double[event_list->index];
        event_list->energies = new double[event_list->index];
        event_list->avgs_4samplesDerivative = new double[event_list->index];
        event_list->grades1  = new int[event_list->index];
        event_list->grades2  = new int[event_list->index];
        event_list->pulse_heights  = new double[event_list->index];
        event_list->ph_ids   = new long[event_list->index];*/
        
        for (int ip = 0; ip<(*pulsesAll)->ndetpulses; ip++) {
          //log_debug("eventlist index %i record %i", ip, i);
          event_list->event_indexes[ip] = 
            ((*pulsesAll)->pulses_detected[ip].Tstart - rec->time)/rec->delta_t;
          
          event_list->energies[ip] = (*pulsesAll)->pulses_detected[ip].energy;
          
          event_list->avgs_4samplesDerivative[ip]  = (*pulsesAll)->pulses_detected[ip].avg_4samplesDerivative;
          event_list->Es_lowres[ip]  = (*pulsesAll)->pulses_detected[ip].E_lowres;
          event_list->grading[ip] = (*pulsesAll)->pulses_detected[ip].grading;
          event_list->grades1[ip]  = (*pulsesAll)->pulses_detected[ip].grade1;
          event_list->grades2[ip]  = (*pulsesAll)->pulses_detected[ip].grade2;
          event_list->pulse_heights[ip]  = (*pulsesAll)->pulses_detected[ip].pulse_height;
          event_list->ph_ids[ip]   = 0;    
          event_list->phis[ip] = record_pulses->pulses_detected[ip].phi;
          event_list->lagsShifts[ip] = record_pulses->pulses_detected[ip].lagsShift;
        }
      }
    }

    //log_debug("Eventlist from record %i", (i + 1) );
    //log_debug("%i, %i, %i",event_list->size, event_list->size_energy, event_list->index);
    for (int j = 0; j < event_list->index; ++j){
      /*log_debug("%f, %f, %f, %f, %i %i, %ld", event_list->event_indexes[j],
                event_list->pulse_heights[j],
                event_list->avgs_4samplesDerivative[j],
                event_list->energies[j],
                event_list->grades1[j],
                event_list->grades2[j],
                event_list->ph_ids[j]);*/
    }
  }// for event_list
}

void scheduler::get_test_event(TesEventList** test_event, TesRecord** record)
{
  if(this->current_record == this->num_records) return;
  //log_trace("Getting eventlist from record %i", (this->current_record + 1));
  *test_event = this->data_array[this->current_record]->event_list;
  *record = this->data_array[this->current_record]->rec;
  this->current_record++;
}

void scheduler::init()
{
  if(threading){
    this->num_cores = std::thread::hardware_concurrency();
    //this->num_cores = 1;
    if(this->num_cores < 2){
      this->max_detection_workers = 1;
    }else{
      this->max_detection_workers = this->num_cores - 1;
    }
    this->detection_workers = new std::thread[this->max_detection_workers];
    for (unsigned int i = 0; i < this->max_detection_workers; ++i){
      this->detection_workers[i] = std::thread (detection_worker);
    }
  }
}

void scheduler::init_v2()
{
  if(threading){
    this->num_cores = std::thread::hardware_concurrency();
    if(this->num_cores < 2){
      this->max_detection_workers = 1;
      this->max_energy_workers = 1;
    }else{
      this->max_energy_workers = (this->num_cores - 1) / 2;
      this->max_detection_workers = 
        (this->num_cores - 1) - this->max_energy_workers;
      /*log_debug("detection %u energy %u", this->max_detection_workers,
                this->max_energy_workers);*/
    }
    this->detection_workers = new std::thread[this->max_detection_workers];
    for (unsigned int i = 0; i < this->max_detection_workers; ++i){
      this->detection_workers[i] = std::thread (detection_worker);
    }
    this->energy_workers = new std::thread[this->max_energy_workers];
    for (unsigned int i = 0; i < this->max_energy_workers; ++i){
      this->energy_workers[i] = std::thread (energy_worker_v2);
    }
  }
}

scheduler::scheduler():
  threading(false),      // true: Activate THREADING, false: No THREADING
  num_cores(0),
  max_detection_workers(0),
  max_energy_workers(0),
  num_records(0),
  is_running_energy(false),
  current_record(0),
  data_array(0)
{
  this->init_v2();
}

scheduler::~scheduler()
{
  if(threading){
    instance = 0;
  }
}

/* ****************************************************************************/
/* Data structures implementation *********************************************/
/* ****************************************************************************/

phidlist::phidlist():
  phid_array(0),
  times(0),
  wait_list(0),
  n_elements(0),
  index(0),
  size(0)
{
  
}

phidlist::phidlist(const phidlist& other):
  phid_array(0),
  times(0),
  wait_list(other.wait_list),
  n_elements(other.n_elements),
  index(other.index),
  size(other.size)
{
  if (other.phid_array && other.size > 0){
    phid_array = new long[other.size];
    for (unsigned int i = 0; i < size; ++i){
      phid_array[i] = other.phid_array[i];
    }
  }
  if (other.wait_list && other.times && other.size > 0){
    times = new double[other.size];
    for (unsigned int i = 0; i < size; ++i){
      times[i] = other.times[i];
    }
  }
}

phidlist::phidlist(PhIDList* other):
  phid_array(0),
  times(0),
  wait_list(other->wait_list),
  n_elements(other->n_elements),
  index(other->index),
  size(other->size)
{
  if (other->phid_array && size > 0){
    phid_array = new long[size];
    for (unsigned int i = 0; i < size; ++i){
      phid_array[i] = other->phid_array[i];
    }
  }
  
  if (wait_list && other->times && size > 0){
    times = new double[size];
    for (unsigned int i = 0; i < size; ++i){
      times[i] = other->times[i];
    }
  }
}

phidlist& phidlist::operator=(const phidlist& other)
{
  if(this != &other){
    wait_list = other.wait_list;
    n_elements = other.n_elements;
    index = other.index;
    size = other.size;
    if(phid_array){
      delete [] phid_array; phid_array = 0;
    }
    if (times){
      delete [] times; times = 0;
    }
    if (other.phid_array && other.size > 0){
      phid_array = new long[other.size];
      for (unsigned int i = 0; i < size; ++i){
        phid_array[i] = other.phid_array[i];
      }
    }
    if (other.wait_list && other.times && other.size > 0){
      times = new double[other.size];
      for (unsigned int i = 0; i < size; ++i){
        times[i] = other.times[i];
      }
    }
  }
  return *this;
}

phidlist::~phidlist()
{
  if (phid_array){
    delete [] phid_array; phid_array = 0;
  }
  if (times){
    delete [] times; times = 0;
  }
}

PhIDList* phidlist::get_PhIDList() const
{
  PhIDList* ret = new PhIDList;
  ret->wait_list = wait_list;
  ret->n_elements = n_elements;
  ret->index = index;
  ret->size = size;
  if (phid_array && size > 0){
    ret->phid_array = new long[size];
    for (unsigned int i = 0; i < size; ++i){
      ret->phid_array[i] = phid_array[i];
    }
  }
  if (wait_list && times && size > 0){
    ret->times = new double[size];
    for (unsigned int i = 0; i < size; ++i){
      ret->times[i] = times[i];
    }
  }
  return ret;
}

tesrecord::tesrecord():
  trigger_size(0),
  time(0),
  delta_t(0),
  adc_array(0),
  adc_double(0),
  pixid(0),
  phid_list(0)
{
  
}

tesrecord::tesrecord(TesRecord* other):
  trigger_size(other->trigger_size),
  time(other->time),
  delta_t(other->delta_t),
  adc_array(0),//trigger_size
  adc_double(0),//trigger_size
  pixid(other->pixid),
  phid_list(0)//MAXIMPACTNUMBER
{
  if(other->adc_array && trigger_size > 0){
    adc_array = new uint16_t[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      adc_array[i] = other->adc_array[i];
    }
  }
  if(other->adc_double && trigger_size > 0){
    adc_double = new double[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      adc_double[i] = other->adc_double[i];
    }
  }
  if(other->phid_list && other->phid_list->size > 0){
    phid_list = new phidlist(other->phid_list);
  }
}

tesrecord::tesrecord(const tesrecord& other):
  trigger_size(other.trigger_size),
  time(other.time),
  delta_t(other.delta_t),
  adc_array(0),
  adc_double(0),
  pixid(other.pixid),
  phid_list(0)
{
  if(other.adc_array && trigger_size > 0){
    adc_array = new uint16_t[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      adc_array[i] = other.adc_array[i];
    }
  }
  if(other.adc_double && trigger_size > 0){
    adc_double = new double[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      adc_double[i] = other.adc_double[i];
    }
  }
  if(other.phid_list && other.phid_list->size > 0){
    phid_list = new phidlist(*other.phid_list);
  }
}

tesrecord& tesrecord::operator=(const tesrecord& other)
{
  if(this != &other){
    trigger_size = other.trigger_size;
    time = other.time;
    delta_t = other.delta_t;
    pixid = other.pixid;
    
    if(adc_array){
      delete [] adc_array; adc_array = 0;
    }

    if(adc_double){
      delete [] adc_double; adc_double = 0;
    }

    if(phid_list){
      delete phid_list; phid_list = 0;
    }

    if(other.adc_array && trigger_size > 0){
      adc_array = new uint16_t[trigger_size];
      for (unsigned int i = 0; i < trigger_size; ++i){
        adc_array[i] = other.adc_array[i];
      }
    }
    if(other.adc_double && trigger_size > 0){
      adc_double = new double[trigger_size];
      for (unsigned int i = 0; i < trigger_size; ++i){
        adc_double[i] = other.adc_double[i];
      }
    }
    if(other.phid_list && other.phid_list->size > 0){
      phid_list = new phidlist(*other.phid_list);
    }
  }
  return *this;
}

tesrecord::~tesrecord()
{
  if (adc_array){
    delete [] adc_array; adc_array = 0;
  }
  if (adc_double){
    delete [] adc_double; adc_double = 0;
  }
  if (phid_list){
    delete phid_list; phid_list = 0;
  }
}

TesRecord* tesrecord::get_TesRecord() const
{
  TesRecord* ret = new TesRecord;
  ret->trigger_size = trigger_size;
  ret->time = time;
  ret->delta_t = delta_t;
  ret->pixid = pixid;

  if(adc_array && trigger_size > 0){
    ret->adc_array = new uint16_t[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      ret->adc_array[i] = adc_array[i];
    }
  }

  if(adc_double && trigger_size > 0){
    ret->adc_double = new double[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
        ret->adc_double[i] = adc_double[i];
    }
  }
  if(phid_list && phid_list->size > 0){
    ret->phid_list = phid_list->get_PhIDList();
  }
  return ret;
}

data::data():
  n_record(0),
  last_record(0),
  all_pulses(0),
  record_pulses(0)
{
  
}

data::data(const data& other):
  n_record(other.n_record),
  last_record(other.last_record),
  all_pulses(0),
  record_pulses(0),
  rec(other.rec),
  rec_init(other.rec_init)
{
  
}

data& data::operator=(const data& other)
{
  printf("operator = date\n");
  if(this != &other){
    n_record = other.n_record;
    last_record = other.last_record;
    rec = other.rec;
    rec_init = other.rec_init;
  }
  return *this;
}

data::~data()
{

}


