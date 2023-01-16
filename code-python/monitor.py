import matplotlib.pyplot as plt
import numpy as np
import psutil
import time

"""
Code to track CPU and memory usage.
Author Andrew Bowell.
Tidied up by Peter Thomas.
"""

class Monitor:

    """
    Routines to monitor CPU and memory resource.
    These routines should be rewritten to use named dtype arrays: it is horrible at the moment.
    Attributes
    ----------
    Methods
    -------
    """
    def __init__(self, n_graph, n_func):
        
        self.n_func = n_func

        # To make this a general routine, these descriptions should be passed in as an array or list
        self.function_descs ="""
        1. Initialise Halo properties class (does also find most central halo). \n
        2. Calculate mass in common. \n
        3. Gather DM mass from progenitors. \n
        4. Give hot gas etc to descendents. \n
        5. Calc and set baryon fraction. \n
        6. Collect galaxy progenitor info.This includes giving stellar mass to descendents. \n
        7. Simple function to add bary mass to halo from subhalo stars. \n
        8. Calculates top-up from berhoozi et al. \n
        9. Actually adds the top up from berhoozi et al. \n
        """
        # Total usage
        self.time_storage_array = np.empty((n_graph, n_func, 2))
        self.memory_storage_array = np.empty((n_graph, n_func, 2))
        self.process  = psutil.Process()
                                           
    def graph_timer_setup(self, n_halos):
        """
        # Reinitialise usage arrays
        """
        self.temp_time_store_array = np.empty((n_halos, self.n_func))
        self.temp_mem_store_array = np.empty((n_halos, self.n_func,2))
        return None

    def start_timer(self):
        """
        Resets CPU and memory counters
        """
        self.graph_x_start_mem = self.process.memory_info().rss
        self.graph_x_start_time = time.perf_counter_ns()
        
    def store_func_time(self, i_func, halo_no):
        """
        Stores time elapsed and memory usage since last call to start_timer
        """
        self.temp_time_store_array[halo_no, i_func-1] = (time.perf_counter_ns() - self.graph_x_start_time)
        new_mem = self.process.memory_info().rss
        # Current memory usage
        self.temp_mem_store_array[halo_no, i_func-1,0] = new_mem
        # Change since last call
        self.temp_mem_store_array[halo_no, i_func-1,1] = (new_mem - self.graph_x_start_mem)
        # Rest counters again
        self.graph_x_start_mem = new_mem
        # Down here to avoid mem_processing
        self.graph_x_start_time  = time.perf_counter_ns() 
        return None
    
    def store_average_graph_times(self, i_graph):
        """
        Averages the CPU/memory usage over all graphs and store in graph i_graph -- why would this ever be sensible?
        What is more, it seems to redfine the meaning of the last dimension.
        Uggh.
        """
        self.time_storage_array[i_graph, :, 0] = np.mean(self.temp_time_store_array, axis = 0)
        self.time_storage_array[i_graph, :, 1] = np.std(self.temp_time_store_array, axis = 0)
        self.memory_storage_array[i_graph, :, 0] = self.temp_mem_store_array[-1, :, 0]       
        self.memory_storage_array[i_graph, :, 1] = np.mean(self.temp_mem_store_array[:,:,1], axis = 0)
        return None

    def save_timing_stats(self, output_file_path, file_name):
        """
        Write timing stats to disk
        """
        file_name = file_name.split('.hdf')[0].split('/')[-1]
        save_file = '{}Timing Data {}'.format(output_file_path, file_name)
        np.save(save_file, self.time_storage_array)
               
        # Calculate the mem percentage and save memory data       
        mem_percentage = (self.memory_storage_array[:,:,0] / psutil.virtual_memory().total) * 100
        # The following is a very clumsy way of adding an extra column
        shape = list(self.memory_storage_array.shape)
        shape[-1] = shape[-1] + 1
        z = np.zeros(tuple(shape))
        z[:,:,:2] = self.memory_storage_array
        z[:,:,-1] = mem_percentage
        np.save(save_file, z)
        return None
       
    def plot_timing_barchart(self,output_file_path, file_name, save_file_path):
        """
        Load in data previously saved to disk and produce some bar plots.
        """
        file_name = file_name.split('.hdf')[0].split('/')[-1]       
        plotting_data = np.load('{}Timing Data {}.npy'.format(output_file_path, file_name))
        # A complicated way of producing natural numbers
        labels = np.arange(1,len(plotting_data[0,:,0])+1,1)
        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
        ax1.bar(labels,np.nanmean(plotting_data[:,:,0],axis= 0))
        ax1.set_xticklabels(labels)
        ax1.set_xticks(labels)
        ax1.set_ylabel('Mean time taken [ns] over all graphs with sub-halos.')
        ax1.grid(True, alpha = 0.6)
        ax1.set_axisbelow(True)
        ax1.set_xlabel('Function Measured')
        ax2.set_axis_off()
        ax2.text(0.1,0.1,self.function_descs)
        plt.tight_layout()
        plt.show()
        plt.savefig('{} Bar graph for {}.jpg'.format(save_file_path, 
                                                           file_name), dpi = 600)
        return None
        
    def plot_memory_barchart(self,output_file_path, file_name, save_file_path):
        """
        Load in data previously saved to disk and produce some bar plots.
        """
        file_name = file_name.split('.hdf')[0].split('/')[-1]
        plotting_data = np.load('{}Memory Data {}.npy'.format(output_file_path,  file_name))
        # A complicated way of producing natural numbers        
        labels = np.arange(1,len(plotting_data[0,:,1])+1,1)
        # Exclude negative memory usage (how is this possible?)
        plotting_data = np.where((plotting_data[:,:,1] < 0),np.nan,plotting_data[:,:,1])
        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
        ax1.bar(labels,np.nanmean(plotting_data,axis= 0))
        ax1.set_xticklabels(labels)
        ax1.set_xticks(labels)        
        ax1.set_ylabel('Mean memory used [B] over all graphs with sub-halos.')
        ax1.grid(True, alpha = 0.6)
        ax1.set_axisbelow(True)
        ax1.set_xlabel('Function Measured')
        ax2.set_axis_off()
        ax2.text(0.1,0.1,self.function_descs)
        plt.tight_layout()
        plt.show()
        plt.savefig('{} Memory Bar graph for {}.jpg'.format(save_file_path, 
                                                           file_name), dpi = 600)
        return None
