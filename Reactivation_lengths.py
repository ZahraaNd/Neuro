import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str, required=True)
    parser.add_argument('--output_file', type=str, required=mTrue)
    args = parser.parse_args()

    data = np.load(args.input_file)[solh2] #input file name here

    
    def spikes_to_tuples(spike_array):
        spike_events = []
        for neuron_id, row in enumerate(spike_array):
            for time, value in enumerate(row):
                if not np.isnan(value):
                    spike_events.append((time, neuron_id))
        return spike_events
  
    def detect_spike_propagation_forw(spike_data, time_window, neighbor_distance=2): # the default neighboring distance was set to 2 meaning it only looks for the first neighbor's spikes
        # Assuming spike_data is a list of tuples (time, neuron_id) for blue dots
        # Sort data by time
        spike_data = sorted(spike_data)
        propagation_events = []
    
        for i, (t, neuron_id) in enumerate(spike_data):
            # Look for spikes in neighboring neurons within the time window
            forward_spikes = []
            current_time = t
    
            # Check subsequent spikes
            j = i + 1
            counter = 0
            while j < len(spike_data) and spike_data[j][0] - current_time <= time_window:
                next_time, next_neuron = spike_data[j]
                # Check if the spike is from a neighboring neuron
                second_numbers = [z[1] for z in forward_spikes]  # neuron ids for neurons already included in the propagation list
                if next_neuron - neuron_id == 1 and next_neuron not in second_numbers:  # or less than neighbor distance if you want to consider non-immediate neighbors
                    if counter == 0:
                        seed_f = neuron_id
                    counter += 1
                    forward_spikes.append((next_time, next_neuron))
                    neuron_id, current_time = next_neuron, next_time
                j = j + 1
    
            # If we found a sequence of neighboring spikes
            if len(forward_spikes) > 1:
                propagation_events.append({
                    'start_time': t,
                    'start_neuron': seed_f,
                    'propagation': forward_spikes,
                    'type': 'forward'
                })
                
        return propagation_events
    
    def detect_spike_propagation_back(spike_data, time_window, neighbor_distance=2):
        # Assuming spike_data is a list of tuples (time, neuron_id) for blue dots
        # Sort data by time
        spike_data = sorted(spike_data)
        propagation_events = []
    
        for i, (t, neuron_id) in enumerate(spike_data):
                
            backward_spikes = []
            current_time = t
            
            j = i + 1
            counter = 0
            while j < len(spike_data) and spike_data[j][0] - current_time <= time_window:
                next_time, next_neuron = spike_data[j]
                # Check if the spike is from a neighboring neuron
                second_numbers = [z[1] for z in backward_spikes]  # neuron ids for neurons already included in the propagation list, to count a neuron spiking multiple times only once
                if next_neuron - neuron_id == -1 and next_neuron not in second_numbers:  # or less than neighbor distance if you want to consider non-immediate neighbors
                    if counter == 0:
                        seed_b = neuron_id
                    counter += 1
                    backward_spikes.append((next_time, next_neuron))
                    neuron_id, current_time = next_neuron, next_time
                j = j + 1
    
            
            if len(backward_spikes) > 1:
                propagation_events.append({
                    'start_time': t,
                    'start_neuron': seed_b,
                    'propagation': backward_spikes,
                    'type': 'backward'
                })
    
        return propagation_events
    
    
    # apply
    # spike_data should be a list of (time, neuron_id) tuples extracted from the blue dots
    # You would need to preprocess the image to get these coordinates
    
    def is_subset(list1, list2):
        return set(list1).issubset(set(list2))
    
    
    def find_subset_dicts(dicts):
        subset_indices = []
        for i, dict1 in enumerate(dicts):
            for j, dict2 in enumerate(dicts):
                if i != j and is_subset(dict1['propagation'], dict2['propagation']):
                    subset_indices.append(i)
                    break
        return subset_indices
    
    
    starts = [[] for i in range(7)] # wave initiation sites
    forward = [[] for i in range(7)]
    backward = [[] for i in range(7)]

    for i in range(7):
    
        spike_array = solh2[0][i][np.array([142, 158, 143, 129, 115, 101, 116, 102, 90, 103, 91, 80, 92, 81, 72, 62, 52, 63, 73, 83, 95, 107, 122, 108, 96, 109, 124, 137, 125, 110, 126])]
        spike_data = spikes_to_tuples(spike_array)

        # data prepping

        propagation_events_forw = detect_spike_propagation_forw(spike_data, time_window= 500, neighbor_distance=1)
        propagation_events_back = detect_spike_propagation_back(spike_data, time_window= 500, neighbor_distance=1)

        # Filter out the dictionaries with subset lists
        filtered_prop_forw = [d for i, d in enumerate(propagation_events_forw) if i not in find_subset_dicts(propagation_events_forw)]
        filtered_prop_back = [d for i, d in enumerate(propagation_events_back) if i not in find_subset_dicts(propagation_events_back)]

        dicts_f = filtered_prop_forw
        dicts_b = filtered_prop_back

        # Second histogram: Start values
        #starts[i] = [d['start_neuron'] for d in dicts_f]

        # Second histogram: Start values
        forward[i] = [len(d['propagation']) for d in dicts_f]

        # Second histogram: Start values
        backward[i] = [len(d['propagation']) for d in dicts_b]

    # First, find the maximum length among all sub-arrays
    max_length_s = max(len(row) for row in starts)
    max_length_f = max(len(row) for row in forward)
    max_length_b = max(len(row) for row in backward)
    
    # Create a zero-padded array
    padded_starts = np.array([np.pad(row, (0, max_length_s - len(row)), mode='constant', constant_values=0) for row in starts])
    padded_forward = np.array([np.pad(row, (0, max_length_f - len(row)), mode='constant', constant_values=0) for row in forward])
    padded_backward = np.array([np.pad(row, (0, max_length_b - len(row)), mode='constant', constant_values=0) for row in backward])

    #return all 3
    np.savez(args.output_file, forward = padded_forward, backward = padded_backward)
        

if __name__ == "__main__":
    main()




