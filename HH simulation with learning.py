#%%
import copy
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import random
import numba as nb # used to compile AMD and FCM calculations
import time
import imageio
import matplotlib as mpl
import os

anchor = 126
# turn features on or off
stdp_on = 0
syn_on = 1
noise_on = 1

# create network
net_dict = {
'N' : 240, #0
'n_exc' : 210, #1
'n_inh' : 30, #2

'p_inh_to_exc' : 1, #3
'p_inh_to_inh' : 1, #4
}

#turn to array
net_params = np.fromiter(net_dict.values(), dtype=float)

# Model Parameters
run_dict = {

'N' : 240, #0 TOTAL NB OF NEURONS

#TIME PARAMS
'h' : 0.05, #1
'T' : 85000, #2 Remember to change it in the numerator of n_steps too, parameter nb 5!!
't_syn' : 100, #3
't_skip' : 50, #4
'n_steps' : int(85000 / 0.05), #5

#HH MODEL PARAMS
'v_exc' : 0, #6
'v_inh' : -75, #7

'C_m' : 1, #8
'g_K' : 3, #9
'g_Na' : 24, #10
'g_L' : 0.02, #11

'V_k' : -90, #12
'V_Na' : 55, #13
'V_L' : -60, #14

'tau_f' : 0.5, #15
'tau_s' : 50, #16 double check with michal's

#STDP PARAMETERS: A: AMPLITUDE, AND TAU'S
'A_plus' : 0.6,#0.01, #17
'A_minus' : 0.05,#0.01, #18
'tau_plus' : 15, #19
'tau_minus' : 100, #20

#WEIGHTS
'w_inh_exc' : 0.1,#0.1, #21 inhibitory neurons
'w_exc_exc' : 2.5, #22
'w_inh_inh' : 0.2,#0.6, #23
'w_exc_inh' : 0.1,#0.5, #24

#CURRENTS
'I_exc' : 0, #25
'I_inh' : 0.5, #26
'I_it' : 4, #4, #27
'I_bias' : 2, #28

#ACH
'g_Ks_exc' : 0.6, #29
'g_Ks_inh' : 0, #30

#NOISE PARAMS
'noise_amp' :15, #31
'noise_freq' : 0.0004, #32 freq in 1/ms x1000 to Hz

#NETWORK STRUCTURE
'n_exc' : 210, #33
'n_inh' : 30, #34

#TEMPORAL STEP SIZE
'pulse_time' : 3000, #35 150 ms in timesteps!

#STDP MAS WEIGHTS
'w_max_exc' : 8, #36
'w_max_inh' : 8, #37

'I_stim' : 2, #38
'gks_var': 0 #39
}

Re = np.zeros((int(run_dict['N']),int(run_dict['N'])))
Rs = np.zeros((int(run_dict['N']),int(run_dict['N'])))
Rf = np.zeros((int(run_dict['N']),int(run_dict['N'])))

Re[:int(run_dict['n_exc']), :int(run_dict['n_exc'])] = 0.15
Re[:int(run_dict['n_exc']), int(run_dict['n_exc']):int(run_dict['N'])] = 0.08
Re[int(run_dict['n_exc']):int(run_dict['N']), int(run_dict['n_exc']):int(run_dict['N'])] = 0.15
Rs[int(run_dict['n_exc']):int(run_dict['N']), :int(run_dict['n_exc'])] = 0.05

run_params = np.fromiter(run_dict.values(), dtype=float)

def create_net(net_params):
    # Create a 10x10 2D square DiGraph
    G = nx.DiGraph()

    nodes = [(i,j) for i in range(23) for j in range(23)]
    sorted_tuples = sorted(nodes, key=lambda x: x[0] + x[1])

    # Add nodes to the DiGraph
    for i in range(len(sorted_tuples)):
        G.add_node(sorted_tuples[i], type='excitatory')

    # Create horizontal and vertical edges
    for i in range(23):
        for j in range(23):
            if i > 0:
                G.add_edge((i, j), (i - 1, j))  # Down
                G.add_edge((i - 1, j), (i, j))
            if j > 0:
                G.add_edge((i, j), (i, j - 1))  # Left
                G.add_edge((i, j - 1), (i, j))

    # Create a list of nodes to remove (below the diagonal)
    nodes_to_remove = [node for node in G.nodes() if - node[0] + 22 < node[1]]
    G.remove_nodes_from(nodes_to_remove)

    # Create a dictionary to map node positions to node numbers
    nodes_to_skip = [(1,1), (3, 1), (4, 1), (6, 1), (8, 1), (9, 1), (11, 1), (13, 1), (14, 1), (16, 1), (18, 1), (19, 1), (1, 3), (3, 3), (5, 3), (6, 3), (8,3), (10, 3), (11, 3), (13, 3), (15, 3), (17, 3), (1, 5), (2, 5), (4, 5), (5, 5), (7, 5), (8, 5), (10, 5), (13, 5), (15, 5), (1, 7), (3, 7), (4, 7), (6, 7), (7,7), (9, 7), (11, 7), (13, 7), (1, 9), (2, 9), (4, 9), (6, 9), (8, 9), (10, 9), (12, 9), (1, 11), (3, 11), (5, 11),(7, 11), (8, 11), (10, 11),(1, 13), (2, 13), (4, 13), (6, 13), (8, 13), (1, 15), (3, 15), (4, 15), (6, 15), (1, 17), (3, 17), (4, 17), (2, 19), (3, 19)]
    G.remove_nodes_from(nodes_to_skip)

    node_labels = {node: i for i, node in enumerate(G.nodes())}

    for node in nodes_to_skip:
        i = node[0]
        j = node[1]
        edges_to_remove = [((i,j),(i+1,j)), ((i-1,j),(i,j)), ((i,j),(i,j+1)), ((i,j-1),(i,j))]
        G.remove_edges_from(edges_to_remove)

    # inhibitory network
    for i in range(int(net_params[1]), int(net_params[0])):
        G.add_node(i, type='inhibitory')

    # Make inhibitory connections
    # Define the percentage of inhibitory neurons that project to excitatory neurons
    p_inh_to_exc = net_params[3]
    # Define the percentage of inhibitory neurons that project to inhibitory neurons
    p_inh_to_inh = net_params[4]

    inhibitory_neurons = list(G.nodes())[int(net_params[1]):int(net_params[0])]
    excitatory_neurons = list(G.nodes())[:int(net_params[1])]

    for node in G.nodes():

        # variable 4: percentage of excitation, currently set to two connections
        if G.nodes[node]['type'] == 'excitatory':

            proj_exc_inh_neurons = random.sample(inhibitory_neurons, 30)
            # Add edges from the excitatory neuron to the selected inhibitory neurons
            for j in proj_exc_inh_neurons:
                G.add_edge(node, j)

        elif G.nodes[node]['type'] == 'inhibitory':

            # Select a random subset of excitatory neurons to project to
            proj_inh_exc_neurons = random.sample(excitatory_neurons, round(len(excitatory_neurons) * p_inh_to_exc))
            proj_inh_inh_neurons = random.sample(inhibitory_neurons, round(len(inhibitory_neurons) * p_inh_to_inh))
            # Add edges from the inhibitory neuron to the selected inhibitory neurons
            for j in proj_inh_inh_neurons:
                if node != j:
                    G.add_edge(node, j)
                    # G.add_edge(j,i)

            # Add edges from the inhibitory neuron to the selected excitatory neurons
            for j in proj_inh_exc_neurons:
                G.add_edge(node, j)
                # G.add_edge(j,i)

    node_labels = {node: i for i, node in enumerate(list(G.nodes())[:int(net_params[1])])}

    # Adjacency matrix
    A = nx.to_numpy_array(G)

    # Specify the node of interest and number of steps
    #node_of_interest = (11,11)
    #steps = 10

    #print(nx.single_source_shortest_path_length(G, node_of_interest).items())
    # Get nodes that are exactly 10 steps away
    #result = [node for node, distance in nx.single_source_shortest_path_length(G, node_of_interest).items() if
    #          distance == steps]

    return A, node_labels

@nb.jit(nb.float64[:,:](nb.float64[:,:], nb.float64[:,:], nb.float64[:,:], nb.float64[:], nb.int64, nb.int64, nb.float64[:,:]), nopython=True)
def stdp(w, A, t_events, run_params, i, n, freq):

    # exc_exc
    # calculate frequency
    last_true_index = -1
    second_last_true_index = -1
    for l, value in enumerate(~np.isnan(t_events[i])):
        if value:
            second_last_true_index = last_true_index
            last_true_index = l
    if last_true_index != -1 and second_last_true_index != -1:
        #print(last_true_index, second_last_true_index)
        f = float(1000.0 / (last_true_index - second_last_true_index))
        #freq[0, 0] = f
        #print(f'step: {n}, node: {i}, freq: {f}')

    # !!! note: this doesn't update weights of post connection in case of unidirectional conn, it assumes bi only
    for j in np.where(A[:,i] == 1)[0]:
        if j < int(run_params[33]):
            train_j = t_events[j][(np.isnan(t_events[j]) == False)]  # pre

            dt = n - train_j[-1] # pre - post, will be negative
            # strengthen this connection and weaken it's opposite
            dw = (run_params[17] * np.exp((-dt * 0.05) / run_params[19]))*(1-w[i,j]/run_params[36])
            #print(dw)
            if f >= 6:
                #print('p')
                if w[j, i] + dw >= int(run_params[36]):
                    w[j, i] = run_params[36]
                else:
                    w[j, i] += dw

                if w[i, j] + dw >= int(run_params[36]):
                    w[i, j] = run_params[36]
                else:
                    w[i, j] += dw

            if 4.8 < f < 6:
                #print('bi')
                if w[j, i] + dw >= int(run_params[36]):
                    w[j, i] = run_params[36]
                else:
                    w[j, i] += dw

            dw = - run_params[18] * np.exp((-dt * 0.05) / run_params[20])
            #print(dw)
            if 4.8 < f < 6:
                #print('bi')

                if w[i, j] + dw <= 0:
                    w[i, j] = 0
                else:
                    w[i, j] += dw

            if 0 < f <= 4.8 :
                #print('d')
                if w[i, j] + dw <= 0:
                    w[i, j] = 0
                else:
                    w[i, j] += dw

                if w[j, i] + dw <= 0:
                    w[j, i] = 0
                else:
                    w[j, i] += dw

    return w

@nb.jit(nb.float64[:,:](nb.float64[:,:], nb.float64[:], nb.float64[:], nb.float64[:], nb.float64[:], nb.float64[:]), nopython = True)
def ode_sys(y, I_syn, I_noise, I_ext, g_Ks, run_params):
    #print('hon')
    dydt = np.zeros((int(run_params[0]),4))
    dydt[:,0] = (1.0/(1 + np.exp((y[:,3] + 53) / 7)) - y[:,0]) / (0.37 + 2.78 * 1.0/(1 + np.exp((y[:,3] + 40.5) / 6)))
    #print('tira')
    dydt[:,1] = (1.0/(1 + np.exp((-y[:,3] - 30) / 10)) - y[:,1]) / (0.37 + 1.85 * 1.0/(1 + np.exp((y[:,3] + 27) / 15)))
    dydt[:,2] = (1.0/(1 + np.exp((-y[:,3] - 39) / 5)) - y[:,2]) / 75
    # dVdt = (1/C_m) * ( - I_Na - I_K - I_Ks - I_leak + I_ext - I_syn  +  I_noise)
    # dydt[3] = (1 / run_params[8]) * (- run_params[10] * math.pow(1 + math.exp((-y[3] - 30) / 9.5), -1) ** 3 * y[0] * (y[3] - run_params[13]) - run_params[9] * y[1] ** 4 * (y[3] - run_params[12]) - g_Ks * y[2] * (y[3] - run_params[12]) - run_params[11] * (y[3] - run_params[14]) + I_ext - I_syn + I_noise)
    # less expensive version
    dydt[:,3] = (1 / run_params[8]) * (- run_params[10] * 1.0/(1 + np.exp((-y[:,3] - 30) / 9.5)) * 1.0/(1 + np.exp((-y[:,3] - 30) / 9.5)) * 1.0/(1 + np.exp((-y[:,3] - 30) / 9.5)) * y[:,0] * (y[:,3] - run_params[13]) - run_params[9] * y[:,1] * y[:,1] * y[:,1] * y[:,1] * (y[:,3] - run_params[12]) - g_Ks * y[:,2] * (y[:,3] - run_params[12]) - run_params[11] * (y[:,3] - run_params[14]) + I_ext - I_syn + I_noise)

    return dydt

@nb.jit(nb.types.Tuple((nb.float64[:,:], nb.float64[:,:], nb.float64[:,:], nb.float64[:,:]))(nb.float64[:,:], nb.float64[:,:], nb.float64[:,:], nb.float64, nb.float64[:], nb.float64[:]), nopython = True)
def I_syn_calc(w, A, t_events, n, y, run_params):

    ef_last_spike = np.zeros((int(run_params[0]), int(run_params[0])))
    es_last_spike = np.zeros((int(run_params[0]), int(run_params[0])))
    post_V = np.ones((int(run_params[0]), int(run_params[0])))
    rp = np.append(np.full((int(run_params[33]), 1), run_params[6]), np.full((int(run_params[34]), 1), run_params[7]))

    for i in range(int(run_params[0])):
        train = t_events[i][(np.isnan(t_events[i]) == False)]

        ef = np.exp(-(n - train[-1]) * run_params[1] / run_params[15]) if train.size > 0 else 0.0
        es = np.exp(-(n - train[-1]) * run_params[1] / run_params[16]) if train.size > 0 and i >= int(run_params[33]) else 0.0

        for j in range(int(run_params[0])):
            ef_last_spike[i,j] = ef
            es_last_spike[i,j] = es

        post_V[i] = y - rp[i]

    I_syn = A * w * post_V * ((Re * ef_last_spike) + (Rs * es_last_spike))

    return I_syn, post_V, ef_last_spike, es_last_spike

def movie(frame_name, movie_name, n):
    frames = []
    # nb_of_runs = 10
    gks = [0.75, 1, 1.25, 1.5]
    sol1 = solh
    stim1 = stimh
    for k in range(3):
        print(k)
        for l in range(len(stim1[n][k])):

            step = [0, 0.5, 1]
            for m in range(1, len(step)):
                # remove first index of w_ev if you're running single trials
                we = w_evh[n][k][l]
                sol = sol1[n][k]
                stim = stim1[n][k]

                # Create a mesh grid
                x = np.linspace(0, 22, 23)
                y = np.linspace(0, 22, 23)
                X, Y = np.meshgrid(x, y)

                mesh_points = list(zip(X.flatten(), Y.flatten()))
                x_excluded = []
                y_excluded = []

                for x, y in mesh_points:
                    if (x, y) not in node_labels:
                        x_excluded.append(x)
                        y_excluded.append(y)

                # Create a 2D weight array (for example, random weights)
                weight_left_f = np.zeros((X.shape[0] - 1, X.shape[1] - 1))  # Random weights
                weight_bottom_f = np.zeros((X.shape[0] - 1, X.shape[1] - 1))  # Random weights
                weight_left_b = np.zeros((X.shape[0] - 1, X.shape[1] - 1))  # Random weights
                weight_bottom_b = np.zeros((X.shape[0] - 1, X.shape[1] - 1))  # Random weights

                for i in range(weight_left_f.shape[0]):  # Loop through rows
                    for j in range(weight_left_f.shape[1]):  # Loop through columns

                        try:
                            weight_bottom_f[i, j] = we[int(node_labels[(i, j)]), int(node_labels[(i + 1, j)])]
                            weight_bottom_b[i, j] = we[int(node_labels[(i + 1, j)]), int(node_labels[(i, j)])]

                        except KeyError:
                            pass

                        try:
                            weight_left_f[i, j] = we[int(node_labels[(i, j)]), int(node_labels[(i, j + 1)])]
                            weight_left_b[i, j] = we[int(node_labels[(i, j + 1)]), int(node_labels[(i, j)])]

                        except KeyError:
                            pass

                # weight_bottom = weight_bottom[(np.isnan(weight_bottom) == False)]

                # Normalize weights to [0, 1] for color mapping

                minn = 0
                maxx = 8
                norm_weights_left_f = (weight_left_f - minn) / (maxx - minn)
                norm_weights_bottom_f = (weight_bottom_f - minn) / (maxx - minn)
                norm_weights_left_b = (weight_left_b - minn) / (maxx - minn)
                norm_weights_bottom_b = (weight_bottom_b - minn) / (maxx - minn)
                # Create a plot
                # plt.figure(figsize=(8, 8))

                # Plot edges with colors from the 2D weight array
                for i in range(weight_left_f.shape[0]):  # Loop through rows
                    for j in range(weight_left_f.shape[1]):  # Loop through columns
                        # Get the color for the right edge and bottom edge based on weights
                        left_color_f = plt.cm.Blues(norm_weights_left_f[i, j])[:3]  # Right edge color
                        bottom_color_f = plt.cm.Blues(norm_weights_bottom_f[i, j])[:3]  # Bottom edge color

                        # Define the edges of the mesh
                        x_bottom_edge = [X[i, j], X[i + 1, j]]  # Right edge
                        y_bottom_edge = [Y[i, j], Y[i + 1, j]]

                        x_left_edge = [X[i, j], X[i, j + 1]]  # Bottom edge
                        y_left_edge = [Y[i, j], Y[i, j + 1]]

                        # Plot the right edge
                        if weight_left_f[i, j] != 0:
                            plt.plot(y_left_edge, x_left_edge, color=left_color_f, linewidth=2)

                        if weight_bottom_f[i, j] != 0:
                            # Plot the bottom edge
                            plt.plot(y_bottom_edge, x_bottom_edge, color=bottom_color_f, linewidth=2)

                mask = np.ones_like(X, dtype=bool)
                for x, y in zip(x_excluded, y_excluded):
                    mask[np.where((X == x) & (Y == y))] = False

                # Set labels
                #plt.xlabel('X axis')
                #plt.ylabel('Y axis')
                #plt.title(f'dc: {dc[k]}, step : {l}')
                plt.title(f'gks_exc = 0.6, gks_inh: {gks[k]}')

                #w_asym = [5.5, 6, 6.5, 7, 7.5]
                #plt.title(f'wf = {w_asym[n]}, wb: {w_asym[k]}')

                plt.axis('equal')  # Equal scaling of axes
                plt.scatter(X[mask].flatten(), Y[mask].flatten(), color='blue', s=5, zorder=2, edgecolors='black')
                #plt.scatter(list(node_labels.keys())[int(stim[-1])][0], list(node_labels.keys())[int(stim[-1])][1], color='red', marker = '*', s=5, zorder=2, edgecolors='black')

                # Reconstructed nodes and edges
                A_reco = np.zeros((210, 210))
                A_exc = A[:210, :210]
                idx = np.array([])
                reac = np.array([])

                spiking = [[x for x in row if run_params[35]/20 * (l + step[m-1]) <= x <= run_params[35]/20 * (l + step[m])] for row in prep(sol, run_params)[0:210]]
                freq = [len(spiking[m]) for m in range(len(spiking))]
                for j in range(len(spiking)):
                    if len(spiking[j]) != 0:
                        idx = np.append(idx, j)
                for s in idx:
                    for t in idx:
                        if int(A_exc[int(s)][int(t)]) == 1:
                            A_reco[int(s)][int(t)] = 1
                            reac = np.append(reac, int(s))
                            reac = np.append(reac, int(t))

                reac = np.unique(reac)

                specific_nodes_x = []
                specific_nodes_y = []
                for node in reac:
                    specific_nodes_x.append(list(node_labels.keys())[int(node)][0])
                    specific_nodes_y.append(list(node_labels.keys())[int(node)][1])

                color_lookup = {}
                for i in range(len(reac)):
                    color_lookup.update({int(reac[i]): freq[int(reac[i])]})

                if reac.size != 0:
                    low = sorted(color_lookup.values())[0]
                    high = sorted(color_lookup.values())[-1]
                    norm = mpl.colors.Normalize(vmin=float(low), vmax=float(high), clip=True)
                    mapper = mpl.cm.ScalarMappable(norm=norm, cmap='Reds')
                    plt.scatter(specific_nodes_x, specific_nodes_y, c=[mapper.to_rgba(i) for i in color_lookup.values()],
                                s=25, zorder=3, edgecolors='black')

                    # sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=weight_left.min(), vmax=weight_left.max()))
                    # sm.set_array([])  # Only needed for older versions of Matplotlib
                    # plt.colorbar(sm, label='Weight Value')

                #plt.scatter([list(node_labels.keys())[int(stim[l])][0]], [list(node_labels.keys())[int(stim[l])][1]], color='red', s=50, marker = '*', zorder=4, edgecolors='black')

                # Show the plot
                filename = f'{frame_name}_{k},{l}.png'
                plt.savefig(filename)
                if k == 1:
                    frames.append(filename)
                plt.clf()

    with imageio.get_writer(f'{movie_name}.gif', mode='I', duration=1500) as writer:
        for filename in frames:
            image = imageio.imread(filename)
            writer.append_data(image)

    for k in range(nb_of_runs):
        for l in range(len(stim)):
            os.remove(f'{frame_name}_{k},{l}.png')

    '''
    plt.figure(figsize=(8, 8))
    for i in range(weight_left_f.shape[0]):  # Loop through rows
        for j in range(weight_left_f.shape[1]):  # Loop through columns
            # Get the color for the right edge and bottom edge based on weights
            left_color_b = plt.cm.viridis(norm_weights_left_b[i, j])[:3]  # Right edge color
            bottom_color_b = plt.cm.viridis(norm_weights_bottom_b[i, j])[:3]  # Bottom edge color

            # Define the edges of the mesh
            x_bottom_edge = [X[i, j], X[i + 1, j]]  # Right edge
            y_bottom_edge = [Y[i, j], Y[i + 1, j]]

            x_left_edge = [X[i, j], X[i, j + 1]]  # Bottom edge
            y_left_edge = [Y[i, j], Y[i, j + 1]]

            # Plot the right edge
            if weight_left_b[i, j] != 0:
                plt.plot(y_left_edge, x_left_edge, color=left_color_b, linewidth=2)

            if weight_bottom_b[i, j] != 0:
                # Plot the bottom edge
                plt.plot(y_bottom_edge, x_bottom_edge, color=bottom_color_b, linewidth=2)

    # Set labels

    plt.xlabel('X axis')
    plt.ylabel('Y axis')
    plt.title('Mesh Edges Colored by Weight Array - backward: '+str(k))
    plt.axis('equal')  # Equal scaling of axes
    # sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=weight_left.min(), vmax=weight_left.max()))
    # sm.set_array([])  # Only needed for older versions of Matplotlib
    # plt.colorbar(sm, label='Weight Value')

    # Show the plot
    plt.show()
    '''

nb_of_runs = 2
path_length = 30 #changed
low_inh_length = 7
high_inh_length = 0
keep_dc_high = 0
adaptation_on = 0
static = 1

end_at_bound = 0
timed_end = 1
alpha = 0.000004


A, node_labels = create_net(net_params)
run_params = np.fromiter(run_dict.values(), dtype=float)

c = ['b' for i in range(210)]+['r' for i in range(30)]

w0 = np.ones((int(run_params[0]), int(run_params[0])))
w0[0:int(run_params[33]), 0:int(run_params[33])] = run_params[22]
w0[int(run_params[33]):int(run_params[0]), 0:int(run_params[33])] = run_params[21]
w0[0:int(run_params[33]), int(run_params[33]):int(run_params[0])] = run_params[24]
w0[int(run_params[33]):int(run_params[0]), int(run_params[33]):int(run_params[0])] = run_params[23]
w0 = w0 * A

# used for different modalities of using the function
#pre_path = np.array([198, 178, 166, 148, 167, 149, 136, 123, 109, 124, 137, 125, 110, 126, 111, 99, 87, 77, 67, 57, 48])
pre_path = np.array([142, 158, 143 ,129, 115, 101, 116, 102, 90, 103, 91, 80, 92, 81, 72, 62, 52, 63, 73, 83, 95, 107, 122, 108, 96, 109, 124, 137, 125, 110, 126])

# select path
p = pre_path[0: int(path_length)+1]
if low_inh_length !=0:
    p = np.concatenate((pre_path[0: int(path_length)], np.full((low_inh_length, ), pre_path[int(path_length)])))
if high_inh_length !=0:
    p = np.concatenate((pre_path[0: int(path_length)], np.full((high_inh_length, ), pre_path[int(path_length)])))

# by hand strengthen path
for i in range(len(p) - 1):
    w0[p[i], p[i + 1]] = 7.5
    w0[p[i + 1], p[i]] = 7.5

p2 = np.full(p.shape, anchor)

force_path = 1
specified_start = 0 # can't be on with a forced path, includes choosing the start point only

path_length = 0
low_inh_length = 1
high_inh_length = 0

@nb.jit(nb.types.Tuple((nb.float64[:, :], nb.float64, nb.float64[:, :], nb.float64[:], nb.float64[:],
                            nb.float64[:, :, :, :], nb.float64[:, :, :], nb.int32, nb.float64[:, :]))(nb.float64[:, :],
                                                                                                      nb.int32,
                                                                                                      nb.int32,
                                                                                                      nb.float64[:],
                                                                                                      nb.int64,
                                                                                                      nb.float64[:, :],
                                                                                                      nb.float64[:],
                                                                                                      nb.int32[:],
                                                                                                      nb.int32[:],
                                                                                                      nb.float64[:, :,
                                                                                                      :, :]),
            nopython=True)
def solve_RK4(w, counter, reward_node, occupancy, start, A, run_params, stimh, prepath, w_ev):

    print(f'run {counter}')
    frequency = np.zeros((210, int(run_params[5])))

    I_ext = np.append(np.full((int(run_params[33]), 1), run_params[25]),
                      np.full((int(run_params[34]), 1), run_params[26]))
    g_Ks = np.append(np.full((int(run_params[33]), 1), run_params[29]),
                     np.full((int(run_params[34]), 1), run_params[30]))

    w_ev[0] = w

    stim = np.zeros(1)
    if specified_start == 1:
        stim[0] = int(start)
    elif force_path == 1:
        stim[0] = int(prepath[0])
    else:
        stim[0] = 198

    Ec = 0
    Ic = 0

    noise = np.zeros((int(run_params[0]), int(run_params[5])))

    # noise of the specified frequency for each neuron
    for i in range(int(run_params[0])):
        for n in range(int(run_params[5])):
            rand_num = random.random()
            # Check if the neuron receives input
            if rand_num < run_params[32] and noise_on == 1:
                for j in range(n, (n + int(1 / run_params[1])) % run_params[5]):
                    noise[i, j] = run_params[31]

    y = np.zeros((int(run_params[5]), int(run_params[0]), 4))
    t_events = np.full((int(run_params[0]), 1), np.nan)

    k1 = np.zeros((int(run_params[0]), 4))
    k2 = np.zeros((int(run_params[0]), 4))
    k3 = np.zeros((int(run_params[0]), 4))

    # variable 8: initial conditions of neurons
    for i in range(int(run_params[0])):
        y[0, i, 0] = np.round(random.uniform(0.45, 0.9), 2)  # h0
        y[0, i, 1] = np.round(random.uniform(0.3, 0.6), 2)  # n0
        y[0, i, 2] = np.round(random.uniform(0.03, 0.07), 2)  # z0
        y[0, i, 3] = np.round(random.uniform(-65, -80), 2)  # V0

    dont_move = 0
    I_edit = 0
    I_stop = -40
    stimulate = 1
    n_last = 0
    for n in range(int(run_params[5]) - 1):

        I_syn = np.zeros((4, int(run_params[0]), int(run_params[0])))
        I_syn_s = np.zeros((4, int(run_params[0])))

        # w = w - alpha/w
        if n >= int(run_params[3] / run_params[1]) and syn_on == 1:
            I_syn[0] = I_syn_calc(w, A, t_events, n, y[n, :, 3], run_params)[0]
            I_syn[1] = I_syn_calc(w, A, t_events, n + 0.5, y[n, :, 3] + k1[:, 3] / 2.0, run_params)[0]
            I_syn[2] = I_syn_calc(w, A, t_events, n + 0.5, y[n, :, 3] + k2[:, 3] / 2.0, run_params)[0]
            I_syn[3] = I_syn_calc(w, A, t_events, n + 1, y[n, :, 3] + k3[:, 3], run_params)[0]

        if n == run_params[35]:
            w_ev[counter, 1] = w

        if n % run_params[35] == 0 and n > int(run_params[3] * 20):

            # document
            occupancy[int(stim[-1])] += 1
            w_ev[counter, int(n / run_params[35])] = w
            print(f'step: {int(n / run_params[35])}')

            # switch stimulation

            if force_path == 1:

                if low_inh_length != 0 and n // run_params[35] >= path_length:
                    print(n // run_params[35])
                    #g_Ks[int(run_params[33]):int(run_params[0])] = run_params[39]
                    I_stop = 0
                    if n // run_params[35] == path_length + low_inh_length:
                        break
                elif high_inh_length != 0 and run_params[35] >= path_length:
                    I_stop = 0
                    if n // run_params[35] == path_length + high_inh_length:
                        break
                else:
                    if n // run_params[35] == path_length:
                        break

                stim = np.append(stim, prepath[int(n // run_params[35])])
                print(stim)

        #I_ext[int(stim[-3])] = run_params[25]
        # keep dc high, keeps dc high along stimh, or p
        if keep_dc_high == 1:
            for g in stimh:
                I_ext[int(g)] = run_params[38]

        # move stimulating current one step
        I_ext[int(stim[-1])] = run_params[27] if stimulate == 1 else I_edit
        #I_ext[int(stim[-2])] = run_params[25] + I_stop
        #print(I_ext[int(198)])

        for i in range(int(run_params[0])):
            I_syn_s[0] = np.sum(I_syn[0], axis=0)
            I_syn_s[1] = np.sum(I_syn[1], axis=0)
            I_syn_s[2] = np.sum(I_syn[2], axis=0)
            I_syn_s[3] = np.sum(I_syn[3], axis=0)

        k1 = run_params[1] * ode_sys(y[n], I_syn_s[0], noise[:, n], I_ext, g_Ks, run_params)
        k2 = run_params[1] * ode_sys(y[n] + k1 / 2.0, I_syn_s[1], noise[:, n], I_ext, g_Ks, run_params)
        k3 = run_params[1] * ode_sys(y[n] + k2 / 2.0, I_syn_s[2], noise[:, n], I_ext, g_Ks, run_params)
        k4 = run_params[1] * ode_sys(y[n] + k3, I_syn_s[3], noise[:, n], I_ext, g_Ks, run_params)

        y[n + 1] = y[n] + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        # variable 7: threshold of spike recording
        a = np.full((int(run_params[0]), 1), np.nan)
        for i in range(int(run_params[0])):
            if y[n, i, 3] < 5.0 and y[n + 1, i, 3] > 5.0:
                # print('8')
                a[i, 0] = n
                if stdp_on == 1 and n > int((run_params[3] + run_params[4]) / run_params[1]) and a[i, 0] == n and i < int(run_params[33]):
                    w = stdp(w, A, t_events, run_params, i, n, frequency)
        t_events = np.column_stack((t_events, a))
        for i in range(4):
            Ec += np.sum(I_syn[i][0:int(run_params[33])][0:int(run_params[33])])
            Ic += np.sum(I_syn[i][int(run_params[33]):int(run_params[0])][0:int(run_params[33])])

    EI_b = abs(Ec / Ic) if Ic != 0 else 0
    return t_events, EI_b, w, stim, occupancy, w_ev, y, reward_node, frequency


def prep(sol, run_params):

    copi = [[] for i in range(len(sol))]
    for i in range(len(sol)):
        copi[i] = [x*0.05 for x in sol[i] if x != np.nan and x*0.05 > (run_params[3] + run_params[4])]

    return copi

#w_h = [5.5, 6.5, 7.5, 8.5]
#w_l = [4.5, 3.5, 2.5, 1.5]
gks = [0]

nb_of_runs = len(gks)
nb_of_trials = len(gks)

# prep high
eih = [np.zeros(nb_of_runs) for j in range(nb_of_trials)]
wh = [[[] for i in range(nb_of_runs)] for j in range(nb_of_trials)]
solh = [[[] for i in range(nb_of_runs)] for j in range(nb_of_trials)]
stimh = [[np.zeros(1, dtype=int) for i in range(nb_of_runs)] for j in range(nb_of_trials)]
w_evh = np.zeros((nb_of_trials, nb_of_runs, len(p) + 1, int(run_params[0]), int(run_params[0])))
occupancyh = np.zeros(run_dict['N'])
start_pts = [195, 195, 195, 196, 196, 196, 197, 197, 197, 198, 198, 198, 199, 199, 199, 200, 200, 200, 201, 201, 201, 202, 202, 202, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 11, 12, 13, 14, 15,16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]
reward_node = 0

#initialization
for j in range(nb_of_trials):
    for i in range(nb_of_runs):
        stimh[j][i][-1] = 198
        wh[j][i] = w0

for z in range(len(gks)):
    run_params[29] = gks[z]# exc
    for g in range(len(gks)):
        run_params[30] = gks[g] # inh

        w_new = w0
        t_start = time.time()
        # stimh[i-1][-1], this instead of start_pts[i] as argument to do back_to_back
        solz = solve_RK4(w_new, g, reward_node, occupancyh, start_pts[g], A, run_params, p, p2, w_evh[z]) # p is where dc is changed, the actual path, p2 is what is followed during the run
        solh[z][g], eih[z][g], wh[z][g], stimh[z][g], occupancyh, w_evh[z], vt, reward_node, f = solz
        t_end = time.time()
        duration = t_end - t_start
        print(f'run {g} duration: {int(duration // 3600):02d}:{int((duration % 3600) // 60):02d}:{int(duration % 60):02d}')

        plt.eventplot(prep(solh[z][g], run_params), colors=c)
        plt.title(f'high inh, run: {g}')
        plt.show()

        print('reached'+str(stimh[z][g][-1]))
