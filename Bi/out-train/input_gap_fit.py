# Input parameters for teach_sparse, to be parsed by run_teach_sparse.py
# Most parameters are default values and don't need to be changed.
# To iterate over multiple values (i.e. for convergence tests), make a list.
# doc: https://libatoms.github.io/QUIP/teach_sparse.html#

# External input (executable, files, external potential)
external = {
        'gap_fit': '~/.local/bin/gap_fit',
        'at_file' : 'train.xyz',  # training database
        #'core_param_file' : 'pairpot_Si_Bi.xml',  # external potential file
        #'core_ip_args' : '{IP Glue}',  # external potential type
        }

validation = {
        'quip' : '~/.local/bin/quip',
        'testxyz': 'test.xyz',
        }


# Global hyperparameters, keyword: global
# intrinsic error of atomic/bond energy, used to regularise
# the sparse covariance matrix
general = {
        'sparse_jitter' : '1e-8',  # small number used for regularisation
        'do_copy_at_file': 'False',  # copy at_file into the GAP XML file or not
        'gp_file': 'gap_Bi-Si.xml',  # name of output xml file
        'rnd_seed': '999',  # seed
        # default errors in training db energies, forces, virials, hessians
        'default_sigma' : '{0.001 0.04 0.1 0.0}',
        'virial_parameter_name': 'virial_fit',
        'e0_method':'isolated',
        'default_kernel_regularisation_local_property':'0.001', 
        'default_local_property_sigma':'0.001',
        'hessian_displacement':'1.0e-2',
        'hessian_delta':'1.0e-2',
        'config_type_sigma':'{isolated_atom:0.0001:0.04:0.01:0.0:Bi_Si_dimer:0.01:4.0:1.0:0.0}',
        }



# DESCRIPTORS
descriptors = {}

descriptor = 'distance_2b'
descriptors[descriptor] =[
        {
        'cutoff':'5.0', 
        'covariance_type':'ard_se',
        'delta':'10.0',
        'theta_uniform':'1.0',
        'sparse_method':'uniform',
        'n_sparse':'20',
        'print_sparse_index':'sparse_indices_2b.out',
        'Z1':'14',
        'Z2':'83',
        'add_species':'F'
        },        
        {
        'cutoff':'5.0', 
        'covariance_type':'ard_se',
        'delta':'10.0',
        'theta_uniform':'1.0',
        'sparse_method':'uniform',
        'n_sparse':'20',
        'print_sparse_index':'sparse_indices_2b.out',
        'Z1':'83',
        'Z2':'83',
        'add_species':'F'
        },
        ]

descriptor = 'eam_density'
descriptors[descriptor] =[
        {
        'cutoff' : '5.0',
        'pair_function' : 'FSgen' ,
        'order' : '3',
        'rmin': '0.0',
        'covariance_type' : 'ard_se',
        'delta' : '1.0',
        'theta_uniform' : '1.0',
        'sparse_method' : 'uniform',
        'n_sparse' : '20',
        'print_sparse_index' : 'sparse_indices_eam.out', 
        'Z':'83',
        'add_species': 'F' 
        },        
        ]

descriptor = 'angle_3b'
descriptors[descriptor] =[
        {
        'cutoff' : '4.5',
        'cutoff_transition_width' : '0.6',
        'n_sparse' : '600',
        'covariance_type' : 'ard_se',
        'delta' : '0.1',
        'theta_uniform' : '1.0',
        'sparse_method' : 'uniform',
        'print_sparse_index' : 'sparse_indices_3b.out',
        'Z':'14',
        'Z1':'83',
        'Z2':'83',
        'add_species':'F'
        },    
        {
        'cutoff' : '4.5',
        'cutoff_transition_width' : '0.6',
        'n_sparse' : '600',
        'covariance_type' : 'ard_se',
        'delta' : '0.1',
        'theta_uniform' : '1.0',
        'sparse_method' : 'uniform',
        'print_sparse_index' : 'sparse_indices_3b.out',
        'Z':'83',
        'Z1':'14',
        'Z2':'83',
        'add_species':'F'
        },
        {
        'cutoff' : '4.5',
        'cutoff_transition_width' : '0.6',
        'n_sparse' : '600',
        'covariance_type' : 'ard_se',
        'delta' : '0.1',
        'theta_uniform' : '1.0',
        'sparse_method' : 'uniform',
        'print_sparse_index' : 'sparse_indices_3b.out',
        'Z':'83',
        'Z1':'14',
        'Z2':'14',
        'add_species':'F'
        },
        {
        'cutoff' : '4.5',
        'cutoff_transition_width' : '0.6',
        'n_sparse' : '600',
        'covariance_type' : 'ard_se',
        'delta' : '0.1',
        'theta_uniform' : '1.0',
        'sparse_method' : 'uniform',
        'print_sparse_index' : 'sparse_indices_3b.out',
        'Z':'14',
        'Z1':'83',
        'Z2':'14',
        'add_species':'F'
        },
        {
        'cutoff' : '4.5',
        'cutoff_transition_width' : '0.6',
        'n_sparse' : '600',
        'covariance_type' : 'ard_se',
        'delta' : '0.1',
        'theta_uniform' : '1.0',
        'sparse_method' : 'uniform',
        'print_sparse_index' : 'sparse_indices_3b.out',
        'Z':'83',
        'Z1':'83',
        'Z2':'83',
        'add_species':'F'
        },
        ]
