 format  1  0

MONOMERS
N_monomer
              2
kuhn
  1.0000000E+00  1.0000000E+00

CHAINS
N_chain
              1
N_block
              2
block_monomer
              1              2
block_length
  2.50000000E-01  7.5000000E-01

COMPOSITION
ensemble
              0
phi_chain
  1.0000000E+00

INTERACTION
interaction_type
        'chi'
chi
  2.0000000E+01

UNIT_CELL
dim
              3
crystal_system
        'cubic'
N_cell_param
              1
cell_param                              
   1.9254998725

DISCRETIZATION
ngrid
             10           10           10
chain_step
  1.0000000E-02

BASIS
group_name
      'I m -3 m'

OUTPUT_GROUP
output_filename
          'out/group_f'

OUTPUT_GROUP_CPP
output_filename
          'out/group_c'

OUTPUT_WAVES
output_filename
          'out/waves'

ITERATE
input_filename
           'in/w.bf'
output_prefix
               'out/'
max_itr
             50
error_max
  1.0000000E-08
domain
              T
itr_algo
           'NR'
N_cut
            50

FINISH
