### Tensoring Mulimodules over the torus algebra
import numpy

#updated Dec, 2013: added function to WeightedTree class that checks if a tree is negative definite


# the following two functions are used for multiplication in the torus algebra,
# where elements are represented by strings consisting of the digits 1, 2, or 3
# (for example, \rho_{12} is represented by the string '23'
#
# The first function inputs two algebra elements and returns the product, or
# None if the product is zero. The second function does the same for more than two inputs

def torus_product(str1, str2):
    '''str1,2 should be '', '1', '2', '3', '12, '23', or '123'
    returns concatenation if it makes sense, None otherwise
    '''
    if str1 == '':
        return str2
    if str2 == '':
        return str1
    if (str1, str2) in [('1', '2'), ('1', '23'), ('2', '3'), ('12', '3')]:
        return str1 + str2
    return None

def torus_product_list(elements):
    '''elements is list, each element should be '', '1', '2', '3', '12, '23', or '123'
    returns concatenation if it makes sense, None otherwise
    '''
    result = elements.pop(0)
    while len(elements) > 0:
        next_element = elements.pop(0)
        result = torus_product(result, next_element)
        if result == None:
            return None
    return result



#####      The class Multimodule is used to store bordered Heegaard Floer multimodules over the torus algebra
##### It stores the generators and differential operations, as well as some other data relevant for gluing pieces together
##### It contains several useful funtions, such as canceling differentials (remove_differentials)
##### There will be a function, tensor, which takes two Multimodules as arguments and returns their box tensor product

class Multimodule:
    def __init__(self, generators, boundary_types, fiber_directions, idempotence, operations):
        '''generators: list of strings, the names of the generators
        boundary types: list of 'D' or 'A', specifying type of each boundary
        fiber direction: list of 1 or 2 for each boundary, specifying which alpha arc is the fiber at each boundary component
                                (Note: this only makes sense in the context of S^1 bundles making up a graph manifold)
        idempotence: dictionary, keys are generators, entries are list of either 1 or 2 for each boundary component
                    1 if that boundary is type A and has alpha1 occupied or is type D and has alpha2 occupied
                    2 if that boundary is type A and has alpha2 occupied or is type D and has alpha1 occupied
        operations: list of lists- first and last entries are initial and final generators, in between
                    are entries for each boundary component: strings for type D outputs and lists of strings for type A inputs
        '''
        self.generators = generators
        self.boundary_types = boundary_types
        self.fiber_directions = fiber_directions
        self.idempotence = idempotence
        self.operations = operations


    def print_ops(self):
        for op in self.operations:
            print op

    def operations_starting_with(self, generator):
        '''returns a list of all operations with a given initial generator'''
        result = []
        for op in self.operations:
            if op[0] == generator:
                result.append(op)
        return result

    def operations_ending_with(self, generator):
        '''returns a list of all operations with a given final generator'''
        result = []
        for op in self.operations:
            if op[-1] == generator:
                result.append(op)
        return result
    
    def reduce_mod_2(self):
        '''removes repeated operations in pairs'''
        self.operations.sort()
        for i in range(len(self.operations)-1):
            if self.operations[i] == self.operations[i+1]:
                self.operations[i] = 'REMOVE'
                self.operations[i+1] = 'REMOVE'
        while 'REMOVE' in self.operations:
            self.operations.remove('REMOVE')

    def relabel_generators(self):
        if len(self.generators) < 26:
            names = 'abcdefghijklmnopqrstuvwxyz'
        else:
            names = []
            for i in range(len(self.generators)):
                names.append('g' + str(i))
        translate = {}
        new_generators = []
        new_idempotence = {}
        for i in range(len(self.generators)):
            translate[self.generators[i]] = names[i]
            new_generators.append(names[i])
            new_idempotence[names[i]] = self.idempotence[self.generators[i]]

        self.generators = new_generators
        self.idempotence = new_idempotence
        for op in self.operations:
            op[0] = translate[ op[0] ]
            op[-1] = translate[ op[-1] ]
            

### the following functions are related to simplifying the differential
    def find_differentials(self):
        '''returns a list of all operations with no Reeb chords.
        #assumes all boundary components are type D'''
        result = []
        for op in self.operations:
            keep = True
            for i in range(1, len(op)-1):
                if not op[i] in ['', []]:
                    keep = False
            if keep:
                result.append(op)

        return result

    def find_removable_differential(self):
        '''returns an operation which is a differential and such that
        there are no other operations between the same two generators.
        if no such differential exists, returns None'''
        
        differentials = self.find_differentials()
        for d in differentials:
            initial_gen = d[0]
            ops_from_init_gen = self.operations_starting_with(initial_gen)
            ops_from_init_gen.remove(d)
            if not d[-1] in [op[-1] for op in ops_from_init_gen]:
                return d

        return None


    def remove_differentials(self):
        '''removes all removable differentials, generating new operations via the zig-zag rule
        '''

        #find a differential to remove, if any exist
        differential = self.find_removable_differential()

        while not differential == None:
            initial_generator = differential[0]
            final_generator = differential[-1]


            #find list of all other operations out of initial_generator and into final_generator
            incoming_ops = self.operations_ending_with(final_generator)
            outgoing_ops = self.operations_starting_with(initial_generator)
            incoming_ops.remove(differential)
            outgoing_ops.remove(differential)


            #find new operations resulting from zig-zag method
            new_operations = []
            for op1 in incoming_ops:
                for op2 in outgoing_ops:
                    new_op = [op1[0]]
                    for i in range(len(self.boundary_types)):
                        if self.boundary_types[i] == 'D':                       #that is, if the ith boundary component is type D
                            new_op.append(torus_product(op1[i+1], op2[i+1]))        #then the algebra element for the ith boundary in the new operation is the product of the elements for the original two operations
                        else:                                                   #if the ith boundary is type A, we concatenate the lists of inputs for the new operation
                            new_op.append( op1[i+1] + op2[i+1] )                #if the ith boundary is type A, we concatenate the lists of inputs for the new operation
                    new_op.append(op2[-1])
                    if not None in new_op:                              #None appears in new_op if some product of algebra elements which appears is zero
                        new_operations.append(new_op)

            #it may be that some new ops contain initial_generator or final_generator (i.e. if there are self differentials)
            #  these should be removed, since initial_generator and final_generator will be cancelled
            new_operations2 = []
            for op in new_operations:
                if not (op[0] in [initial_generator, final_generator] or op[-1] in [initial_generator, final_generator]):
                    new_operations2.append(op)
            new_operations = new_operations2

            #all original operations which do not involve initial_generator or final generator survive
            # these are now added to the new list of operations
            for op in self.operations:
                if not (op[0] in [initial_generator, final_generator] or op[-1] in [initial_generator, final_generator]):
                    new_operations.append(op)

            #give the multimodule the updated list of operations        
            self.operations = new_operations

            #remove initial_generator and final_generator, and the corresponding entries for generator_signs
            i = self.generators.index(initial_generator)
            j = self.generators.index(final_generator)
            self.generator_signs.pop(max(i, j))
            self.generator_signs.pop(min(i, j))
            self.generators.remove(initial_generator)
            self.generators.remove(final_generator)
            
            self.reduce_mod_2()           

            differential = self.find_removable_differential()   #pick new removable differential and repeat while loop            


    def is_valid_D_sequence(self, list_of_ops, boundary):
        '''for a given list of ops, checks if the (nonempty) list of operations has nonzero
        products for all type D boundary components except the specified ones'''

        for i in range(1,boundary+1):
            if type(list_of_ops[0][i]) == str:  ## i.e., if this is a type D bdy component
                rho_list = [op[i] for op in list_of_ops]
                while '' in rho_list:
                    rho_list.remove('')
                for j in range(len(rho_list)-1):
                    if not (rho_list[j][-1], rho_list[j+1][0]) in [('1', '2'), ('2', '3')]:
                        return False

        for i in range(boundary+2,len(list_of_ops[0])-1):
            if type(list_of_ops[0][i]) == str:  ## i.e., if this is a type D bdy component
                rho_list = [op[i] for op in list_of_ops]
                while '' in rho_list:
                    rho_list.remove('')
                for j in range(len(rho_list)-1):
                    if not (rho_list[j][-1], rho_list[j+1][0]) in [('1', '2'), ('2', '3')]:
                        return False

        return True

        
    def d_squared(self):
        '''prints d^2, when all boundaries are type D
        this should always print the empty list for valid type D multimodules'''

        double_ops = []

        for op1 in self.operations:
            for op2 in self.operations:
                if op1[-1] == op2[0]:
                    outputs = [ torus_product(op1[i], op2[i]) for i in range(1, len(op1)-1) ]
                    if not None in outputs:
                        double_ops.append( [op1[0]] + outputs + [op2[-1]] )

        double_ops.sort()
        for i in range(len(double_ops)-1):
            if double_ops[i] == double_ops[i+1]:
                double_ops[i] = 'REMOVE'
                double_ops[i+1] = 'REMOVE'
        while 'REMOVE' in double_ops:
            double_ops.remove('REMOVE')

        print double_ops

    def rank(self):
        ''' for a multimodule with no boundary components (i.e. a chain complex)
        returns the number of generators after canceling differentials'''

        if not len(self.boundary_types) == 0:
            print 'Error: this is a module, not a chain complex'
        self.remove_differentials()
        return len(self.generators)

    def homology(self):
        ''' for a multimodule with no boundary components (i.e. a chain complex)
        representing a closed three-manioflds, returns the absolute value of the
        signed count of generators, which gives the size of H_1 of the manifold'''

        if not len(self.boundary_types) == 0:
            print 'Error: this is a module, not a chain complex'
        self.remove_differentials()
        return abs(sum(self.generator_signs))
    


############
####    Tensoring
############
                
# when tensoring, this function determines the new operation resulting from a type A operation
# with a corresponding sequence of type D operations
def make_new_operation(A_op, boundary1, list_of_D_ops, boundary2):
    '''list_of_D_ops must be non_empty

    makes a single operation for the tensor product from one operation coming
    from a module which is type A on boundary1, and a VALID list of operations for a
    module which is type D on boundary 2, where the type A inputs match up with
    the sequence of type D outputs, and all other type D products make sense'''

    initial_generator = (A_op[0], list_of_D_ops[0][0])
    final_generator = (A_op[-1], list_of_D_ops[-1][-1])
    new_operation = [initial_generator]

    for i in range(1, boundary1+1):         #include all the boundary algebra elements on the type A side before the given boundary
        new_operation.append(A_op[i])       
    for i in range(boundary1+2, len(A_op)-1):   #include all the boundary algebra elements on the type A side after the given boundary
        new_operation.append(A_op[i])
        
    for i in range(1, boundary2+1):
        new_sequence = [op[i] for op in list_of_D_ops]
        if type(new_sequence[0]) == list:       #if this boundary is type A
            result = []
        if type(new_sequence[0]) == str:        #if this boundary is type D
            result = ''
        for term in new_sequence:
            result += term
        new_operation.append(result)
                    
    for i in range(boundary2+2, len(list_of_D_ops[0]) - 1):
        new_sequence = [op[i] for op in list_of_D_ops]
        if type(new_sequence[0]) == list:       #if this boundary is type A
            result = []
        if type(new_sequence[0]) == str:        #if this boundary is type D
            result = ''
        for term in new_sequence:
            result += term
        new_operation.append(result)

    new_operation.append(final_generator)

    return new_operation

# if an operation has no inputs wrt the type A boundary, there are new operations coming from
# that operation tensor any generator (with appropriate idempotent) on the type D side.
# this function finds the new op given a type A op and appropriate type D generator
def make_new_operation_from_Adiff(A_op, boundary1, module2, D_gen, boundary2):
    initial_generator = (A_op[0], D_gen)
    final_generator = (A_op[-1], D_gen)
    new_op = [ initial_generator ]

    #For the unused boundaries on the type A side, just copy outputs from the differential
    for i in range(0, boundary1):
        new_op.append( A_op[i+1] )
    for i in range(boundary1+1, len(A_op)-2):
        new_op.append( A_op[i+1] )

    #For the unused boundaries on type D side, fill in with appropriate type of "empty" output
    for i in range(0, boundary2):
        if module2.boundary_types[i] == 'A':      #if boundary component i on module1 is type A
            new_op.append( [] )         
        else:
            new_op.append( '' )

    for i in range(boundary2+1, len(module2.boundary_types)):
        if module2.boundary_types[i] == 'A':      #if boundary component i on module1 is type A
            new_op.append( [] )         
        else:
            new_op.append( '' )

    new_op.append(final_generator)

    return new_op


# if an operation has no outputs wrt the type D boundary being tensored, there are new operations coming from
# that operation tensor any generator (with appropriate idempotent) on the type A side.
# this function finds the new op given a type D op and appropriate type A generator
def make_new_operation_from_Ddiff(module1, A_gen, boundary1, D_op, boundary2):
    initial_generator = (A_gen, D_op[0])
    final_generator = (A_gen, D_op[-1])
    new_op = [ initial_generator ]

    #For the unused boundaries on type A side, fill in with appropriate type of "empty" output
    for i in range(0, boundary1):
        if module1.boundary_types[i] == 'A':         #if boundary component i on module1 is type A
            new_op.append( [] )         
        else:
            new_op.append( '' )

    for i in range(boundary1+1, len(module1.boundary_types)):
        if module1.boundary_types[i] == 'A':         #if boundary component i on module1 is type A
            new_op.append( [] )         
        else:
            new_op.append( '' )

    #For the unused boundaries on the type D side, just copy outputs from the differential
    for i in range(0, boundary2):
        new_op.append( D_op[i+1] )
    for i in range(boundary2+1, len(D_op)-2):
        new_op.append( D_op[i+1] )

    new_op.append(final_generator)

    return new_op


    



def tensor(module1, boundary1, module2, boundary2, simplify=True):
    '''module1 and module2 are Multimodule instances
    boundary1 and boundary 2 are integers, indicating which
    boundary components to glue
    module1 boundary component should be type A
    module2 boundary component should be type D

    returns a new Multimodule instance
    differentials will be canceled after tensoring unless simplify=False'''


    #determine new generators and new generator signs
    new_gens = []
    new_gen_signs = []
    for i in range(len(module1.generators)):
        gen1 = module1.generators[i]
        for j in range(len(module2.generators)):
            gen2 = module2.generators[j]
            if module1.idempotence[gen1][boundary1] == module2.idempotence[gen2][boundary2]:
                new_gens.append( (gen1, gen2) )
                new_gen_signs.append( module1.generator_signs[i]*module2.generator_signs[j] )
                



    new_boundary_types = (  module1.boundary_types[:boundary1]
                           + module1.boundary_types[boundary1+1:]
                           + module2.boundary_types[:boundary2]
                           + module2.boundary_types[boundary2+1:]  )

    new_fiber_directions = (  module1.fiber_directions[:boundary1]
                           + module1.fiber_directions[boundary1+1:]
                           + module2.fiber_directions[:boundary2]
                           + module2.fiber_directions[boundary2+1:]  )

    new_idempotence = {}
    for gen in new_gens:
        new_idempotence[gen] = (   module1.idempotence[gen[0]][:boundary1]
                                 + module1.idempotence[gen[0]][boundary1+1:]
                                 + module2.idempotence[gen[1]][:boundary2]
                                 + module2.idempotence[gen[1]][boundary2+1:]  )

    
    new_operations = []
    # most new operations come from a type A operations with inputs (on the approprate boundary) of [r1, r2, ..., rn]
    # and a sequence of n type D operations with outputs r1, r2, ... rn
    #
    # the stragegy for finding these operations is to list all sequences of input for type A operations (as an intermediate step, we also list partial seequences)
    # to each sequence we associate all sequences of type D operations with the correct outputs
    # for each such sequence, we make a new operations
    
    Reeb_chains = {}    #dictonary, keys are (partial) chains of inputs for type A operations (lists), with idempotence of initial generator at beginning
                        #entries are lists of matching sequences of type D ops
    for op in module1.operations:
        Reeb_sequence = op[boundary1 + 1]       #this is the sequence of inputs for the type A operation
        initial_gen = op[0]
        initial_idem = module1.idempotence[initial_gen][boundary1]
        for i in range(1, len(Reeb_sequence)+1):
            Reeb_chains[ tuple([initial_idem] + Reeb_sequence[0:i]) ] = []
    #Reeb_chains now has an entry for every partial or complete sequences of inputs for a type A operation

#    # we also need to consider "type A differentials", where there are no inputs for the given boundary component
#    Reeb_chains[ (1, '') ] = []
#    Reeb_chains[ (2, '') ] = []

    for op in module2.operations:
        initial_gen = op[0]
        initial_idem = module2.idempotence[initial_gen][boundary2]
        Reeb_chord = op[1+boundary2]
        final_gen = op[-1]

        if (initial_idem, Reeb_chord) in Reeb_chains.keys():
            Reeb_chains[ (initial_idem, Reeb_chord) ].append( [op] )
    #all of the Reeb_chain entries with a single Reeb chord (i.e. (1, '23') ) are now filled in with all matching D op sequences

    partial_chains = Reeb_chains.keys()
    partial_chains.sort(key = len)

#    partial_chains.remove((1, ''))
#    partial_chains.remove((2, ''))

    for pc in partial_chains:
        pclist = list(pc)
        for path in Reeb_chains[pc]:        #path is a list of type D operations
            (init_gen, final_gen) = (path[0][0], path[-1][-1])
            for op in module2.operations_starting_with( path[-1][-1] ):
                if tuple(pclist + [op[boundary2 + 1]]) in partial_chains and module2.is_valid_D_sequence(path + [op], boundary2):
                    Reeb_chains[ tuple(pclist + [op[boundary2 + 1]]) ].append( path + [op] )
    #all entries for nonempty sequences of inputs now have all matching paths of type D ops

    for op in module1.operations:
        initial_Agen = op[0]
        initial_idem = module1.idempotence[initial_Agen][boundary1]
        Reeb_sequence = op[boundary1 + 1]
        final_Agen = op[-1]

        if Reeb_sequence == []:
            for Dgen in module2.generators:
                if initial_idem == module2.idempotence[Dgen][boundary2]:
                    new_operations.append( make_new_operation_from_Adiff(op, boundary1, module2, Dgen, boundary2) )
        else:
            for path in Reeb_chains[ tuple([initial_idem] + Reeb_sequence) ]:
                new_operations.append( make_new_operation(op, boundary1, path, boundary2))

    ### we also get new operations with a generator on type A side tensored with a differential on type D side
    differentials = []
    for op in module2.operations:
        if op[boundary2 + 1] == '':
            differentials.append(op)
    for diff in differentials:
        initial_D_gen = diff[0]
        for initial_A_gen in module1.generators:
            if module1.idempotence[initial_A_gen][boundary1] == module2.idempotence[initial_D_gen][boundary2]:
                ### make a new op from initial_A_gen tensor diff
                new_op = make_new_operation_from_Ddiff(module1, initial_A_gen, boundary1, diff, boundary2)
                new_operations.append( new_op)
    #now all new operations have been found

    new_module = Multimodule(new_gens, new_boundary_types, new_fiber_directions, new_idempotence, new_operations)
    new_module.reduce_mod_2()
    
    new_module.generator_signs = new_gen_signs


    if simplify:
        new_module.remove_differentials()
    #new_module.relabel_generators()

    return new_module

##################
###############  End of tensor product code
##################
##################
##################
###############  Saved multimodules
##################


######### CFAA(Id)
generators = ['w_1', 'w_2', 'x', 'y', 'z_1', 'z_2']
boundary_types = ['A', 'A']
fiber_directions = [1,2]
idempotence = {'w_1':[2,1],
               'w_2':[2,1],
               'x':[1,1],
               'y':[2,2],
               'z_1':[1,2],
               'z_2':[1,2]}
operations = [['w_1', [], ['1'], 'y'],
              ['w_1', ['2'], ['12'], 'x'],
              ['w_1', [], [], 'w_2'],
              ['w_1', ['23'], ['12'], 'w_2'],
              ['w_1', ['2'], ['123'], 'z_2'],
              ['w_1', ['2'], ['3', '2', '1'], 'z_2'],

              ['z_1', ['1'], [], 'y'],
              ['z_1', ['12'], ['2'], 'x'],
              ['z_1', [], [], 'z_2'],
              ['z_1', ['12'], ['23'], 'z_2'],
              ['z_1', ['123'], ['2'], 'w_2'],

              ['y', ['2'], ['2'], 'x'],
              ['y', ['23'], ['2'], 'w_2'],
              ['y', ['2'], ['23'], 'z_2'],

              ['x', ['3'], [], 'w_2'],
              ['x', [], ['3'], 'z_2'] ]
CFAA_id = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFAA_id.generator_signs = [1, -1, 1, -1, 1, -1]
###########


######### CFD( solid torus )    meridian alpha arc is alpha_2
generators = ['p']
boundary_types = ['D']
fiber_directions = [1]
idempotence = {'p':[2]}
operations = [ ['p', '23', 'p'] ]

CFD_2 = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFD_2.generator_signs = [1]
###########


######### CFD( solid torus )    meridian alpha arc is alpha_1
generators = ['p']
boundary_types = ['D']
fiber_directions = [2]
idempotence = {'p':[1]}
operations = [ ['p', '12', 'p'] ]

CFD_1 = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFD_1.generator_signs = [1]
###########

CFA_2 = tensor(CFAA_id,0, CFD_2, 0)
CFA_1 = tensor(CFAA_id,1, CFD_1, 0)
#Note: tensoring CFA_1 with CFD_1  or CFA_2 with CFD_2 gives S^3
#     tensoring CFA_1 with CFD_2 gives S^2 x S^1

########### CFDDD( pair-of-pants X S^1 )
generators = ['V', 'W', 'X', 'Y', 'Z']
boundary_types = ['D','D','D']
fiber_directions = [2,2,1]
idempotence = {'V': [1, 1, 1],
               'W': [1, 1, 1],
               'X': [2, 1, 1],
               'Y': [2, 2, 2],
               'Z': [1, 2, 1]}
operations = [ ['V','1'  ,'123','3'  , 'Y'],
               ['V','123','123','123', 'Y'],
               ['V', ''  , '3' , ''  , 'Z'],
               ['V', '3' , ''  , ''  , 'X'],
               
               ['W', '1' , '1' , '3' , 'Y'],
               ['W', '3' , ''  ,'12' , 'X'],
               ['W','123', '1' ,'123', 'Y'],
               
               ['X', '2' , ''  ,'12' , 'V'],
               ['X', ''  , '3' , '1' , 'Y'],
               ['X', '2' , ''  , ''  , 'W'],
               
               ['Y', ''  , '2' , '2' , 'X'],
               ['Y', '2' , ''  , '2' , 'Z'],
               
               ['Z', ''  , '2' , ''  , 'W'],
               ['Z', '3' , ''  , '1' , 'Y']   ]

CFDDD = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFDDD.generator_signs = [1, -1, 1, -1, 1]
#############


########## CFDA( positive Dehn twist about alpha1 on type A side )
        ## fiber is alpha_1
        ## This is tau_l in LOT bimodules pg 147
generators = ['x', 'y', 'z']
boundary_types = ['D','A']
fiber_directions = [1,1]
idempotence = {'x':[2, 2],
               'y':[1, 2],
               'z':[1, 1]}
operations = [ ['x', '2', ['2', '1'], 'y'],
               ['x', '2', ['2', '12'], 'z'],
               ['x', '23', ['2', '123'], 'x'],
               ['y', '1', [], 'x'],
               ['y', '', ['2'], 'z'],
               ['y', '3', ['23'], 'x'],
               ['z', '12', ['1'], 'y'],
               ['z', '12', ['12'], 'z'],
               ['z', '123', ['123'], 'x'],
               ['z', '3', ['3'], 'x']  ]
CFDA_twist1 = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFDA_twist1.generator_signs = [1, 1, 1]

        ## Negative Dehn twist about alpha1 on type A side
generators = ['x', 'y', 'z']
boundary_types = ['D','A']
fiber_directions = [1,1]
idempotence = {'x':[1, 1],
               'y':[1, 2],
               'z':[2, 2]}
operations = [ ['x', '3', ['3'], 'z'],
               ['x', '', ['1'], 'y'],
               ['x', '12', ['12'], 'x'],
               ['x', '123', ['123'], 'z'],
               ['x', '1', ['12', '1'], 'z'],
               ['y', '12', ['2'], 'x'],
               ['y', '123', ['23'], 'z'],
               ['y', '1', ['2', '1'], 'z'],
               ['z', '2', [], 'y']  ]
CFDA_twist1_inv = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFDA_twist1_inv.generator_signs = [1, -1, 1]
###########


######### CFDA(positive Dehn twist about alpha2 on type A side )
        # fiber is alpha_2  
generators = ['x', 'y', 'z']
boundary_types = ['D','A']
fiber_directions = [2,2]
idempotence = {'x':[1, 1],
               'y':[2, 1],
               'z':[2, 2]}
operations = [ ['x', '1', ['1'], 'z'],
               ['x', '123', ['12'], 'y'],
               ['x', '123', ['123'], 'z'],
               ['x', '3', ['3', '2'], 'y'],
               ['x', '3', ['3', '23'], 'z'],
               ['y', '', ['3'], 'z'],
               ['y', '2', [], 'x'],
               ['z', '23', ['2'], 'y'],
               ['z', '23', ['23'], 'z'] ]
CFDA_twist2 = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFDA_twist2.generator_signs = [1, -1, 1]

        ## negative Dehn twist about alpha2 on type A side
generators = ['x', 'y', 'z']
boundary_types = ['D','A']
fiber_directions = [2,2]
idempotence = {'x':[1, 1],
               'y':[2, 1],
               'z':[2, 2]}
operations = [ ['x', '1', ['1'], 'z'],
               ['x', '1', ['12'], 'y'],
               ['x', '123', ['123'], 'z'],
               ['x', '12', ['123', '2'], 'x'],
               ['x', '3', [], 'y'],
               ['y', '23', ['3'], 'z'],
               ['y', '2', ['3','2'], 'x'],
               ['z', '', ['2'], 'y'],
               ['z', '2', ['23', '2'], 'x'],
               ['z', '23', ['23'], 'z'] ]
CFDA_twist2_inv = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFDA_twist2_inv.generator_signs = [1, 1, 1]
################


##### CFDD_id  ########
generators = ['a', 'b']
boundary_types = ['D', 'D']
fiber_directions = [2, 1]
idempotence = {'a': [2, 2], 'b': [1, 1]}
operations = [['a', '2', '2', 'b'],
              ['b', '1', '3', 'a'],
              ['b', '123', '123', 'a'],
              ['b', '3', '1', 'a']]
CFDD_id = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFDD_id.generator_signs = [1, -1]


##################
###############  End of saved multimodules
##################
##################
##################
###############  Functions for plumbing trees
##################

# Objects in the following class, WeightedTree, are weighted acyclic graphs.
class WeightedTree:
    def __init__(self, vertex_weights, edges):
        '''vertex_weights is a list of integers, one for each vertex
        edges is a list of pairs (i, j), indicating that there is an edge
        connecting v_i to v_j (where vertices are v_0, v_1, ..., v_n)

        an edge may have the form (i, '*'), indicating that there is to be
        an open boundary component at v_i

        assumes graph is a tree'''

        self.weights = vertex_weights
        self.edges = edges
        self.edges_from_vertex = []
        for i in range(len(self.weights)):
            edges_from_i = []
            for edge in edges:
                if edge[0] == i:
                    edges_from_i.append(edge[1])
                if edge[1] == i:
                    edges_from_i.append(edge[0])
            self.edges_from_vertex.append(edges_from_i)
        #edges_from_vertex[i] gives an (ordered) list of vertices (represented by integers) connected to v_i

        self.matrix = []            #this is -1 times the adjacency matrix
        for i in range(len(self.weights)):
            line_i = [0]*len(self.weights)
            line_i[i] = -self.weights[i]
            for edge in self.edges:
                if edge[0] == i:
                    line_i[edge[1]] = -1
                if edge[1] == i:
                    line_i[edge[0]] = -1
            self.matrix.append(line_i)

        self.det = int(round(numpy.linalg.det(self.matrix)))
        if len(self.matrix) == 1:
            self.minor_det = 1
        else:
            self.minor_det = numpy.linalg.det( [line[1:] for line in self.matrix[1:]] )
            self.minor_det = int(round(self.minor_det))
        

    def is_negative_definite(self):
        '''returns true if the tree represents a negative definite plumbing
            (that is, if self.matrix, which is negative the adjacency matrix, is positive definite'''
        e = numpy.linalg.eig(self.matrix)   #the list of eigenvalues of self.matrix
        min_lambda = min(e[0])
        if min_lambda > 0:                  #if the largest eigenvalue is positive, then the matrix is positive definite: return True
            return True
        else:
            return False
                

    def remove_vertex(self, vertex):
        '''vertex = integer
        removes the given vertex and returns a list of the connected
        components of the remaining graph. the vertex that was connected
        to the given vertex is always vertex 0 in the returned sub_graphs
        '''
        result = []
        for v in self.edges_from_vertex[vertex]:
            queue = self.edges_from_vertex[v]
            queue.remove(vertex)
            sub_tree_vertices = [v]
            while len(queue) > 0:
                current_vertex = queue.pop(0)
                if not current_vertex in sub_tree_vertices:
                    queue += self.edges_from_vertex[current_vertex]
                    sub_tree_vertices.append(current_vertex)
            sub_tree_edges = []
            for edge in self.edges:
                if edge[0] in sub_tree_vertices and edge[1] in sub_tree_vertices:
                    sub_tree_edges.append(edge)
            translate = {}
            sub_tree_weights = []
            for i in range(len(sub_tree_vertices)):
                translate[sub_tree_vertices[i]] = i
                sub_tree_weights.append( self.weights[ sub_tree_vertices[i] ] )
            for i in range(len(sub_tree_edges)):
                edge = sub_tree_edges[i]
                sub_tree_edges[i] = (translate[edge[0]], translate[edge[1]])

            result.append( WeightedTree(sub_tree_weights, sub_tree_edges) )
        return result
                
            
     
# functions which compute HFhat for plumbings along WeightedTrees
def plumb_recursive(graph, vertex, fiber_alpha):
    '''takes a weighted graph and a chosen vertex
    returns the module CFD obtained by plumbing according to the graph,
    and adding one extra boundary component at the chosen vertex,
    such that fiber_alpha (1 or 2) is the alpha arc correspondng
    to the fiber over that vertex'''
   
    valence = len(graph.edges_from_vertex[vertex])
    weight = graph.weights[vertex]

    # determine twisting based on weight
    if fiber_alpha == 1:
        twist = CFDA_twist1
        twist_inv = CFDA_twist1_inv
    elif fiber_alpha == 2:
        twist = CFDA_twist2
        twist_inv = CFDA_twist2_inv
    if weight == 0:
        twisting = tensor(twist, 1, twist_inv, 0)
    if weight > 0:
        twisting = twist
        for i in range(weight - 1):
            twisting = tensor(twist, 1, twisting, 0)
    if weight < 0:
        twisting = twist_inv
        for i in range(-weight - 1):
            twisting = tensor(twist_inv, 1, twisting, 0)
    # found twist
    
    #base case
    if valence == 0:
        if fiber_alpha == 1:
            result = tensor(twisting, 1, CFD_2, 0)
            return result
        else:
            result = tensor(twisting, 1, CFD_1, 0)
            return result

    #find vertex_module, the multimodule associated to the S^1 bundle at the chosen vertex
    if valence == 1:
        vertex_module = CFDD_id
    if valence > 1:
        vertex_module = CFDDD
        #boundary 0 of vertex_module should always have alpha_fiber 1
        for i in range(valence - 2):
            add_on = tensor(CFAA_id, 0, CFDDD, 2)
            vertex_module = tensor(add_on, 0, vertex_module, 0)

    chosen_boundary = vertex_module.fiber_directions.index(fiber_alpha)
    vertex_module = tensor(twisting, 1, vertex_module, chosen_boundary)
    #now boundary 0 of vertex_mfd is the one I want to keep to the end
    #it has the correct alpha_fiber

    #construct a list of components of the graph minus the distinguished vertex
    sub_trees = graph.remove_vertex(vertex)
    #for each component, recursively call this function to determin CFD
    while len(vertex_module.boundary_types) > 1:
        sub_tree = sub_trees.pop(0)
        desired_fiber_alpha = vertex_module.fiber_directions[1]
        cap = plumb_recursive(sub_tree, 0, desired_fiber_alpha)
        
        if desired_fiber_alpha == 1:
            vertex_module = tensor(CFAA_id, 0, vertex_module, 1)
            vertex_module = tensor(vertex_module, 0, cap, 0)
        if desired_fiber_alpha == 2:
            vertex_module = tensor(CFAA_id, 1, vertex_module, 1)
            vertex_module = tensor(vertex_module, 0, cap, 0)

    return vertex_module


def plumb(graph):
    '''takes a weighted tree (an instance of WeightedTree), and returns
    the chain complex CFhat for the 3-manifold obtained by plumbing along the tree
    all differentials are removed, so len(result.generators) gives the rank of HFhat
    '''
    
    #special case for graph = pt
    if len(graph.weights) == 1:
        weight = graph.weights[0]
        result = CFD_1
        if weight > 0:
            for i in range(weight):
                result = tensor(CFDA_twist2, 1, result, 0)
        if weight < 0:
            for i in range(-weight):
                result = tensor(CFDA_twist2_inv, 1, result, 0)
        result = tensor(CFA_2, 0, result, 0)
        
    else:       #that is, for more than one vertex
        #first we find an end of the tree, a velence 1 vertex
        valences = [len(graph.edges_from_vertex[i]) for i in range(len(graph.weights))]
        first_end = valences.index(1)
        weight = graph.weights[first_end]

        #then we call plumb_recursive for the tree with that 
        result = CFD_1
        if weight > 0:
            for i in range(weight):
                result = tensor(CFDA_twist2, 1, result, 0)
        if weight < 0:
            for i in range(-weight):
                result = tensor(CFDA_twist2_inv, 1, result, 0)
        result = tensor(CFAA_id, 1, result, 0)

        new_graph = graph.remove_vertex(first_end)[0]
        cap = plumb_recursive(new_graph, 0, 2)

        
        result = tensor(result, 0, cap, 0)

    return result
    
