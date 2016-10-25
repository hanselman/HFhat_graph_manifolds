# Computing HF_hat of graph manifolds
#
#   Last upated 1/10/14




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


def list_sum(list_of_lists):
    '''takes a list of lists and returns the single list obtained by concatenating them'''
    
    result = list_of_lists.pop(0)
    while len(list_of_lists) > 0:
        next_list = list_of_lists.pop(0)
        result = result + next_list
    return result


def connected_components(list_of_ops):
    '''inputs list of operations
    returns list of lists of operations, each inner list
    is a connected component of the graph given by operations
    as arrows connecting generators'''
    
    components = []
    for op in list_of_ops:
        gen1 = None
        gen2 = None
        for i in range(len(components)):
            if op[0] in components[i]:
                gen1 = i
            if op[-1] in components[i]:
                gen2 = i
        if gen1 == None and gen2 == None:
            components.append([op[0], op[-1]])
        elif gen1 == None:
            components[gen2].append(op[0])
        elif gen2 == None:
            components[gen1].append(op[-1])
        elif gen1 != gen2:
            comp1 = components.pop(gen1)
            comp2 = components.pop(gen2 - (gen1<gen2))
            components.append(comp1+comp2)

    component_ops = [[] for comp in components]
    for op in list_of_ops:
        for i in range(len(components)):
            if op[0] in components[i]:
                component_ops[i].append(op)

    return component_ops

        
        

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
        self.name = 'no_name'   #default name, can be replaced by a string when a module is defined


    def print_ops(self):
        for op in self.operations:
            print op

    def print_ops2(self):
        for comp in connected_components(self.operations):
            for op in comp:
                print op
            print


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
            if not op[0] in translate.keys():
                print op
            op[0] = translate[ op[0] ]
            op[-1] = translate[ op[-1] ]



### the following functions are related to simplifying the differential
    def find_differentials(self):
        '''returns a list of all operations with no Reeb chords.
        assumes all boundary components are type D'''
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

        #should only be used when all boundaries are type D (??)'''

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
                            new_op.append( op1[i+1] + op2[i+1] )
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

##          Note: this is the first difference from the tree version... not having to deal with generator_signs
            #remove initial_generator and final_generator
            self.generators.remove(initial_generator)
            self.generators.remove(final_generator)
            self.reduce_mod_2()           

            differential = self.find_removable_differential()
            

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


    

    def hochschild(self, i, j):
        '''performs Hochschild homology, gluing the i boundary (type A)
        to the j boundary (type D)'''

##        if not self.is_bounded():
##            print 'WARNING: Performing Hochschild homology on unbounded module'

        #first we select which generators survive the hochschild homology: those whose idempotents for the i and j boundaries match
        new_generators = []
        for gen in self.generators:
            if self.idempotence[gen][i] == self.idempotence[gen][j]:
                new_generators.append(gen)

        #we compute the new operations by finding sequences of operations such that
                #   -the first operation has no type A inputs for the ith boundary
                #   -the last operation has no type D outputs for the jth boundary
                #   -for each n, the list of i-boundary inputs of the first n operations is the beginning
                #       of the list of j-boundary outputs of the first n operations
                
        new_op_lists = []           #list: each element is a list of operations produce a new operation in hochschild homology
        for gen1 in new_generators:
            queue = []              #list: each element is a list of operations, which may be the first part of a list that belongs in new_op_lists
            for op in self.operations_starting_with(gen1):
                if op[i+1] == []:                   #since this is the first operation, we only consider operations with no A inputs for the ith boundary
                    if op[j+1] == '':               #if the op also has no outputs, then this single operation will produce a new operation
                        new_op_lists.append([op])
                    else:                           #otherwise, it may be the start of a list of operations that produces an operation, so we add the partial list to queue
                        queue.append([op])
            counter = 0
            while len(queue) > 0 and counter == 0:
                list_of_ops = queue.pop(0)      #take a partial list of ops from the queue to try to extend

                if len(list_of_ops) > 30:       #prevents getting stuck in infinite loops and gives an indicating message in the event that the module is not relatively bounded in the appropriate way
                    print 'loop timeout'
                    counter = 1

                A_sequence = list_sum([ op[i+1] for op in list_of_ops])         #the list of i-boundary type A inputs so far
                D_sequence = list_sum([ [op[j+1]] for op in list_of_ops])       #the list of j-boundary type D outputs so far
                #at this point, A_sequence should be the first part of D_sequence
                for n in range(len(A_sequence)):
                    if not A_sequence[n] == D_sequence[n]:
                        print "ERROR: A, D sequence mismatch"       #for debugging. this error should not occur, since list_of_ops was taken from queue

                D_sequence = D_sequence[len(A_sequence):]           #the j-boundary type D outputs that are still unused as type A inputs
                gen2 = list_of_ops[-1][-1]                          #the generator the partial list of opertions ends on
                for op in self.operations_starting_with(gen2):
                    if op[i+1] == D_sequence and op[j+1] == '':         #if an op starting at gen2 has no j-boundary outputs and uses up all the previous outputs as i-boundary inputs, then adding this op completes the partial list to a list producing a new operation
                        new_op_lists.append( list_of_ops + [op] )
                    if op[j+1] != '' and op[i+1] == D_sequence[:len(op[i+1])]:      #if there is a j-boundary output and the i-boundary inputs use the first (some number of) terms of the previous outputs, we add this operation and put back in the queue
                        queue.append( list_of_ops + [op] )

            
        #now every new operation is represented by a list of operations in new_op_lists
        #we convert these to new opeartaions
        new_operations = []
        for list_of_ops in new_op_lists:
            new_op = [list_of_ops[0][0]]                    #initial generator of new op
            for n in range(len(self.boundary_types)):
                if n == i or n == j:                        #we skip the ith and jth boundaries, since these no longer exist after gluing
                    pass
                elif self.boundary_types[n] == 'A':                                 #for type A boundaries (other than i), we concatenate the inputs for for the ops in list_of_ops
                    new_op.append( list_sum([ op[n+1] for op in list_of_ops ]) )
                elif self.boundary_types[n] == 'D':                                 #for type D boundaries (other than j), we multiply the outputs for the ops in list_of_ops 
                    new_op.append( torus_product_list( [op[n+1] for op in list_of_ops ]) )
            new_op.append(list_of_ops[-1][-1])              #final generator of new op
            if not None in new_op:                          #None will be in new_op if one of the type D products is zero. If not, we add this new operation
                new_operations.append(new_op)
                
        self.operations = new_operations            #assign the updated operations to the multimodule
        self.generators = new_generators            #assign the updated generators to the multimodule
                
            
        new_idempotence = {}                        #we remove the ith and jth boundary entry from the idempotents list for each generator
        for gen in new_generators:
            new_i = self.idempotence[gen]
            del new_i[max(i, j)]
            del new_i[min(i, j)]
            new_idempotence[gen] = new_i
        self.idempotence = new_idempotence          

        del self.fiber_directions[max(i,j)]         #we remove the ith and jth boundary entry from the list of fiber directions
        del self.fiber_directions[min(i,j)]         

        del self.boundary_types[max(i,j)]           #we remove the ith and jth boundary entry from the list of boundary types
        del self.boundary_types[min(i,j)]

        self.reduce_mod_2()
        self.remove_differentials()
        

    def is_there_a_cycle(self, gen, ancestors):
        '''a recursive function used in is_bounded
        will return True if gen is '''
        if gen in ancestors:
            return True
        for gen2 in [op[-1] for op in self.operations_starting_with(gen)]:
            if self.is_there_a_cycle(gen2, ancestors+[gen]):
                return True
        return False

    def is_bounded(self):
        '''considering the set of operations of self as an (unlabeled) directed graph
        this function returns True if there are no cycles. this means that the module is bounded'''
        for gen in self.generators:
            if self.is_there_a_cycle(gen, []):
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

    def copy(self):
        '''returns a new multimodule instance with the same
        generators, operations, boundary_types, idempotence, fiber_directions, and name'''
        result = Multimodule(self.generators[:], self.boundary_types[:], self.fiber_directions[:],dict(self.idempotence), [op[:] for op in self.operations[:]])
        result.name = self.name
        return result

    def self_glue(self, i, j):
        '''i < j should be indices for two type D boundaries
        we add, without canceling any differentials, CFAA_id and self_gluer to the j boundary
        and CFAA_id to the i boundary. This ensures boundedness.
        
        Returns the result of hochschild homology with respect to
        the new type A i-boundary and type D j-boundary.'''
        n = len(self.boundary_types)
        result = tensor(CFAA_id, 0, self, j, False)
        result = tensor(result, 0, self_gluer, 0, False)
        result = tensor(CFAA_id, 0, result, i - (i>j), False)

        result.relabel_generators()
        result.hochschild(0,n-1)

        return result


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

    if not (module1.boundary_types[boundary1] == 'A' and module2.boundary_types[boundary2] == 'D'):
        print 'ERROR: inputs to tensor() must be a type A and a type D boundary, respectively'

    #determine new generators
    new_gens = []
    for i in range(len(module1.generators)):
        gen1 = module1.generators[i]
        for j in range(len(module2.generators)):
            gen2 = module2.generators[j]
            if module1.idempotence[gen1][boundary1] == module2.idempotence[gen2][boundary2]:
                new_gens.append( (gen1, gen2) )
                



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

    if simplify:
        new_module.remove_differentials()
    #new_module.relabel_generators()

    new_module.name = (module1.name, boundary1, module2.name, boundary2)

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
CFAA_id.name = 'CFAA_id'
###########

####### CFDD_id  ########
generators = ['a', 'b']
boundary_types = ['D', 'D']
fiber_directions = [2, 1]
idempotence = {'a': [2, 2], 'b': [1, 1]}
operations = [['a', '2', '2', 'b'],
              ['b', '1', '3', 'a'],
              ['b', '123', '123', 'a'],
              ['b', '3', '1', 'a']]
CFDD_id = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFDD_id.name = 'CFDD_id'
###########

########   Self-gluer   ###########

generators = ['afi', 'afj', 'afk', 'afl', 'ahi',
              'ahj', 'ahk', 'ahl', 'ang', 'amg',
              'ebi', 'ebj', 'ebk', 'ebl', 'edi',
              'edj', 'edk', 'edl', 'enc', 'emc']
boundary_types = ['D','D']
fiber_directions = [1,2]
idempotence = {'afi':[2, 2],
               'afj':[2, 1],
               'afk':[2, 2],
               'afl':[2, 1],
               'ahi':[2, 2],
               'ahj':[2, 1],
               'ahk':[2, 2],
               'ahl':[2, 1],
               'ang':[2, 2],
               'amg':[2, 2],
               'ebi':[1, 2],
               'ebj':[1, 1],
               'ebk':[1, 2],
               'ebl':[1, 1],
               'edi':[1, 2],
               'edj':[1, 1],
               'edk':[1, 2],
               'edl':[1, 1],
               'enc':[1, 1],
               'emc':[1, 1]}
operations = [ ['afi', '2', '', 'edi'],
               ['afj', '', '1', 'afi'],
               ['afj', '2', '', 'edj'],
               ['afk', '', '2', 'afj'],
               ['afk', '2', '', 'edk'],
               ['afl', '', '123', 'afi'],
               ['afl', '', '3', 'afk'],
               ['afl', '', '1', 'ang'],
               ['afl', '2', '', 'edl'],

               ['ahj', '', '1', 'ahi'],
               ['ahk', '', '2', 'ahj'],
               ['ahl', '', '123', 'ahi'],
               ['ahl', '', '3', 'ahk'],

               ['ang', '2', '2', 'enc'],
               ['amg', '', '23', 'afi'],
               ['amg', '', '', 'ahi'],
               ['amg', '23', '', 'ahk'],
               ['amg', '2', '2', 'emc'],
               ['amg', '2', '2', 'ebl'],

               ['ebi', '3', '', 'afi'],
               ['ebj', '3', '', 'afj'],
               ['ebj', '', '1', 'ebi'],
               ['ebj', '12', '', 'enc'],
               ['ebk', '3', '', 'afk'],
               ['ebk', '1', '', 'ang'],
               ['ebk', '', '2', 'ebj'],
               ['ebl', '3', '', 'afl'],
               ['ebl', '', '123', 'ebi'],
               ['ebl', '', '3', 'ebk'],
               ['ebl', '', '', 'enc'],

               ['edi', '1', '', 'ahi'],
               ['edj', '1', '', 'ahj'],
               ['edj', '', '1', 'edi'],
               ['edk', '1', '', 'ahk'],
               ['edk', '', '2', 'edj'],
               ['edl', '1', '', 'ahl'],
               ['edl', '', '123', 'edi'],
               ['edl', '', '3', 'edk'],
               ['edl', '', '12', 'enc'],

               ['enc', '1', '3', 'ang'],
               ['enc', '3', '1', 'ang'],
               ['enc', '123', '123', 'ang'],
               ['emc', '1', '3', 'amg'],
               ['emc', '3', '1', 'amg'],
               ['emc', '123', '123', 'amg'],
               ['emc', '3', '', 'ahj'],
               ['emc', '', '3', 'edi'],         #end of definite operations

               #['amg', '', '', 'ang'],
               #['emc', '', '', 'enc'],


               ['ebi', '123', '', 'ang'],       #12389,10,11
               ['emc', '123', '', 'ahl'],       #123789,10

               ['emc', '', '123', 'ebi'],       #456789,10
               ['ahl', '', '123', 'ang'] ]       #45689,10,11

               #['ahi', '2', '2', 'enc'],        #25899,10,11
               #['amg', '2', '2', 'ebl'] ]        #257899,10
               #]

self_gluer = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
self_gluer.remove_differentials()
self_gluer.name = 'self_gluer'
###########

######### CFD( solid torus )    meridian alpha arc is alpha_2
generators = ['p']
boundary_types = ['D']
fiber_directions = [1]
idempotence = {'p':[2]}
operations = [ ['p', '23', 'p'] ]

CFD_2 = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFD_2.name = 'CFD_2'

######### CFD( solid torus )    meridian alpha arc is alpha_1

generators = ['p']
boundary_types = ['D']
fiber_directions = [2]
idempotence = {'p':[1]}
operations = [ ['p', '12', 'p'] ]

CFD_1 = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFD_1.name = 'CFD_1'


CFA_2 = tensor(CFAA_id,0, CFD_2, 0)
CFA_1 = tensor(CFAA_id,1, CFD_1, 0)
# Note: tensoring CFA_1 with CFD_1  or CFA_2 with CFD_2 gives S^3
#       tensoring CFA_1 with CFD_2 gives S^2 x S^1

########### CFDDD( pair-of-pants X S^1 )

generators = ['v', 'w', 'x', 'y', 'z']
boundary_types = ['D','D','D']
fiber_directions = [2,2,1]
idempotence = {'v': [1, 1, 1],
               'w': [1, 1, 1],
               'x': [2, 1, 1],
               'y': [2, 2, 2],
               'z': [1, 2, 1]}
operations = [ ['v','1'  ,'123','3'  , 'y'],
               ['v','123','123','123', 'y'],
               ['v', ''  , '3' , ''  , 'z'],
               ['v', '3' , ''  , ''  , 'x'],
               
               ['w', '1' , '1' , '3' , 'y'],
               ['w', '3' , ''  ,'12' , 'x'],
               ['w','123', '1' ,'123', 'y'],
               
               ['x', '2' , ''  ,'12' , 'v'],
               ['x', ''  , '3' , '1' , 'y'],
               ['x', '2' , ''  , ''  , 'w'],
               
               ['y', ''  , '2' , '2' , 'x'],
               ['y', '2' , ''  , '2' , 'z'],
               
               ['z', ''  , '2' , ''  , 'w'],
               ['z', '3' , ''  , '1' , 'y']   ]

CFDDD = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFDDD.name = 'CFDDD'


#mirror CFDDD
boundary_types = ['D','D','D']
fiber_directions = [1,1,2]
idempotence = {'v': [2, 2, 2],
               'w': [2, 2, 2],
               'x': [1, 2, 2],
               'y': [1, 1, 1],
               'z': [2, 1, 2]}
operations = [ ['y','3'  ,'123','1'  , 'v'],
               ['y','123','123','123', 'v'],
               ['z', ''  , '1' , ''  , 'v'],
               ['x', '1' , ''  , ''  , 'v'],
               
               ['y', '3' , '3' , '1' , 'w'],
               ['x', '1' , ''  ,'23' , 'w'],
               ['y','123', '3' ,'123', 'w'],
               
               ['v', '2' , ''  ,'23' , 'x'],
               ['y', ''  , '1' , '3' , 'x'],
               ['w', '2' , ''  , ''  , 'x'],
               
               ['x', ''  , '2' , '2' , 'y'],
               ['z', '2' , ''  , '2' , 'y'],
               
               ['w', ''  , '2' , ''  , 'z'],
               ['y', '1' , ''  , '3' , 'z']   ]

CFDDD_mirror = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFDDD_mirror.name = 'CFDDD_mirror'

##CFAD_1_2 = tensor(CFAA_id, 1, CFDDD, 0)
##CFAD_1_2 = tensor(CFA_1, 0, CFAD_1_2, 2)
##
##CFAD_2_1 = tensor(CFAA_id, 0, CFDDD_mirror, 0)
##CFAD_2_1 = tensor(CFA_2, 0, CFAD_2_1, 2)

CFDD_2 = tensor(CFA_1, 0, CFDDD, 2)         #a type D bimodule used to switch boundary parametrization. The fiber is alpha_2 on both sides
CFDD_2.name = 'CFDD_2'
CFDD_1 = tensor(CFA_2, 0, CFDDD_mirror, 2)  #a type D bimodule with fiber alpha_1 on both sides
CFDD_1.name = 'CFDD_1'


######### CFD( fig 8 complement )    meridian alpha arc is alpha_1

generators = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
boundary_types = ['D']
fiber_directions = [2]
idempotence = {'a':[1], 'b':[2], 'c':[1], 'd':[2], 'e':[2], 'f':[1], 'g':[2], 'h':[1], 'i':[1]}
operations = [ ['a', '1', 'd'],
               ['b', '2', 'a'],
               ['c', '3', 'b'],
               ['c', '1', 'e'],
               ['f', '123', 'd'],
               ['g', '2', 'f'],
               ['h', '3', 'g'],
               ['h', '123', 'e'],
               ['i', '12', 'i'] ]

CFD_fig8 = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
CFD_fig8.name = 'CFD_fig8'

#we can change the parametrization so that the meridian is alpha_2 by tensoring with CFDD_2
CFD_fig8_flip = tensor(CFAA_id, 1, CFDD_1, 1)
CFD_fig8_flip = tensor(CFD_fig8_flip, 0, CFD_fig8, 0)
CFD_fig8_flip.name = 'CFD_fig8_flip'


##########  Dehn Twists  ##########
##########               #########
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
CFDA_twist1.name = 'CFDA_twist1'

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
CFDA_twist1_inv.name = 'CFDA_twist1_inv'



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
CFDA_twist2.name = 'CFDA_twist2'

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
CFDA_twist2_inv.name = 'CFDA_twist2_inv'

########  Plus_genus  ########
#plus_genus is a DD bimodule corresponding to the trivial S^1 bundle over a genus 1 surface with two boundary components
#fiber directions are [2,1]
plus_genus = CFDDD_mirror.self_glue(1,2)
plus_genus = tensor(CFAA_id, 0, plus_genus, 0)
plus_genus = tensor(plus_genus, 0, CFDDD, 1)
plus_genus.name = 'plus_genus'


########  minus_genus  ########
#minus_genus is a DD bimodule which adds -1 to the genus of a vertex corresponding to a given multimodule
#according to Neumann's calculus for plumbing, subtracting 1 from the genus of a vertex is equivalentl to glueing on the tree
#          (-2)
#         /  
# --- (-2)
#         \
#          (-2)
# to make a bimodule, we splice this tree (with one boundary component) to a copy of CFDDD
# the resulting bimodule can be glued to a vertex_module to decrease the genus
# the fiber directions are [2,1]
#
# the following form of this module was computed using this graph manifolds program
generators = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
boundary_types = ['D','D']
fiber_directions = [2,1]
idempotence = {'a': [2, 1], 'c': [2, 2], 'b': [2, 2], 'e': [2, 1], 'd': [2, 1], 'g': [1, 1], 'f': [1, 1], 'h': [2, 1]}
operations = [['a', '', '1', 'c'], ['a', '23', '1', 'b'], ['b', '', '2', 'a'], ['c', '23', '2', 'a'], ['d', '', '12', 'e'], ['d', '23', '', 'e'], ['e', '', '12', 'd'], ['e', '23', '', 'd'], ['f', '1', '123', 'c'], ['f', '1', '3', 'b'], ['f', '123', '123', 'b'], ['f', '123', '3', 'c'], ['f', '3', '', 'h'], ['g', '1', '3', 'c'], ['g', '123', '123', 'c'], ['g', '3', '12', 'h'], ['h', '2', '', 'g'], ['h', '2', '12', 'f']]
minus_genus = Multimodule(generators, boundary_types, fiber_directions, idempotence, operations)
minus_genus.name = 'minus_genus'






##################
###############  End of saved multimodules
##################
##################
##################
###############  Functions for plumbing trees
##################



# Objects in the following class, WeightedTree, are weighted acyclic graphs.
class BiweightedGraph:
    def __init__(self, weights1, weights2, edges, edge_signs):
        '''weights1 is list of length #vertices, specifying genus
        weight2 is list of length #vertices, specifying euler number
        edges is list of pairs (i, j), representing an edge btwn vertices i and j
        
        an edge may have the form (i, '*'), indicating that there is an open
        boundary component at the vertex i
        '''

        if len(weights1) != len(weights2):
            print "error, weights mismatch"
            return
        self.num_vertices = len(weights1)
        self.vertices = {}
        for i in range(len(weights1)):
            self.vertices[i] = (weights1[i], weights2[i])
        self.edges = edges
        self.edge_signs = edge_signs

        #edges_from_vertex[i] gives an (ordered) list of vertices (represented by integers) connected to v_i
        self.edges_from_vertex = []
        for i in range(len(weights1)):
            edges_from_i = []
            for j in range(len(edges)):
                if edges[j][0] == i:
                    edges_from_i.append(j)
                if edges[j][1] == i:
                    edges_from_i.append(j)
            self.edges_from_vertex.append(edges_from_i)
                

def compute_vertex_module(genus, euler_number, edges):
    '''constructs the type D multimodule for the S^1 bundle of given euler number
    over the surface of given genus, with (#edges) boundary components,
    labeled by edges.
    fiber direction is 2 on all boundaries except for the last if there is more than one

    genus: integer >= 0
    euler_number: integer
    edges: list of integers'''

    ## produce trivial bundle on disk with right number of boundary components
    if len(edges) == 1 or len(edges) == 0:
        base = CFD_1.copy()
    elif len(edges) == 2:
        base = CFDD_id.copy()
    else:
        base = CFDDD.copy()
        while len(base.boundary_types) < len(edges):
            base = tensor(CFAA_id, 0, base, len(base.boundary_types)-1)
            base = tensor(base, 0, CFDDD, 0)
    #at this point, base has the right boundary components.
    #for >1 boundary, the last one is fiber 1, all others are fiber 2
    #in the case of one boundary component, it is fiber 2


    ## add genus ##
    counter = 0
    while counter < genus:
        counter += 1
        x = tensor(CFAA_id, 0, plus_genus, 1)
        base = tensor(x, 0, base, 0)

    ## subtract genus ##
    counter = 0
    while counter > genus:
        counter -= 1
        x = tensor(CFAA_id, 0, minus_genus, 1)
        base = tensor(x, 0, base, 0)

        

    ## add twist
    twist = CFDA_twist2
    twist_inv = CFDA_twist2_inv

    if euler_number > 0:
        for i in range(euler_number):
            base = tensor(twist, 1, base, 0)
    if euler_number < 0:
        for i in range(-euler_number):
            base = tensor(twist_inv, 1, base, 0)

    if len(edges) == 0:
        base = tensor(CFA_2, 0, base, 0)

    base.remove_differentials()

    #base.relabel_generators()
    
    base.boundary_labels = edges    #for each boundary of base, there is a corresponding integer representing an edge

    return base


def plumb_edge(module1, module2, edge_label, edge_sign):
    '''module1 and module2 are type D multimodules
    each one has a boundary component labeled edge_label
    returns single multimodule obtained from gluing these two boundary
    components, with CFAA_id in between
    edge_sign is +1 or -1'''

    ###this is a quick fix, because one sign convention is wrong (edge signs or +/- dehn twists)
    ###switching either one should fix it, so I'm choosing for now to flip the edge_sign convention
    edge_sign = -edge_sign
    ###
    
    i = module1.boundary_labels.index(edge_label)
    j = module2.boundary_labels.index(edge_label)
    fiber1 = module1.fiber_directions[i]
    fiber2 = module2.fiber_directions[j]

    n = len(module1.boundary_types)     #number of bdy cpnts in module1

    if (fiber1, fiber2, edge_sign) == (1,1,-1):
        # 1:1
        result = tensor(CFAA_id, 0, module1, i)
        result = tensor(result, 0, module2, j)
    elif (fiber1, fiber2, edge_sign) == (1,1,1):
        # 1:2-2:2-2:1           (colon means gluing tow modules, dash connects boundaries of bimodule. A 1:1 gluing is negative edge, a 2-2 gluing is a positive edge, and 1:2 or 2:1 is plumbing
        result = tensor(CFAA_id, 0, module1, i)
        result = tensor(result, 0, CFDD_2, 0)   #currently: mod1 : AA : DD_2, relevant bdy is last
        result = tensor(CFAA_id, 1, result, n-1)
        result = tensor(result, 0, CFDD_2, 0)   #currently: mod1::DD_2::DD_2
        result = tensor(CFAA_id, 0, result, n-1)
        result = tensor(result, 0, module2, j)       
    elif (fiber1, fiber2, edge_sign) == (1,2,-1):
        # 1:1-1:2
        result = tensor(CFAA_id, 0, module1, i)
        result = tensor(result, 0, CFDD_1, 0)
        result = tensor(CFAA_id, 0, result, n-1)
        result = tensor(result, 0, module2, j)
    elif (fiber1, fiber2, edge_sign) == (1,2,1):
        # 1:2-2:2
        result = tensor(CFAA_id, 0, module1, i)
        result = tensor(result, 0, CFDD_2, 0)
        result = tensor(CFAA_id, 1, result, n-1)
        result = tensor(result, 0, module2, j)
    elif (fiber1, fiber2, edge_sign) == (2,1,-1):
        # 2:1-1:1
        result = tensor(CFAA_id, 0, module1, i)
        result = tensor(result, 0, CFDD_1, 0)
        result = tensor(CFAA_id, 1, result, n-1)
        result = tensor(result, 0, module2, j)        
    elif (fiber1, fiber2, edge_sign) == (2,1,1):
        #2:2-2:1
        result = tensor(CFAA_id, 0, module1, i)
        result = tensor(result, 0, CFDD_2, 0)
        result = tensor(CFAA_id, 1, result, n-1)
        result = tensor(result, 0, module2, j)
    elif (fiber1, fiber2, edge_sign) == (2,2,-1):
        #2:1-1:1-1:2
        result = tensor(CFAA_id, 0, module1, i)
        result = tensor(result, 0, CFDD_1, 0)
        result = tensor(CFAA_id, 1, result, n-1)
        result = tensor(result, 0, CFDD_1, 0)
        result = tensor(CFAA_id, 0, result, n-1)
        result = tensor(result, 0, module2, j) 
    elif (fiber1, fiber2, edge_sign) == (2,2,1):
        #2:2
        result = tensor(CFAA_id, 0, module1, i)
        result = tensor(result, 0, module2, j)
    else:
        print 'plumb_edge ERROR', fiber1, fiber2, edge_sign

    result.boundary_labels = ( module1.boundary_labels[:i] +
                               module1.boundary_labels[i+1:] +
                               module2.boundary_labels[:j] +
                               module2.boundary_labels[j+1:] )
    return result


def plumb_self_edge(module, edge_label, edge_sign):
    '''module is type D multimodule, which has two boundary components labeled by edge_label
    returns the multimodule obtained from self-gluing along these two boundary
    components, using the multimodule self_glue function
    edge_sign is +1 or -1'''

    ###this is a quick fix, because one sign convention is wrong (edge signs or +/- dehn twists)
    ###switching either one should fix it, so I'm choosing for now to flip the edge_sign convention
    edge_sign = -edge_sign
    ###

    boundary_indices = []
    for i in range(len(module.boundary_labels)):
        if module.boundary_labels[i] == edge_label:
            boundary_indices.append(i)
    i = boundary_indices[0]
    j = boundary_indices[1]     #i and j are the indices of the two boundaries to be glued
    fiber1 = module.fiber_directions[i]
    fiber2 = module.fiber_directions[j]

    n = len(module.boundary_types)     #number of bdy cpnts in module

    if (fiber1, fiber2, edge_sign) == (1,1,-1):
        # 1:1
        result = module
        result = result.self_glue(i, j)
    elif (fiber1, fiber2, edge_sign) == (1,1,1):
        # 1:2-2:2-2:1
        result = tensor(CFAA_id, 0, module, i)
        result = tensor(result, 0, CFDD_2, 0)   #currently: mod1 : AA : DD_2, relevant bdy is last
        result = tensor(CFAA_id, 1, result, n-1)
        result = tensor(result, 0, CFDD_2, 0)   #currently: mod1::DD_2::DD_2
        result = result.self_glue(j-1, n-1)      
    elif (fiber1, fiber2, edge_sign) == (1,2,-1):
        # 1:1-1:2
        result = tensor(CFAA_id, 0, module, i)
        result = tensor(result, 0, CFDD_1, 0)
        result = result.self_glue(j-1, n-1)
    elif (fiber1, fiber2, edge_sign) == (1,2,1):
        # 1:2-2:2
        result = tensor(CFAA_id, 0, module, i)
        result = tensor(result, 0, CFDD_2, 0)
        result = result.self_glue(j-1, n-1)
    elif (fiber1, fiber2, edge_sign) == (2,1,-1):
        # 2:1-1:1
        result = tensor(CFAA_id, 0, module, i)
        result = tensor(result, 0, CFDD_1, 0)
        result = result.self_glue(j-1, n-1)      
    elif (fiber1, fiber2, edge_sign) == (2,1,1):
        #2:2-2:1
        result = tensor(CFAA_id, 0, module, i)
        result = tensor(result, 0, CFDD_2, 0)
        result = result.self_glue(j-1, n-1)
    elif (fiber1, fiber2, edge_sign) == (2,2,-1):
        #2:1-1:1-1:2
        result = tensor(CFAA_id, 0, module, i)
        result = tensor(result, 0, CFDD_1, 0)
        result = tensor(CFAA_id, 1, result, n-1)
        result = tensor(result, 0, CFDD_1, 0)
        result = result.self_glue(j-1, n-1)
    elif (fiber1, fiber2, edge_sign) == (2,2,1):
        #2:2
        result = module1
        result = result.self_glue(i, j)

    result.boundary_labels = ( module.boundary_labels[:i] +
                               module.boundary_labels[i+1:j] +
                               module.boundary_labels[j+1:] )
    return result


def plumb(graph):
    '''given a graph (instance of BiWeightedGraph), computes HF hat of the corresponding plumbed 3 manifold
    returns a multimodule. if the graph represents a closed 3 manifold, this multimodule will be CF-hat.
    If there open boundary components (edges with at '*'), then the result will be the appropriate type D bordered invariant'''

    #we first compute a module corresponding to each vertex, and store them in a list
    vertex_modules = [compute_vertex_module(graph.vertices[i][0], graph.vertices[i][1], graph.edges_from_vertex[i]) for i in range(graph.num_vertices)]
    
    #next we plumb these modules together according to the edges until we have one connected multimodule
    while len(vertex_modules) > 1:
        vertex_module1 = vertex_modules.pop(0)

        i = 0                                                   #
        edge = vertex_module1.boundary_labels[i]        #edge is an integer, giving the edge in graph.edges which corresponds to the first boundary of vertex_module1
                                                        #which glues to another vertex (the corresponding edge does not conatin '*')
        while '*' in graph.edges[edge]:                         # added to allow for '*' edges
            i += 1                                              # ensures edge (int) is a non-boundary edge
            edge = vertex_module1.boundary_labels[i]            # 

        if edge in vertex_module1.boundary_labels[i+1:]:        # if the edge we have chosen is a self edge (it corresponds to two boundary components on this vertex module), then we self_plumb
            vertex_modules.append( plumb_self_edge(vertex_module1, edge, graph.edge_signs[edge]) )
        else:                                                   # if the edge we have chosen is not a self edge, we find the other vertex_module it attaches to and plumb the two modules
            other_module_index = None
            for i in range(len(vertex_modules)):
                if edge in vertex_modules[i].boundary_labels:
                    other_module_index = i
            vertex_module2 = vertex_modules.pop(other_module_index)
            vertex_modules.append(plumb_edge(vertex_module1, vertex_module2, edge, graph.edge_signs[edge]))

    #at this point, vertex_modules is a list with a single element
    vertex_module = vertex_modules[0]
    num_boundary_components = sum([ ('*' in edge) for edge in graph.edges])
    while len(vertex_module.boundary_labels) > num_boundary_components:
        i = 0
        edge = vertex_module.boundary_labels[i]
        while '*' in graph.edges[edge]:
            i += 1
            edge = vertex_module.boundary_labels[i]
        #now edge is an edge conencted to vertex_module with no '*'. Thus it must be a self edge
        vertex_module = plumb_self_edge(vertex_module, edge, graph.edge_signs[edge])
        
    return vertex_module

##genus = [0,0,0,0]
##weights = [-1,-2,-3,-6]
##edges = [(0,1),(0,2),(0,3)]
##edge_signs = [1,1,1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [0,0,0,0]
##weights = [-1,-2,-4,-4]
##edges = [(0,1),(0,2),(0,3)]
##edge_signs = [1,1,1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [0,0,0,0]
##weights = [-1,-3,-3,-3]
##edges = [(0,1),(0,2),(0,3)]
##edge_signs = [1,1,1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [0,0,0,0,0,0,0]
##weights = [-2,-2,-2,-2,-2,-2,-2]
##edges = [(0,1),(1,2),(0,3),(3,4),(0,5),(5,6)]
##edge_signs = [1,1,1,1,1,1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [0,0,0,0,0,0,0,0]
##weights = [-2,-2,-2,-2,-2,-2,-2,-2]
##edges = [(0,1),(0,2),(2,3),(3,4),(0,5),(5,6),(6,7)]
##edge_signs = [1,1,1,1,1,1,1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [0,0,0,0,0,0,0,0,0]
##weights = [-2,-2,-2,-2,-2,-2,-2,-2,-2]
##edges = [(0,1),(0,2),(2,3),(0,4),(4,5),(5,6),(6,7),(7,8)]
##edge_signs = [1,1,1,1,1,1,1,1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)


##genus = [0]
##weights = [-1]
##edges = [(0,0)]
##edge_signs = [1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##
##genus = [0]
##weights = [0]
##edges = [(0,0)]
##edge_signs = [1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [0]
##weights = [1]
##edges = [(0,0)]
##edge_signs = [1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [0]
##weights = [-1]
##edges = [(0,0)]
##edge_signs = [-1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [0]
##weights = [0]
##edges = [(0,0)]
##edge_signs = [-1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [0]
##weights = [1]
##edges = [(0,0)]
##edge_signs = [-1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)



##genus = [1,2,0,0,0,0,0,0,0,0]
##weights = [3,-2,4,5,-1,-2,-2,-2,-2,-2]
##edges = [(0,1),(0,2),(0,3),(3,4),(4,5),(4,6),(3,7),(7,8),(8,9)]
##edge_signs = [-1,-1,-1,-1,-1,-1,-1,-1,-1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)


##genus = [1, 0, -1]
##weights = [0,0,0]
##edges = [(0,1),(1,2)]
##edge_signs = [-1,-1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)


##genus = [1, 0, 0, 0, 0, 0]
##weights = [0,0,0,0,2,-2]
##edges = [(0,1),(1,2),(2,3),(3,4),(3,5)]
##edge_signs = [1,-1,-1,-1,-1]
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)
##
##genus = [-1]
##weights = [0]
##edges = []
##edge_signs = []
##test = BiweightedGraph(genus, weights, edges, edge_signs)
##x = plumb(test)
##print len(x.generators)

##genus = [0,0,0,0]
##weights = [0,2,-2,0]
##edges = [(0,1),(0,2),(0,3),(3,'*'),(3,'*')]
##edge_signs = [-1,-1,-1,-1,-1]
##graph = BiweightedGraph(genus, weights, edges, edge_signs)
##result = plumb(graph)
