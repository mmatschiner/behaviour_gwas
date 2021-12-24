# m_matschiner Tue Mar 10 10:51:36 CET 2020

# The Tree class.
class Tree(object):

    def __init__(self, newick_string):
        self.newick_string = newick_string
        self.edges = []

    def get_newick_string(self):
        return self.newick_string

    def parse_newick_string(self):
        working_newick_string = self.newick_string
        number_of_internal_nodes = 0
        number_of_edges = 0

        # Remove comments from the tree string.
        pattern = re.compile("\[.*?\]")
        hit = "placeholder"
        while hit != None:
            hit = pattern.search(working_newick_string)
            if hit != None:
                working_newick_string = working_newick_string.replace(hit.group(0),"")

        # Check whether a branch above the root is present, and if so, remove it.
        if working_newick_string[0:2] == "((" and working_newick_string[-1] == ")" and working_newick_string[-2] != ")":
            level = 0
            newick_string_tail_start_pos = 0
            newick_string_tail = ""
            for pos in range(len(working_newick_string)):
                if working_newick_string[pos] == "(":
                    level += 1
                if working_newick_string[pos] == ")":
                    level -= 1
                if level == 1 and pos > 1:
                    newick_string_tail_start_pos = pos
                    newick_string_tail = working_newick_string[pos+1:]
                    break
            if newick_string_tail.count(",") == 0:
                working_newick_string = working_newick_string[1:newick_string_tail_start_pos+1]

        # Parse the bifurcating part of the tree.
        if ":" in working_newick_string:
            pattern = re.compile("\(([a-zA-Z0-9_\.\-]+?):([\d\.Ee-]+?),([a-zA-Z0-9_\.\-]+?):([\d\.Ee-]+?)\)")
        else:
            print("ERROR: It appears that the tree string does not include branch lengths!")
            sys.exit(1)
        hit = "placeholder"
        while hit != None:
            hit = pattern.search(working_newick_string)
            if hit != None:
                number_of_internal_nodes += 1
                number_of_edges += 2
                internal_node_id = "internalNode" + str(number_of_internal_nodes) + "X"
                edge1 = Edge("edge" + str(number_of_edges-1) + "X")
                edge1.set_node_ids([internal_node_id, hit.group(1)])
                edge1.set_length(float(hit.group(2)))
                edge2 = Edge("edge" + str(number_of_edges) + "X")
                edge2.set_node_ids([internal_node_id, hit.group(3)])
                edge2.set_length(float(hit.group(4)))
                self.edges.append(edge1)
                self.edges.append(edge2)
                working_newick_string = working_newick_string.replace(hit.group(0), internal_node_id)

        # Make sure the remaining string includes a single node and use this node id to determine root edges.
        pattern_rooted = re.compile("^internalNode\d+X$")
        hit_rooted = pattern_rooted.search(working_newick_string)
        if hit_rooted == None:
            print('ERROR: The newick tree string could not be parsed!')
            print('The remaining unparsed part of the newick string: \"{}\"'.format(working_newick_string))
            sys.exit(1)
        else:
            root_node_id = hit_rooted.group(0)
            for edge in self.get_edges():
                if edge.get_node_ids()[0] == root_node_id:
                    edge.set_is_root_edge(True)
                else:
                    edge.set_is_root_edge(False)

    def parse_extended_newick_string(self):
        working_newick_string = self.newick_string
        number_of_internal_nodes = 0
        number_of_edges = 0

        # Check whether a branch above the root is present, and if so, remove it.
        if working_newick_string[0:2] == "((" and working_newick_string[-1] == ")" and working_newick_string[-2] != ")":
            level = 0
            newick_string_tail_start_pos = 0
            newick_string_tail = ""
            for pos in range(len(working_newick_string)):
                if working_newick_string[pos] == "(":
                    level += 1
                if working_newick_string[pos] == ")":
                    level -= 1
                if level == 1 and pos > 1:
                    newick_string_tail_start_pos = pos
                    newick_string_tail = working_newick_string[pos+1:]
                    break
            if newick_string_tail.count(",") == 0:
                working_newick_string = working_newick_string[1:newick_string_tail_start_pos+1]

        # Parse the bifurcating part of the tree.
        if ":" in working_newick_string:
            if "[" in working_newick_string:
                pattern = re.compile("\(([a-zA-Z0-9_\.\-]+?)\[\&dmv=\{([\d\.]+?)\}\]:([\d\.Ee-]+?),([a-zA-Z0-9_\.\-]+?)\[\&dmv=\{([\d\.]+?)\}\]:([\d\.Ee-]+?)\)")
            else:
                print("ERROR: It appears that the tree string does not include annotation!")
                sys.exit(1)
        else:
            print("ERROR: It appears that the tree string does not include branch lengths!")
            sys.exit(1)
        hit = "placeholder"
        while hit != None:
            hit = pattern.search(working_newick_string)
            if hit != None:
                number_of_internal_nodes += 1
                number_of_edges += 2
                internal_node_id = "internalNode" + str(number_of_internal_nodes) + "X"
                edge1 = Edge("edge" + str(number_of_edges-1) + "X")
                edge1.set_node_ids([internal_node_id, hit.group(1)])
                edge1.set_dmv(float(hit.group(2)))
                edge1.set_length(float(hit.group(3)))
                edge2 = Edge("edge" + str(number_of_edges) + "X")
                edge2.set_node_ids([internal_node_id, hit.group(4)])
                edge2.set_dmv(float(hit.group(5)))
                edge2.set_length(float(hit.group(6)))
                self.edges.append(edge1)
                self.edges.append(edge2)
                working_newick_string = working_newick_string.replace(hit.group(0), internal_node_id)

        # Make sure the remaining string includes a single node and use this node id to determine root edges.
        pattern_rooted = re.compile("^internalNode\d+X$")
        hit_rooted = pattern_rooted.search(working_newick_string)
        if hit_rooted == None:
            print('ERROR: The newick tree string could not be parsed!')
            print('The remaining unparsed part of the newick string: \"{}\"'.format(working_newick_string))
            sys.exit(1)
        else:
            root_node_id = hit_rooted.group(0)
            for edge in self.get_edges():
                if edge.get_node_ids()[0] == root_node_id:
                    edge.set_is_root_edge(True)
                else:
                    edge.set_is_root_edge(False)

    def get_edges(self):
        return self.edges

    def get_number_of_edges(self):
        return len(self.edges)

    def get_number_of_extant_edges(self):
        number_of_extant_edges = 0
        for edge in self.edges:
            if edge.get_is_extant():
                number_of_extant_edges += 1
        return number_of_extant_edges

    def set_extant_progeny_ids(self):
        for edge in self.get_edges():
            if edge.get_is_extant():
                species_id = edge.get_node_ids()[1]
                this_edge = edge
                species_id_added_to_root_edge = False
                while species_id_added_to_root_edge == False:
                    this_edge.add_extant_progeny_id(species_id)
                    if this_edge.get_is_root_edge():
                        species_id_added_to_root_edge = True
                    else:
                        for other_edge in self.get_edges():
                            if other_edge.get_node_ids()[1] == this_edge.get_node_ids()[0]:
                                parent_edge = other_edge
                                break
                    this_edge = parent_edge

    def set_times(self):
        # Get the durations between root and extant edges.
        # These should be approximately similar, if not produce a warning.
        total_edge_lengths = []
        for edge in self.get_edges():
            if edge.get_is_extant():
                total_edge_length = edge.get_length()
                if edge.get_is_root_edge():
                    total_edge_lengths.append(total_edge_length)
                else:
                    root_edge_found = False
                    this_edge = edge
                    while root_edge_found == False:
                        for other_edge in self.get_edges():
                            if other_edge.get_node_ids()[1] == this_edge.get_node_ids()[0]:
                                parent_edge = other_edge
                                break
                        if parent_edge == None:
                            print('ERROR: The parent edge for edge {} could not be found'.format(this_edge.get_id()))
                            sys.exit(1)
                        total_edge_length += parent_edge.get_length()
                        if parent_edge.get_is_root_edge():
                            root_edge_found = True
                            total_edge_lengths.append(total_edge_length)
                        else:
                            this_edge = parent_edge
        max_total_edge_length = max(total_edge_lengths)
        if max_total_edge_length - min(total_edge_lengths) > 0.1:
            print('WARNING: The tree appears not to be ultrametric. Some terminal branches will be extended so that they all end at the same time.')
            print('')
        # Extend terminal edges if necessary.
        for edge in self.get_edges():
            if edge.get_is_extant():
                total_edge_length = edge.get_length()
                if edge.get_is_root_edge():
                    edge.set_length(edge.get_length() + max_total_edge_length - total_edge_length)
                else:
                    root_edge_found = False
                    this_edge = edge
                    while root_edge_found == False:
                        for other_edge in self.get_edges():
                            if other_edge.get_node_ids()[1] == this_edge.get_node_ids()[0]:
                                parent_edge = other_edge
                                break
                        if parent_edge == None:
                            print('ERROR: The parent edge for edge {} could not be found'.format(this_edge.get_id()))
                            sys.exit(1)
                        total_edge_length += parent_edge.get_length()
                        if parent_edge.get_is_root_edge():
                            root_edge_found = True
                            edge.set_length(round(edge.get_length() + max_total_edge_length - total_edge_length,8))
                        else:
                            this_edge = parent_edge
        # First specify the edges for which the parents still need to be identified.
        for edge in self.get_edges():
            if edge.get_is_root_edge():
                edge.set_parent_needs_times(False)
            else:
                edge.set_parent_needs_times(True)
        # Set the times of all edges.
        for edge in self.get_edges():
            if edge.get_is_extant() == True:
                edge.set_termination(0.0)
                edge.set_origin(edge.get_length())
                this_edge = edge
                while this_edge.get_parent_needs_times():
                    for other_edge in self.get_edges():
                        if other_edge.get_node_ids()[1] == this_edge.get_node_ids()[0]:
                            parent_edge = other_edge
                            break
                    if parent_edge == None:
                        print('ERROR: The parent edge for edge {} could not be found'.format(this_edge.get_id()))
                        sys.exit(1)
                    parent_edge.set_termination(this_edge.get_origin())
                    parent_edge.set_origin(this_edge.get_origin() + parent_edge.get_length())
                    this_edge.set_parent_needs_times = False
                    this_edge = parent_edge

    def get_origin(self):
        for edge in self.get_edges():
            if edge.get_is_root_edge():
                return edge.get_origin()

    def info(self):
        info_string = ''
        info_string += 'Tree'.ljust(20)
        info_string += '\n'
        info_string += 'Number of edges:'.ljust(20)
        info_string += str(self.get_number_of_edges())
        info_string += '\n'
        return info_string


# The Edge class.
class Edge(object):

    def __init__(self, id):
        self.id = id
        self.node_ids = []
        self.extant_progeny_ids = []
        self.origin = None
        self.termination = None
        self.parent_needs_times = True

    def get_id(self):
        return self.id

    def set_node_ids(self, node_ids):
        self.node_ids = node_ids

    def get_node_ids(self):
        return self.node_ids

    def get_is_extant(self):
        if self.node_ids[1][0:12] == 'internalNode':
            return False
        else:
            return True

    def set_length(self, length):
        self.length = length

    def get_length(self):
        return self.length

    def set_is_root_edge(self, is_root_edge):
        self.is_root_edge = is_root_edge

    def get_is_root_edge(self):
        return self.is_root_edge

    def add_extant_progeny_id(self, extant_progeny_id):
        self.extant_progeny_ids.append(extant_progeny_id)

    def get_extant_progeny_ids(self):
        return self.extant_progeny_ids

    def set_termination(self, termination):
        self.termination = termination

    def get_termination(self):
        return self.termination

    def set_origin(self, origin):
        self.origin = origin

    def get_origin(self):
        return self.origin

    def set_parent_needs_times(self, parent_needs_times):
        self.parent_needs_times = parent_needs_times

    def get_parent_needs_times(self):
        return self.parent_needs_times

    def set_dmv(self, dmv):
        self.dmv = dmv

    def get_pop_size(self, generation_time):
        return self.dmv * (1000000/float(generation_time))

    def info(self):
        info_string = ''
        info_string += 'Edge id:'.ljust(28)
        info_string += self.id
        info_string += '\n'
        info_string += 'Edge node 1 id:'.ljust(28)
        info_string += self.node_ids[0]
        info_string += '\n'
        info_string += 'Edge node 2 id:'.ljust(28)
        info_string += self.node_ids[1]
        info_string += '\n'
        info_string += 'Edge length:'.ljust(28)
        info_string += str(self.length)
        info_string += '\n'
        if self.dmv != None:
            info_string += 'Edge dmv:'.ljust(28)
            info_string += str(self.dmv)
            info_string += '\n'            
        info_string += 'Edge origin:'.ljust(28)
        info_string += str(self.origin)
        info_string += '\n'
        info_string += 'Edge termination:'.ljust(28)
        info_string += str(self.termination)
        info_string += '\n'
        info_string += 'Edge is extant:'.ljust(28)
        info_string += str(self.get_is_extant())
        info_string += '\n'
        info_string += 'Edge is root edge:'.ljust(28)
        info_string += str(self.is_root_edge)
        info_string += '\n'
        info_string += 'Edge extant progeny ids:'.ljust(28)
        for extant_progeny_id in self.extant_progeny_ids:
            info_string += '{}, '.format(extant_progeny_id)
        info_string = info_string[:-2]
        info_string += '\n'
        # info_string += 'Edge parent needs times:'.ljust(20)
        # info_string += str(self.get_parent_needs_times())
        # info_string += '\n'
        return info_string


def get_generations_per_branch_length_unit(
        branch_length_units=None,
        generation_time=None):
    """
    Method to calculate the number of generations per branch length
    unit, given the branch length unit and a generation time.
    """
    if branch_length_units == "gen":
        generations_per_branch_length_unit = 1
    elif branch_length_units == "myr":
        generations_per_branch_length_unit = 10**6/generation_time
    else:
        generations_per_branch_length_unit = 1/generation_time
    return generations_per_branch_length_unit


def parse_species_tree(
        species_tree=None,
        branch_length_units="gen",
        Ne=None,
        generation_time=None,
        migration_matrix=None,
        geneFlowPeriod=None):
    """
    Method to parse species trees in Newick
    (https://en.wikipedia.org/wiki/Newick_format) format.

    Trees are assumed to be rooted and ultrametric and branch lengths
    must be included and correspond to time, either in units of millions
    of years ("myr"), years ("yr"), or generations ("gen"; default).
    Leafs must be named. An example for an accepted tree string in
    Newick format is:
    (((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)
    The tree string can end with a semi-colon, but this is not required.

    - An estimate of the effective population size Ne should be
        specified.
    - If and only if the branch lengths are not in units of
        generations, the generation time should be specified.
    """

    # Make sure a species tree is specified.
    if type(species_tree) is not str:
        raise ValueError("A species tree must be specified.")

    # Make sure that branch length units are either "myr", "yr", or "gen".
    allowed_branch_lenth_units = ["myr", "yr", "gen"]
    if branch_length_units not in allowed_branch_lenth_units:
        err = 'The specified units for branch lengths ('
        err += '"{}") are not accepted. '.format(branch_length_units)
        err += 'Accepted units are "myr" (millions of years), "yr" (years), '
        err += 'and "gen" (generations).'
        raise ValueError(err)

    # Make sure that the population size is either None or positive.
    if Ne is not None:
        try:
            Ne = float(Ne)
        except ValueError:
            raise ValueError("Population size Ne must be numeric.")
        if Ne <= 0:
            raise ValueError("Population size Ne must be > 0.")

    # Make sure that the generation time is either None or positive.
    if generation_time is not None:
        try:
            generation_time = float(generation_time)
        except ValueError:
            raise ValueError("Generation time must be numeric.")
        if generation_time <= 0:
            raise ValueError("Generation time must be > 0.")

    # Make sure that the generation time is specified if and only if
    # branch lengths are not in units of generations.
    if branch_length_units == "gen":
        if generation_time is not None:
            err = 'With branch lengths in units of generations ("gen"), '
            err += 'a generation time should not be specified additionally.'
            raise ValueError(err)
    else:
        if generation_time is None:
            err = 'With branch lengths in units of '
            err += '"{}", a generation time must be '.format(branch_length_units)
            err += 'specified additionally.'
            raise ValueError(err)

    # Make sure that a population size is specified.
    if Ne is None:
        raise ValueError("Ne should be specified.")

    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        branch_length_units, generation_time
        )
    # Get the gene flow period limit parameter in the time units of the tree
    geneFlowPeriodInTreeUnits = geneFlowPeriod / generations_per_branch_length_unit

    # Read the input file.
    species_tree_lines = species_tree.splitlines(False)
    if len(species_tree_lines) > 1:
        raise ValueError("The species tree has multiple lines.")
    tree_patterns = re.search('\\(.+\\)', species_tree_lines[0])
    if tree_patterns is None:
        raise ValueError("No species tree string found.")
    species_tree_string = tree_patterns.group(0)

    # Parse the newick tree string.
    root = newick.loads(species_tree_string)[0]

    # Set node depths (distances from root).
    max_depth = 0
    for node in root.walk():
        depth = 0
        moving_node = node
        while moving_node.ancestor:
            depth += moving_node.length
            moving_node = moving_node.ancestor
        node.depth = depth
        if depth > max_depth:
            max_depth = depth

    # Set node heights (distances from present).
    for node in root.walk():
        node.height = max_depth - node.depth

    # Get a list of species IDs.
    species_ids = sorted(root.get_leaf_names())
    
    # Create a new matrix which will hold the migration rates for time 0.0 - to be supplied as a param to msprime.simulate
    timeZero_migration_matrix = [[0.0 for col in range(len(species_ids))] for row in range(len(species_ids))]

    # Determine at which time which populations should merge.
    sources = []
    destinations = []
    divergence_times = []
    allDescendantsOfDestinations = []
    allDescendantsOfSources = []
    for node in root.walk():
        if node.is_leaf is False:
            name_indices = []
            leaf_names_of_nodeDescendants = []
            for descendants in node.descendants:
                leaf_names = (descendants.get_leaf_names())
                name_indices.append(species_ids.index(sorted(leaf_names)[0]))
                leaf_names_of_nodeDescendants.append(sorted(leaf_names))
            
            s = sorted(zip(name_indices, leaf_names_of_nodeDescendants))
            name_indices, leaf_names_of_nodeDescendants = map(list, zip(*s))
            new_destination = name_indices[0]
            name_indices.remove(name_indices[0])
            newDescendantsOfDestination = leaf_names_of_nodeDescendants[0]
            leaf_names_of_nodeDescendants.remove(leaf_names_of_nodeDescendants[0])
            for x in range(len(name_indices)):
                sources.append(name_indices[x])
                allDescendantsOfSources.append(leaf_names_of_nodeDescendants[x])
                destinations.append(new_destination)
                allDescendantsOfDestinations.append(newDescendantsOfDestination)
                divergence_times.append(node.height)

    # Sort the lists source_sets, destinations, and divergence_times
    # according to divergence_time.
    s = sorted(zip(divergence_times, sources, destinations, allDescendantsOfSources, allDescendantsOfDestinations))
    divergence_times, sources, destinations, allDescendantsOfSources, allDescendantsOfDestinations = map(list, zip(*s))
    
   
   #DEBUG STUFF:
#    print("Batmin: " + str(species_ids.index("Batmin")))
#    print("Batleo: " + str(species_ids.index("Batleo")))
#    print("Batgra: " + str(species_ids.index("Batgra")))
#    print("Batfer: " + str(species_ids.index("Batfer")))
#    print("Bathor: " + str(species_ids.index("Bathor")))
#    print("Batfas: " + str(species_ids.index("Batfas")))
#    print("Batvit: " + str(species_ids.index("Batvit")))
    
    
    for x in range(len(divergence_times)):
        # Start of gene flow ('Start' in the sense of looking backwards in coalescent terms)
        geneFlowStart =  divergence_times[x] - geneFlowPeriodInTreeUnits
        
        if geneFlowStart < 0.0:
            geneFlowStart = 0.0
        else:
            descendantsToRemove = [] # These should not have any gene-flow, because they do not exist within the geneFlowPeriod for this node
            for s in allDescendantsOfSources[x]:
                speciesID = species_ids.index(s)
                for y in range(len(divergence_times)):
                    if speciesID == sources[y]:
                        if divergence_times[y] <= geneFlowStart:
                            descendantsToRemove.append(s)
            
            for r in descendantsToRemove:
                allDescendantsOfSources[x].remove(r)
            
            descendantsToRemove = []
            
            for s in allDescendantsOfDestinations[x]:
                speciesID = species_ids.index(s)
                for y in range(len(divergence_times)):
                    if speciesID == sources[y]:
                        if divergence_times[y] <= geneFlowStart:
                            descendantsToRemove.append(s)
                            
            for r in descendantsToRemove:
                allDescendantsOfDestinations[x].remove(r)
    

    # Define the species/population tree for msprime.
    population_configurations = []
    for _ in range(len(root.get_leaves())):
        population_configurations.append(
            msprime.PopulationConfiguration(
                initial_size=Ne))
    demographic_events = []
    demographic_event_times = []
    for x in range(len(divergence_times)):
        demographic_events.append(
            msprime.MassMigration(
                time=divergence_times[x]*generations_per_branch_length_unit,
                source=sources[x],
                destination=destinations[x]))
        demographic_event_times.append(divergence_times[x])
        
        
        # Start of gene flow ('Start' in the sense of looking backwards in coalescent terms)
        geneFlowStart =  divergence_times[x] - geneFlowPeriodInTreeUnits
    #DEBUG STUFF:
#        print("geneFlowPeriodInTreeUnits: ")
#        print(geneFlowPeriodInTreeUnits)
#        print("divergence_times[x]: ")
#        print(divergence_times[x])
#        print("geneFlowStart: ")
#        print(geneFlowStart)
        
        if geneFlowStart < 0.0:
            geneFlowStart = 0.0
            for d in allDescendantsOfDestinations[x]:
                for s in allDescendantsOfSources[x]:
                    timeZero_migration_matrix[species_ids.index(s)][species_ids.index(d)] = migration_matrix[species_ids.index(s)][species_ids.index(d)]
                    timeZero_migration_matrix[species_ids.index(d)][species_ids.index(s)] = migration_matrix[species_ids.index(d)][species_ids.index(s)]
            
 #DEBUG STUFF:
#        print("geneFlowStart: ")
#        print(geneFlowStart)
#
#        print("")
        else:
        # Now loop over all descendants of this node which should have gene-flow among them and set that
            for d in allDescendantsOfDestinations[x]:
                for s in allDescendantsOfSources[x]:
        #DEBUG STUFF:
    #               if divergence_times[x] > 1.5:
    #                    print(d)
    #                    print(s)
    #                    print(species_ids[8])
    #                    print(species_ids[67])
    #                    print(species_ids[68])
    #                    sys.exit(1)
                    if migration_matrix[species_ids.index(s)][species_ids.index(d)] > 0.0:
                        demographic_events.append(
                        msprime.MigrationRateChange(
                            time=geneFlowStart*generations_per_branch_length_unit,
                            rate=migration_matrix[species_ids.index(s)][species_ids.index(d)],
                            matrix_index=(species_ids.index(s), species_ids.index(d))))
                        demographic_event_times.append(geneFlowStart)
                        
                    if migration_matrix[species_ids.index(d)][species_ids.index(s)] > 0.0:
                        demographic_events.append(
                        msprime.MigrationRateChange(
                            time=geneFlowStart*generations_per_branch_length_unit,
                            rate=migration_matrix[species_ids.index(d)][species_ids.index(s)],
                            matrix_index=(species_ids.index(d), species_ids.index(s))))
                            
                        demographic_event_times.append(geneFlowStart)
                    
                    # Of course better if the rates in internal branches were averages, not just randomply selected 'surviving' ones
        
        for y in range(len(root.get_leaves())):
            if y != sources[x]:
                demographic_events.append(
                    msprime.MigrationRateChange(
                        time=divergence_times[x]*generations_per_branch_length_unit,
                        rate=0,
                        matrix_index=(sources[x], y)))
                demographic_event_times.append(divergence_times[x])
                
                demographic_events.append(
                    msprime.MigrationRateChange(
                        time=divergence_times[x]*generations_per_branch_length_unit,
                        rate=0,
                        matrix_index=(y, sources[x])))
                demographic_event_times.append(divergence_times[x])
                
    
    # Sort the list demographic_events according to demographic_event_times (the "zip" approach does't work on the type in demographic_events)
    if len(demographic_event_times) == len(demographic_events):
        desiredOrder = numpy.argsort(demographic_event_times)
        sorted_demographic_events = []
        for x in range(len(desiredOrder)):
            sorted_demographic_events.append(demographic_events[desiredOrder[x]])
    else:
        print("Something went horribly wrong: len(demographic_event_times) != len(demographic_events)")
        sys.exit(1)

    # Return a tuple of population_configurations, demographic_events, and the initial migration_matrix at time 0.0
    return population_configurations, sorted_demographic_events, timeZero_migration_matrix

# Import libraries.
import msprime
import sys
import newick
import random
from random import randint
import datetime
import re
import os
import numpy
import argparse

# Get the command line arguments/parameters
argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="Build the msprime model and run simulations given a tree and a migration matrix based on f4-ratio statistics produced by Dsuite. ", add_help=True)
argparser.add_argument("tree", type = str, help="Path to a .newick tree file")
argparser.add_argument("f4ratios", type=str, help="Path to file containing f4-ratio statitics matrix derived from Dsuite results")
argparser.add_argument("out", type=str, help="Prefix to the output file which will contain the simulation results")
argparser.add_argument("-g", "--gen_time", type=float, help="Average generation time in years.", default=3.0)
argparser.add_argument("-u", "--mut_rate", type=float, help="Mutation rate per bp per generation",default=3.5e-9)
argparser.add_argument("-r", "--rho", type=float, help="Recombination rate per bp per generation (in the future could use a recombination map); set to 0 to simulate a non-recombining region", default=2.2e-8)
argparser.add_argument("-N", "--Ne", type=int, help="Effective population size (for now this is the same throughout the tree.", default=20000)
argparser.add_argument("-l","--chr_length", type=int, help="Length of the 'chromosome' to simulate", default=100)
argparser.add_argument("-p", "--gene_flow_period", help="How long gene-flow persist after species/populations split (in generations).", type=int, default=1000000)
argparser.add_argument("-s", "--scaling_factor", help="A scaling factor between the f4-ratios and the msprime migration units", type=float, default=1.0e-5)
argparser.add_argument("-n", "--num_indiv", help="The number of individuals to simulate from each population/species", type=int, default=1)
argparser.add_argument("--ploidy", help="The ploidy of simulated individuals", type=int, default=2)
argparser.add_argument("-d", "--debugger", help="Save the msprime.DemographyDebugger output to FILENAME", type=str, default=None,metavar="FILENAME")
argparser.add_argument("--species_list", type=str, help="Path to a file with a list of species to subsample")
args = argparser.parse_args()

# Read the species tree string.
with open(args.tree, "r") as f:
    species_tree_string = f.read()

# Parse the species tree.
root = newick.loads(species_tree_string)[0]
species_ids = sorted(root.get_leaf_names())

# Read the species tree string.
sub_sample_species_list = []
with open(args.species_list, "r") as f:
    for line in f.readlines():
        sub_sample_species_list.append(line.strip())

# Check if the specified migration rate string is a number or a file name.
print('Preparing the migration matrix...', end='', flush=True)
if os.path.isfile(args.f4ratios):
    migrations = []
    # Read the migration rate matrix.
    with open(args.f4ratios) as f:
        dsuite_file_lines = f.readlines()
        dsuite_file_species_ids = dsuite_file_lines[0].split()
        migration_matrix = []
        for species1 in species_ids:
            row = []
            if species1 in dsuite_file_species_ids:
                dsuite_species1_index = dsuite_file_species_ids.index(species1)
            else:
                dsuite_species1_index = None
          
            for species2 in species_ids:
                if species2 in dsuite_file_species_ids:
                    dsuite_species2_index = dsuite_file_species_ids.index(species2)
                else:
                    dsuite_species2_index = None
   
                if species1 == species2:
                    row.append(0)
                else:
                    if dsuite_species1_index is None or dsuite_species2_index is None:
                        row.append(0)
                    else:
                        row_list = dsuite_file_lines[dsuite_species1_index+1].split()
                        migration = float(row_list[dsuite_species2_index+1]) * args.scaling_factor
                        row.append(migration)
                        migrations.append(migration)

            migration_matrix.append(row)
    print(" done. Mean migration rate is " + str(sum(migrations)/len(migrations)) + ".")

    
#print(migration_matrix)
    
# Parse the species tree with msprime and generate population configurations and demographic evens.
parsed_tuple = parse_species_tree(
    species_tree=species_tree_string,
    branch_length_units="myr",
    Ne=args.Ne,
    generation_time=args.gen_time,
    migration_matrix=migration_matrix,
    geneFlowPeriod=args.gene_flow_period
    )
population_configurations = parsed_tuple[0]
demographic_events = parsed_tuple[1]
timeZero_migration_matrix = parsed_tuple[2]
for n in population_configurations:
    n.sample_size = args.ploidy * args.num_indiv

if not args.debugger is None:
    dd = msprime.DemographyDebugger(
            population_configurations=parsed_tuple[0],
            demographic_events=parsed_tuple[1],
            migration_matrix=timeZero_migration_matrix)
    dd_outfile = open(args.debugger, 'w')
    dd.print_history(dd_outfile)
#sys.exit(0)

# Simulate tree sequences.
print('Simulating with msprime...', end='', flush=True)
ts = msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=timeZero_migration_matrix,
        demographic_events=demographic_events,
        mutation_rate=args.mut_rate,
        length=args.chr_length,
        recombination_rate=args.rho,
        random_seed=random.randint(1, 10000000)
        )
print(" done.")

# Subsample the tree sequence to include only those species listed in the species list.
sub_sample_indices = []
for sub_sample_species in sub_sample_species_list:
    if sub_sample_species in species_ids:
        species_ids_index=species_ids.index(sub_sample_species)-1
        for x in range(args.num_indiv):
            for y in range(args.ploidy):
                sub_sample_indices.append(args.ploidy*args.num_indiv*species_ids_index + args.ploidy*x + y + 1)
    else:
        print("ERROR: Species id " + sub_sample_species + " could not be found in the list of species ids!")
        sys.exit(1)
sub_ts = ts.simplify(sub_sample_indices)

# Run c-genie code.
tree_sequences = []
for tree in sub_ts.trees():
    this_tree_to_append = []
    this_tree_to_append.append(tree.interval[1]-tree.interval[0])
    this_tree_to_append.append(tree.newick())
    tree_sequences.append(this_tree_to_append)

# Analyse msprime results.
c_gene_sizes = []
tract_sizes = []
tree_strings = []
topology_strings = []
for length_and_tree in tree_sequences:
    length = int(length_and_tree[0])
    tree_string = length_and_tree[1].rstrip(';')
    topology_string = re.sub(r':\d+\.\d+','', tree_string)
    c_gene_sizes.append(length)
    tree_strings.append(tree_string)
    topology_strings.append(topology_string)

# Compare each tree with the one before to identify topology breakpoints.
previous_topology_string = topology_strings[0]
previous_topology_length = c_gene_sizes[0]
num_topology_changes = 0
for x in range(1, len(topology_strings)):
    # Only parse trees if the topology strings differ.
    if topology_strings[x] != previous_topology_string:
        # Make sure the two trees are really different.
        tree1 = Tree(tree_strings[x-1])
        tree1.parse_newick_string()
        tree1.set_extant_progeny_ids()
        tree2 = Tree(tree_strings[x])
        tree2.parse_newick_string()
        tree2.set_extant_progeny_ids()
        clades1 = []
        for edge in tree1.get_edges():
            extant_progeny_ids = edge.get_extant_progeny_ids()
            extant_progeny_ids.sort()
            if len(extant_progeny_ids) > 1:
                clades1.append(extant_progeny_ids)
        clades1.sort()
        clades2 = []
        for edge in tree2.get_edges():
            extant_progeny_ids = edge.get_extant_progeny_ids()
            extant_progeny_ids.sort()
            if len(extant_progeny_ids) > 1:
                clades2.append(extant_progeny_ids)
        clades2.sort()
        if clades1 != clades2:
            num_topology_changes = num_topology_changes + 1
            tract_sizes.append(previous_topology_length)
            previous_topology_length = c_gene_sizes[x]
        else:
            previous_topology_length += c_gene_sizes[x]
    else:
        previous_topology_length += c_gene_sizes[x]
    previous_topology_string = topology_strings[x]

# If there have not been any topology changes, set the tract lengths to equal the lengths of the simulated chromosomes.
if num_topology_changes == 0:
    tract_sizes = [args.chr_length]

# Remove the last c-gene size as it is truncated not by recombination, but by the chromosome end.
c_gene_sizes = c_gene_sizes[:-1]

# Write output files.
c_gene_size_string = ""
for c_gene_size in c_gene_sizes:
    c_gene_size_string += str(c_gene_size) + "\n"
with open(args.out + "_c-genes.txt", "w") as c_genes_file:
    c_genes_file.write(c_gene_size_string)
tract_size_string = ""
for tract_size in tract_sizes:
    tract_size_string += str(tract_size) + "\n"
with open (args.out + "_tracts.txt", "w") as tracts_file:
    tracts_file.write(tract_size_string)

# Screen output.
print("Wrote file " + args.out + "_c-genes.txt.")
print("Wrote file " + args.out + "_tracts.txt.")
