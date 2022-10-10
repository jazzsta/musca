from openalea.mtg import mtg, io
from past.utils import old_div
import numpy as np
from openalea.mtg import algo
import time

g = io.read_mtg_file("Med_cul10_1996_testcorrect.mtg")


def compute_path(g,i,j, with_gca=True):
    """
    Compute topological path between two vertices

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'i', 'j' - vertices between which computing the path

    :Returns:
        topological path between verteces 'i' and 'j', excluding the two input vertices.

    :Remark:
        - When a branching point, in common between the ancestors of the vertices, is encountered, it is removed from the path.

    """


    list1 = g.Ancestors(i)
    list2 = g.Ancestors(j)

    _gca_id = algo.lowestCommonAncestor(g, [i,j])
    position1 = list1.index(_gca_id)
    position2 = list2.index(_gca_id)

    path = list1[:position1]+[_gca_id]+list(reversed(list2[:position2]))
    # uniqueAncestors = list(set(g.Ancestors(i)) ^ set(g.Ancestors(j))) #Ancestors not in common between i and j

    # if i in uniqueAncestors:
    #     uniqueAncestors.remove(i)
    # if j in uniqueAncestors:
    #     uniqueAncestors.remove(j)
    # path = uniqueAncestors
    return(path)

def bary_base(g, selected_scale):
    """
    Compute barycentres and basis of vertices of mtg at selected scale.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'selected_scale' - scale at which calculating baricentres and basis

    :Remark:
        - it calls the function 'barycentre' and 'base'
    """

    if selected_scale < g.max_scale():
        bary_base(g, selected_scale+1)

    if selected_scale == 1:
        pass
    else:
        for i in g.vertices(scale=selected_scale):
            barycentre(g, i)
            base(g, i)
    print("barycentres and basis computed in scale "+str(selected_scale))

def vertices_semilength(g,selected_scale):
    """
    Compute euclidean distance between base and baricentre of each vertex at the selected scale.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'selected_scale' - scale at which computing the semilength of vertices
    :Returns:
        g with property 'semilength': euclidean distance between barycentre and base of each vertex.

    :Remark:
        - It calls the function: geom_distance
    """
    print("Compute vertices_semilength at scale ", selected_scale)
    for i in g.vertices(scale=selected_scale):
        node_i = g.node(i)
        #print "calcola semilength per vertice",i
        node_i.semilength = geom_distance(g,i,i,option="base_bary")
    print("vertices_semilength computed ")

def barycentre(g, i):
    """
    Compute barycenter of a vertex.
    This is given by:
    (i) for a vertex at the finest scale: average of the spatial coordinates of a vertex and its parent;
    (ii) for vertex at coarse scale: the average of the spatial coordinates of the barycenters of its components, weighted by their lengths.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'i' - vertex for which calculating baricentre

    :Returns:
        stores barycentre for the vertex 'i' in: barycentre_XX, barycentre_YY, barycentre_ZZ.

    """
    # compute baricentre of the vertex, depending if it is at max_scale or a coarse scale
    if g.scale(i)==g.max_scale():
        # for finest scale
        #print "fine scale"
        if g.node(g.parent(i)):
            # for vertices other than the base
            g.node(i).barycentre_XX = old_div((g.node(i).XX + g.node(g.parent(i)).XX), 2)
            g.node(i).barycentre_YY = old_div((g.node(i).YY + g.node(g.parent(i)).YY), 2)
            g.node(i).barycentre_ZZ = old_div((g.node(i).ZZ + g.node(g.parent(i)).ZZ), 2)
            #print "internal vertex"
        else:
            # for the base: find the basal coordinates in the complex
            g.node(i).barycentre_XX = old_div((g.node(i).XX + g.node(g.complex(i)).XX), 2)
            g.node(i).barycentre_YY = old_div((g.node(i).YY + g.node(g.complex(i)).YY), 2)
            g.node(i).barycentre_ZZ = old_div((g.node(i).ZZ + g.node(g.complex(i)).ZZ), 2)
            #print "extremity"

    else:
        # for a corse scale
        X,Y,Z,tot_length = 0,0,0,0
        for k in g.components(i):
            node = g.node(k)
            # X += node.XX * node.length
            # Y += node.YY * node.length
            # Z += node.ZZ * node.length
            X += node.barycentre_XX * node.length
            Y += node.barycentre_YY * node.length
            Z += node.barycentre_ZZ * node.length
            tot_length += node.length
        g.node(i).barycentre_XX = old_div(X,tot_length)
        g.node(i).barycentre_YY = old_div(Y,tot_length)
        g.node(i).barycentre_ZZ = old_div(Z,tot_length)
        #print "coarse scale"

def base(g, i):
    """
    Compute base of a vertex.
    This is given by:
    (i) for a vertex at the finest scale:
            for the basal vertex: spatial coordinates (XX,YY,ZZ) stored in the complex (for the basal vertex);
            for other vertices: the spatial coordinates of the parent;
    (ii) for a vertex at a coarse scale:
            for the basal vertex: the spatial coordinates (XX,YY,ZZ) stored in the same vertex;
            for other vertices: the spatial coordinates of the predecessor of the first vertex of its components.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'i' - vertex for which defining the base

    :Returns:
        stores base for the vertex 'i' in: base_XX, base_YY, base_ZZ.

    """
    # compute basal coordinates of the vertex, depending if it is at max_scale or a coarse scale
    if g.scale(i)==g.max_scale():
        # for finest scale
        #print "vertex at fine scale"
        if not(g.node(g.parent(i))):
            # for the base: you will find the basal coordinates in the complex
            g.node(i).base_XX = g.node(g.complex(i)).XX #(g.node(i).XX + g.node(g.parent(i)).XX) / 2
            g.node(i).base_YY = g.node(g.complex(i)).YY
            g.node(i).base_ZZ = g.node(g.complex(i)).ZZ
            #print "base"
        else:
            # for vertices other than the base, basal coordinates are the (apical) coordinates of the parent
            g.node(i).base_XX = g.node(g.parent(i)).XX
            g.node(i).base_YY = g.node(g.parent(i)).YY
            g.node(i).base_ZZ = g.node(g.parent(i)).ZZ
            #print "not the base"
    else:
        # for a corse scale
        #print "vertex at coarse scale"
        if not(g.node(g.parent(i))):
            # for the base: find the basal coordinates in the (XX,YY,ZZ) of the same vertex
            g.node(i).base_XX = g.node(i).XX
            g.node(i).base_YY = g.node(i).YY
            g.node(i).base_ZZ = g.node(i).ZZ
            #print "base"
        else:
            # for other non base vertices:
            # find the basal element of this coarse scale vertex: it is the vertex whose parent is not a component of the coarse scale vertex
            # and take coordinates of its parent
            components_i = g.components(i)
            for j in components_i:
                if g.parent(j) not in components_i:
                    g.node(i).base_XX = g.node(g.parent(j)).XX
                    g.node(i).base_YY = g.node(g.parent(j)).YY
                    g.node(i).base_ZZ = g.node(g.parent(j)).ZZ
                    #print "not base"

def euclidean_distance(Xs,Ys,Zs):
    """
    Compute euclidean distance between couples of spatial coordinates by means of the pitagora theorem applied to the X, Y, Z coordinates.

    :Parameters:
        - 'Xs, Ys, Zs' - couples of x,y,z spatial coordinates

    :Returns:
        euclidean distance between the points defined by the X,Y,Z couples of spatial coordinates
    """
    import math
    #print Xs, Ys, Zs
    euclidean_distance = math.sqrt((Xs[0] - Xs[1])**2 + (Ys[0] - Ys[1])**2 + (Zs[0] - Zs[1])**2)
    return(euclidean_distance)

def geom_distance(g,i,j,option="barycentres"):
    """
    Compute euclidean distance between bases and/or baricentre of one/two vertices based on their spatial coordinates.
    This is given by the straight line that connects them,
    and is calculated by means of the pitagora theorem applied to the X, Y, Z coordinates.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'i', 'j' - vertices between which computing the distance
        - 'option' - 'barycentres'(default) the distance is calculated among barycentres
                     'basis',  the distance is calculated among basis
                     'base_bary',  the distance is calculated among barycentre of i and base of j.

    :Returns:
        - 'distance'- euclidean distance between barycentres/basis/base and barycentre of 'i' and 'j'.

    :Remark:
        - It calls the function: euclidean_distance

    """
    node_i = g.node(i)
    node_j = g.node(j)
    #print ("i",i,"; j",j)
    if option == "barycentres":
        Xs = node_i.barycentre_XX, node_j.barycentre_XX
        Ys = node_i.barycentre_YY, node_j.barycentre_YY
        Zs = node_i.barycentre_ZZ, node_j.barycentre_ZZ
    elif option == "basis":
        Xs = node_i.base_XX, node_j.base_XX
        Ys = node_i.base_YY, node_j.base_YY
        Zs = node_i.base_ZZ, node_j.base_ZZ
    elif option == "base_bary":
        Xs = node_i.barycentre_XX, node_j.base_XX
        Ys = node_i.barycentre_YY, node_j.base_YY
        Zs = node_i.barycentre_ZZ, node_j.base_ZZ
    distance = euclidean_distance(Xs,Ys,Zs)
    return (distance)

def gca(g,i,j, last_ancestors=False):
    """
    Finds Greatest Commen Ancestor (GCA) of i and j.
    Can compute the last ancestors of i and j before GCA.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'i', 'j' - vertices for which calculating GCA
        - last_ancestors (default = FALSE) - if TRUE: the last ancestors of i and j before GCA are also returned.

    :Returns:
        GCA of vertices 'i', 'j'.

    """
    # Greatest Common Ancestors
    ancestors_in_common = list(set(g.Ancestors(i)).intersection(g.Ancestors(j)))
    root_at_current_scale = g.roots( scale=g.scale(i) )[0]

    # find oldest ancestor between i and GCA
    # if one vertex is the base of the mtg, it is the CGA
    if i == root_at_current_scale:
        GCA = i
        last_ancestor_i = i
    else:
        found_i = False
        for k in g.Ancestors(i):
            if found_i == False:
                #print k,i,j
                if g.parent(k) in ancestors_in_common:
                    found_i = True
                    last_ancestor_i = k
                    # print last_ancestor_i
        GCA = g.parent(last_ancestor_i)

    if last_ancestors:
        if j == root_at_current_scale:
            last_ancestor_j = j
        else:
            found_j = False
            for k in g.Ancestors(j):
                if found_j == False:
                    if g.parent(k) in ancestors_in_common:
                        found_j = True
                        last_ancestor_j = k
                        # print last_ancestor_j
        GCA = GCA, last_ancestor_i, last_ancestor_j
    return(GCA)

def compute_distance(g,i,j):
    """
    Compute geometrical distance between the barycentres of two vertices.
    This is given by the sum of:
    - the euclidean distance between the barycentre and the base of each extremity (i,j),
    - the sum of the euclidean distance between couple of basis along the path that goes
        from each extremity to their Greater Common Ancestor (the latter excluded),
    - the euclidean distance between insertion points of the first ancestors of the two vertices before the GCA.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'i', 'j' - vertices between which computing the path

    :Returns:
        geometrical distance between two vertices

    :Remark:
        - It calls the functions: gca

    """
    node_i = g.node(i)
    node_j = g.node(j)
    if i==j:
        distance = 0
    else:
        # the distance is given by the sum of the semilength of the extremities
        distance = node_i.semilength + node_j.semilength

        # if the vertices are not neighbors add also...
        #print j,g.Sons(i),i,g.Sons(j)
        if ( j not in g.Sons(i) ) & ( i not in g.Sons(j) ):
            # retrieve the Greatest common ancestor of vertices i,j
            GCA = gca(g,i,j)

            # from vertex i to GCA, add the distance (from base to base) between each couple of vertices (except between GCA-1 and GCA)
            if i != GCA:
                # this cycle iteratively adds to distance between n-th and n-th-1 parent of i, getting progressively closer to the GCA
                while g.parent(i) != GCA:
                    distance += geom_distance(g,i,g.parent(i),option="basis")
                    i = g.parent(i)

            # as above but for j
            if j != GCA:
                while g.parent(j) != GCA:
                    distance += geom_distance(g,j,g.parent(j),option="basis")
                    j = g.parent(j)

            # in case both i and j are sons of GCA, need to compute distance between their basis (i = GCA-1, j=CGA-1, both inserted in GCA)
            #print i,j,GCA,g.parent(i),g.parent(j)
            if ( j in g.Sons(GCA) ) & ( i in g.Sons(GCA) ):
                # add the distance between basis of first sons of GCA, if
                distance += geom_distance(g,i,j,option="basis")
    return(distance)

def distances_baricentre(g, selected_scale):
    """
    Computes a matrix of distances between vertices at scale 'selected_scale'

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'selected_scale' (int) - topological scale of the plant components between which computing distances

    :Returns:
        - Dists - a matrix of distances

    :Remark:
        - it calls the function 'compute_geometrical_distance'

    """
    # Distances leaves-demands
    Dists = DistanceMatrix(g.vertices(scale=selected_scale))
    for Leaf in g.vertices(scale=selected_scale):
        for Demand in g.vertices(scale=selected_scale):
            # if Leaf != Demand:
            #     Dists[Leaf,Demand] = compute_geometrical_distance(g,Leaf,Demand)
            # else:
            #     Dists[Leaf,Demand] = 0.5 * compute_geometrical_distance(g,Leaf,Demand)
            Dists[Leaf,Demand] = compute_geometrical_distance(g,Leaf,Demand, selected_scale)                # substituted to be able to have distance from component to same component = 0
    return Dists

def compute_geometrical_distance(g,i,j, selected_scale):
    """
    Compute geometrical distance between two vertices by summing length of topological path components.
    This is given by the sum of the lengths of components in the path that connects them,
    plus half the length of each one of the extremities.
    In case a coarse scale is used, the length of extremities is the length of their main axes.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'i', 'j' - vertices between which computing the path

    :Returns:
        geometrical distance between 'i' and 'j', excluding the input vertices.

    :Remark:
        - It includes half of the length of each input vtx.
        - It calls length_axis function

    """
    path = compute_path(g,i,j)
    length = g.property('length')
    geometrical_distance = 0

    if i!=j:                                            #introduced to have distance from component to same component = 0
        if selected_scale == g.max_scale():
            Length_i = length[i]
            Length_j = length[j]
        else:
            Length_i = length_axis(g, i)
            Length_j = length_axis(g, j)

        geometrical_distance = old_div((Length_i + Length_j),2.)

    else:                                                 #introduced to have distance from component to same component = 0
        geometrical_distance = 0                          #introduced to have distance from component to same component = 0
    for k in path :
        #node = g.vertices(k)
        geometrical_distance += length[k]
    return(geometrical_distance)



def length_axis(g, vid):
    """
    Aggregate length of the axis of a coarse scale component.
    It sums the lengths of the elements at fine scale connected by "<" to the first element of the coarse scale component.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'vid' (int) - vertex at coarse scale for which the length of the main axis will be aggregated

    :Returns:
        The length of the main axis of component vid.

    """
    #print "Compute Length_axis", vid
    i = g.components(vid)[0] #selects the first metamer of the coarse scale
    A = True
    Length = 0
    while A:
        if g.Sons(i, EdgeType = "<"):
            Length  += g.node(i).length #print "Prima" #
            i = g.Sons(i)[0]
        else:
            A = False
    g.node(vid).length = Length
    #print Length
    return(Length)

def aggregate_axes_length(g, selected_scale):
    """
    Computes the length of coarse scale components of an mtg
    as the length of their axes

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'selected_scale' (int) - Scale to which aggregating the properties

    :Returns:
        An MTG with the property 'length' aggregated at the coarse scale

    :Remark:
    - it calls the function 'length_axis'

    """
    print("Aggregate_axes_length")
    for i in g.vertices(scale=selected_scale):
        length_axis(g, i)


class DistanceMatrix(object):
    def __init__(self, vtxlist):
        self.idmap = dict([(vid,i) for i,vid in enumerate(vtxlist)])
        size = len(vtxlist)
        self.values = np.zeros((size, size))

    def __getitem__(self, vids):
        v1, v2 = vids
        return self.values[self.idmap[v1]][self.idmap[v2]]

    def __setitem__(self, vids, value ):
        v1, v2 = vids
        self.values[self.idmap[v1]][self.idmap[v2]] = value

def distances(g, selected_scale, geom=False):
    """
    Computes a matrix of distances between vertices at scale 'selected_scale'

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'selected_scale' (int) - topological scale of the plant components between which computing distances

    :Returns:
        - Dists - a matrix of distances

    :Remark:
        - it calls the function 'compute_geometrical_distance'

    """
    # Distances leaves-demands
    print("Compute matrix of distances:")
    startime = time.time()
    Dists = DistanceMatrix(g.vertices(scale=selected_scale))
    nb_vertices = len(g.vertices(scale=selected_scale))
    nb_vertices_done = 0
    score_init = 0
    for Leaf in g.vertices(scale=selected_scale):
        for Demand in g.vertices(scale=selected_scale):
            if geom:
                Dists[Leaf,Demand] = compute_geometrical_distance(g,Leaf,Demand, selected_scale)                # substituted to be able to have distance from component to same component = 0
            else:
                Dists[Leaf,Demand] = compute_distance(g,Leaf,Demand)

        nb_vertices_done += 1
        score = sum((old_div((nb_vertices_done*100),nb_vertices)) >= np.array([ 10,  20,  30,  40, 50,  60, 70, 80, 90, 100]))
        if score > score_init:
            print("computed",(old_div((nb_vertices_done*100),nb_vertices)), "% of vertices (", nb_vertices_done,")")
            score_init = score
    print("Time taken:", time.time()-startime)
    return Dists

def C_allocation_SIMWAL2ses(g, h, selected_scale, Dists, group=False): #work in progress
    """
    Computes carbon allocation among plant components (from Leaf to Demand),
    carbon for growth and in excess at the 'selected_scale'.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'h' (real) - friction parameter
        - 'selected_scale' (int) - topological scale of the plant components between which carbon allocation is computed

    :Returns:
        an MTG with attributes
        'C_allocation' - carbon allocated to a component,
        'C_for_growth' - carbon available for growth in a component,
        'C_excess' - carbon allocated in excess in respect to the demand

    :Remark:
        - it calls the distance matrix 'Dists'

    """
#    beginning_time = time.time()
    petiole_length = 0   # this parameter is later added to the supply-sink distance in order account for petiole length
    supply = g.property('C_supply')
    demand = g.property('C_demand')

    SimW_Influx = dict.fromkeys(g.vertices(), 0) #Init      #Dictionaries to record fluxes (3 lines)
    SimW_Exflux = dict.fromkeys(g.vertices(), 0) #Init
    SimW_Excess = dict.fromkeys(g.vertices(), 0) #Init

    total_supply = sum(supply.values())                     #Compute total Demand & Supply
    total_demand = sum(demand.values())
    SimW_Tot_Influx = 0
    SimW_Tot_Excess = 0


    # Denumerators for SIMWAL: for each Leaf: SUM of all (Demand * distance)
    DEN = {vid:1 for vid in g.vertices()}

    nodes = set(g.vertices(scale=selected_scale))              # !! but the ID of new nodes is higher than the length of the vector of the nodes
    for Leaf in nodes: #supply:
        if True:#g.node(Leaf).leaf_area > 0:
            geo_dem_i_k = 0
            for Demand, dem_k in demand.items():   #Vertex i and its demand.
                if Demand in nodes:
###                    if (dem_k > 0):                      #for vertices with a demand ....all vertices have a demand
                    geo_distance_i_k = Dists[Leaf,Demand] + petiole_length    # with petiole_length > 0, there is a distance between a leaf and its own metamer
                    #print geo_distance_i_k
                    #geo_dem = dem_k * (geo_distance_i_k**(-h))               #substituted to have distance from component to same component = 0
                    geo_dem = dem_k * ( (old_div(1, (1 + geo_distance_i_k)))**h )    #introduced to have distance from component to same component = 0
                    geo_dem_i_k += geo_dem
                    # print "geo_dem", geo_dem
                    #if not geo_dem > -1 and geo_dem <-1 or geo_dem==0:
                    #print "leaf", Leaf, "demand", Demand,"geo_dem", geo_dem
#            print DEN, nodes
            DEN[Leaf] = geo_dem_i_k
            #print "Leaf", Leaf, "Demand", Demand, "Denumerator", DEN[Leaf]

    if (not group) or total_supply < total_demand: # Use SIMWAL if group ==TRUE or total demand
        print("C_allocation_SIMWAL2ses: computes C allocation, for growth and in excess at scale", selected_scale)

        for Demand, dem_j in demand.items():
            if Demand in nodes:
                #print "demand of", Demand, "is", dem_j
                SimW_Flux_in_j = 0
                if dem_j > 0.:               # DELETE: all vertices have a demand      # for vertices with a demand
#                    print "Demand", Demand, "receives from"
                    for Leaf in supply:
                        if Leaf in nodes:
                            if (supply[Leaf]>0):    # for vertices with a positive supply
                                geo_distance_i_j =  petiole_length
#                                if (Leaf!= Demand):
                                    #COMMON to SIMWAL and QUALITREE                        #CALCULATE Numerator of Flux from vid to vjd
                                geo_distance_i_j +=  Dists[Leaf,Demand]
                                #CALCULATE Flux from vid to vjd
                                #SimW_Flux_i_j = (dem_j * (geo_distance_i_j ** (-h)) * supply[Leaf]) / DEN[Leaf]    #substituted to have distance from component to same component = 0
                                SimW_Flux_i_j = old_div((dem_j * ( (old_div(1, (1 + geo_distance_i_j))) **h ) * supply[Leaf]), DEN[Leaf])    #introduced to have distance from component to same component = 0
                                SimW_Flux_in_j   += SimW_Flux_i_j         #CALCULATE Total Flux Into vjd,     which will be obtained by the end of the cycle
                                SimW_Exflux[Leaf] += SimW_Flux_i_j         #CALCULATE Total Efflux out of vid, which will be obtained by the end of the cycle
                                g.node(Leaf).C_leftover -= SimW_Flux_i_j   # used to verify that C leaving a source is <= C_supply of that source
                                #print DEN[Leaf]
###                                print "Leaf", Leaf, "this much", SimW_Flux_i_j, "with a distance of", Dists[Leaf][Demand], "geo_distance", geo_distance_i_j, "DEN",DEN[Leaf]
#                                print "dem", dem_j*10000, "from supply", supply[Leaf]*10000, "geo_distance_i_j", geo_distance_i_j, "obtains", SimW_Flux_i_j*10000
                g.node(Demand).C_allocation += SimW_Flux_in_j     # Total Amount allocated to a Vertex
                #print "allocated", g.node(Demand).C_allocation
                g.node(Demand).C_excess = max(0, SimW_Flux_in_j - dem_j) # carbon allocated in excess in respect to demand
                g.node(Demand).C_for_growth = g.node(Demand).C_allocation - g.node(Demand).C_excess
                #print "g.node(Demand).C_allocation",g.node(Demand).C_allocation, " - g.node(Demand).C_excess", g.node(Demand).C_excess

                SimW_Excess[Demand] = max(0, SimW_Flux_in_j - dem_j)
                SimW_Tot_Influx += SimW_Influx[Demand]
                SimW_Tot_Excess += SimW_Excess[Demand]

    else:                  # Use common Pool Assumption
        print("Common Pool Assumption")
        SimW_Tot_Excess = total_supply - total_demand
        print("excess", SimW_Tot_Excess)
        for Demand, dem_j in demand.items(): #vertex j and its demand
            if Demand in nodes:
                SimW_Influx[Demand]  = dem_j  #Total flow in j for SIMWAL
#    elapsed_time = time.time()- beginning_time
#    print elapsed_time
#    import winsound
    #    winsound.Beep(400,100)

def clean_properties(g, selected_scale):
    """
    Clean properties (C_allocation, C_leftover) at scale different from selected scale.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'selected_scale' (int) - topological scale of the plant components between which carbon allocation is computed
    """
    print("Set C-allocation and C-leftover for vertexes not in scale", selected_scale, "to 0")
    nodes = set(g.vertices(scale=selected_scale))              # !! but the ID of new nodes is higher than the length of the vector of the nodes
    for i in g.vertices():     # clean info for C_allocation and C_leftover in scale not used in this allocation
        if i not in nodes:
            g.node(i).C_allocation = 0
            g.node(i).C_leftover = 0


def disaggregate(g,  coarse_scale, equality=False):
    """
    Down-scale properties (C_allocation, C_excess and C_for_growth)
        from a 'coarse_scale' to the finer scale.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'coarse_scale' (int) - topological scale of the plant components between which carbon allocation is computed
        - equality (boolean) - If True, properties are downscaled, and assigned equally to all elements of a plant component
                             - If False, the assigned properties are prportional to the C_demand of each element

    :Returns:
        an MTG with attributes C_allocation, C_excess and C_for_growth estimated at the scale finer than 'coarse_scale'
        Properties are:
        'C_allocation' - carbon allocated to a component,
        'C_for_growth' - carbon available for growth in a component,
        'C_excess' - carbon allocated in excess in respect to the demand

    """
    print("Disaggregates C-allocated, C-for growth and C-excess from scale ", coarse_scale, " to scale ", coarse_scale+1)
    for j in g.vertices(scale = coarse_scale):
        node_j = g.node(j)
        for i in g.components(j):
            node_i = g.node(i)
            if equality == False:                                                            #Assign Carbon Proportional to Demand of each component
                # if  i in g.components(1531):
                #         print "in ", j
                #         print "C_allocation", node_j.C_allocation
                #         print "C_demand", node_j.C_demand
                #         print "C_for_growth", node_j.C_for_growth
                #         print "C_excess", node_j.C_excess
                #         print""
                #
                #         print "in ", i
                #         print "C_demand", node_i.C_demand
                #         print "C_allocation", node_i.C_allocation
                #         print "C_for_growth", node_i.C_for_growth
                #         print "C_excess", node_i.C_excess
                #         print "------------------------------"
                C_demand_fract_i = old_div(node_i.C_demand, float(node_j.C_demand))     # fraction of carbon demand of j due to i

                node_i.C_allocation = C_demand_fract_i * node_j.C_allocation    #update C_allocation of lower scale
                node_i.C_for_growth = C_demand_fract_i * node_j.C_for_growth   #update C_allocation of lower scale
                if node_j.C_excess > 0:
                    node_i.C_excess = C_demand_fract_i * node_j.C_excess           #Excess of Carbon is allocated proportional to the demand

                    # if  i in g.components(1531):
                    #     print "........After computation:"
                    #     print "C_demand", node_i.C_demand
                    #     print "C_allocation", node_i.C_allocation
                    #     print "C_for_growth", node_i.C_for_growth
                    #     print "C_excess", node_i.C_excess
                    #     print "_____________________________"
                else:
                    node_i.C_excess = 0   # avoids "float division by zero"
            else:                                                                            #Assign Carbon equally to all components
                node_i.C_allocation = (old_div(node_j.C_allocation, float(g.nb_components(j))))                        #update C_allocation of lower scale
                node_i.C_excess = (old_div(node_j.C_excess, float(g.nb_components(j))))                                #Excess of Carbon is allocated equally among components
                node_i.C_for_growth = (old_div(node_j.C_for_growth, float(g.nb_components(j))))                                #Excess of Carbon is allocated equally among components
    print("Disaggregation done")


