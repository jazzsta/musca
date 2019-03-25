from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import os
import numpy as np
from openalea.mtg import *
from openalea.mtg import algo
from random import random
import time

import openalea.plantgl.all as pgl
from alinea.pyratp import pyratp
from alinea.pyratp.skyvault import Skyvault
import alinea.pyratp.grid as ggr
from alinea.pyratp.vegetation import Vegetation
from alinea.pyratp.micrometeo import MicroMeteo
from alinea.pyratp.runratp import runRATP
from alinea.pyratp.RATP2VTK import RATP2VTK,RATPVOXELS2VTK
from alinea.pyratp_wralea.ExtractLight import *
import pdb

"""
All properties are assigned to the top of every internode
The root that is going to be left will be used as the root compartment, instead of an added the root compartment, as I was doing up to now
Equation for allocation adapted ()

Then verify how the model works with these changes
"""

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

    path = list1[:position]+[_gca_id]+list(reversed(list2[:position]))
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
    print("Computed barycentres and basis at scale", selected_scale)
    for i in g.vertices(scale=selected_scale):
        # print i
        barycentre(g, i)
        base(g, i)
    print("barycentres and basis computed")

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

def define_attributes_MappleT(g):
    """
    Define attributes of an MTG

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant

    :Returns:
        MTG with a attributes contained in the list "props"

    :Remark:
        - the 'fruit' attribute is added only if not already present in the MTG

    """
    print("Define attributes")
    props = "stem_demand fruit_demand leaf_demand C_demand C_supply C_allocation C_leftover C_excess volume woody_mass fruit_dry_weight init_fruit_dry_weight init_leaf_dry_weight init_woody_mass" #as in MTG from MappleT
    props = props.split()
    for prop in props:
        g.add_property(prop)
    a = "FALSE"
    for i in g.property_names():
        if i =="fruit":
            a = "TRUE"
    if a == "FALSE":
        g.add_property("fruit")
        #print "added fruit property"
    print("attributes defined", props)

def Volumes_MappleT(g):
    """
    Computes volumes of an MTG (output of MappleT) at metamer scale

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant

    :Returns:
        MTG with the volumes of individual metamers computed

    :Remark:
        - length and radius (of metamers) are in m

    """
    print("Computes Volumes at max scale")
    for vid in g.vertices(scale=max(g.scales())):
        node = g.node(vid)
        # #Formalism to compute Volumes using twice the same radius for apical metamers
        # # Assign the apical diameter for a metamer
        # if len(g.children(vid)) == 0 :                  # if vtx is apical
        #     successor = g.node(vid)                     # consider successor = node. Namely consider top and base diameter as the same for extremities.
        # else:
        #     successor = g.node(g.Successor(vid))        # if vtx is not an extremity, successor = Successor(node)
        # #Volume of trunk parts is calculated as a truncated cone (m^3)
        # node.volume = (1./3) * 3.14 * (node.radius**2 + node.radius*successor.radius + successor.radius**2) * node.length

        ##Formalism to compute Volumes using twice the same radius for basal metamer
        if not(g.parent(vid)):                            # if vtx is first metamer
            Parent = g.node(vid)                     # consider predecessor = node. Namely consider top and base diameter as the same for base node.
                                                        # however, in case of a branching, using this assumption makes the basal radius of the branch as large as the trunk at insertion point (maybe too much!!)
        else:
            Parent = g.node(g.parent(vid))        # if vtx is not an extremity, Parent = parent(node)
        #Volume of trunk parts is calculated as a truncated cone (m^3)
        node.volume = (old_div(1.,3)) * 3.14 * (Parent.radius**2 + Parent.radius*node.radius + node.radius**2 ) * node.length
    print("volumes_MappleT computed")

def proxy_root(g, init_root):
    """
    Create a root for the MTG

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'init_root' - boolean variable. If True, a root element is added to the MTG

    :Returns:
        An MTG with a root

    :Remark:
        - The root is assigned to a new component starting from the first basal component of the plant
        - its length is equal to the mean distance between the root itself and the vegetative shoots (at metamer scale)
        - its woody mass equals a fraction (equal to shoot_to_root_growth) of the biomass of vegetative shoots

        - it calls the function 'compute_geometrical_distance' to initialize the root length

    """
    print("Initialize root at max scale")
    if init_root:
        max_scale = g.max_scale()
        #first_metamer = g.component_roots_at_scale_iter(g.root, scale=max_scale).next()  #substituted to have root on the actual root
        #root_vtx = g.add_child(first_metamer, label = 'R', edge_type='+')                #substituted to have root on the actual root
        #root_node = g.node(root_vtx)                                                     #substituted to have root on the actual root

        root_vtx = g.root                                                                #inserted to have root on the actual root
        root_node = g.node(root_vtx)                                                     #inserted to have root on the actual root
        root_node.label = 'R'                                                            #inserted to have root on the actual root

        tot_dist_shoot_root = 0
        nb_shoots = 0
        tot_shoot_mass = 0
        shoot_to_root_growth = 4.5

        for vid in g.vertices(scale=max_scale):
            node = g.node(vid)
            if node.leaf_area:
                #tot_dist_shoot_root += compute_geometrical_distance(g, g.parent(root_vtx), vid)   # mean distance between shoots and the root
                #tot_dist_shoot_root += compute_geometrical_distance(g, 3, vid, g.max_scale())   # mean distance between shoots and the root
                tot_dist_shoot_root += compute_distance(g,3,vid)
                nb_shoots +=1
                tot_shoot_mass += node.woody_mass + node.leaf_dry_weight

        root_node.woody_mass = old_div(tot_shoot_mass, shoot_to_root_growth)
        #print "(total root biomass is:", root_node.woody_mass, ")"
        root_node.length = old_div(tot_dist_shoot_root, nb_shoots)
        #print "(root length is:", root_node.length, ")"
        root_node.ZZ = g.node(g.roots(scale=1)[0]).ZZ - old_div(root_node.length,2)  # spatial coordinates
        root_node.XX, root_node.YY = 0, 0
        #print "(mean root depth is:", root_node.length/2, ")"
        root_node.leaf_area, root_node.leaf_dry_weight, root_node.fruit, root_node.fruit_dry_weight, node.C_excess, node.C_supply = 0,0,0,0,0,0
    print("root initialized")


def compute_biomasses(g, GDD, init_root=False):
    """
    Computes biomasses of an MTG

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'GDD' (real) - Accumulated Growing Degree Days from full bloom (after cutoff outside the range 4.5 - 35 C)
        - 'init_root' (boolean) - If 'True', proxy_root function is called.

    :Returns:
        An MTG with a root

    :Remark:
        - dry mass of fruits is calculated starting from fresh weights
        - dry mass of old wood is calculated starting from its volume by using a fixed density (wood_density)
        - dry mass of metamers of current year shoots is calculated as follows:
            - whole shoot dry biomass is calculated by using a thermal time dependent allomatric relationship:
                shoot biomass = f(shoot length, GDD)
            - biomass of individual metamer is a fraction of whole shoot biomass and proportional to the metamer volume

            Units are:
            - Geometrical Properties of MTG
            - Volumes in m^3
            - leaf_area is in m^2
            - all biomasses in kg

        - it calls the 'proxy_root' function to initialize a root.

    """
    print("Computes biomasses at max scale")
    import math
    wood_density = 700 #kg/m3
    spec_leaf_surface  = 0.080 # kg/m^2
    # dry_to_fresh_DWratio = 0.15 # ratio to convert fresh to dry weight of Apples
    dry_to_fresh_DWratio = 1 # ratio to convert fresh to dry weight of Apples: I think input fruit are already expressed in dry mass

    for vid in g.vertices(scale=max(g.scales())):
        node = g.node(vid)
        node.C_excess, node.C_supply = 0,0
        if not node.fruit:                          # add value = 0 where there is none or is already = 0
            node.fruit = 0

        # assign a default value of dry weight = 7g to all fruits (good for DOY 180)
        if node.fruit > 0:
            node.fruit = 0.008

        node.fruit_dry_weight = node.fruit * dry_to_fresh_DWratio # convert fresh weight to dry weight
        node.leaf_dry_weight = node.leaf_area * spec_leaf_surface
        if node.label != 'R':
            if node.leaf_area == 0:                          # trunk parts (old wood)
                node.woody_mass = (node.volume * wood_density)    # compute dry mass of old wood
            else:
                shoot_length = g.node(g.complex(vid)).length
                # Shoot dry mass is calculated according to GLM model.
                # Since GLM gives result in grams, this is multiplied by g_to_Kg = 0.001 to obtain mass in kg
                # Since GLM uses lengths in cm, this is multiplied by m_to_cm = 100 to obtain length in m
                # GLM model: Shoot sry mass = f(shoot_length, GDD)
                g_to_Kg = 0.001
                m_to_cm = 100
                shoot_mass =  g_to_Kg * math.exp(-3.07340953689901 + 0.706724496955661 * math.log(shoot_length * m_to_cm ) + 5.78240214303552e-05 * GDD + 0.0200137067735007 * (shoot_length * m_to_cm ) + 9.81192475591949e-05 * math.log(shoot_length * m_to_cm ) * GDD)
                shoot_volume = g.node(g.complex(vid)).volume
                node.woody_mass = (old_div(node.volume, shoot_volume)) * shoot_mass  # mass of M equals its own proportion of volume in respect to the whole shoot, times the mass of the whole shoot
                #print "biomass", vid, node.woody_mass, "vegetative shoot", shoot_mass, "GDD", GDD, "node volume", node.volume
    proxy_root(g, init_root)

    for vid in g.vertices(scale=max(g.scales())):
        node = g.node(vid)
        node.init_fruit_dry_weight = node.fruit_dry_weight
        node.init_leaf_dry_weight = node.leaf_dry_weight
        node.init_woody_mass = node.woody_mass
    print("biomass computed")


def activities_bis_functional(g, GDD, DD, i, day, constant_supply = False):
    """
    Predicts Carbon Demands and Leaf photosynthesis of an MTG
    Assigns C_excess to supplies (optional)
    Initialize C-allocation, C-excess to 0, at max scale

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'GDD' (real) - Accumulated Growing Degree Days from full bloom (after cutoff outside the range 4.5 - 35 C)
        - 'DD' (real) - Growing Degree Days of the current day (after cutoff outside the range 4.5 - 35 C)
        - 'constant_supply' (boolean) - If True, photosynthesis is calculated with a radiative model
                                      - If False, a constant (photo_activity) per-unit-surface photosynthesis is assumed

    :Returns:
        An MTG with predicted carbon demands and photosynthesis at metamer scale

    :Remark:
        It takes values of leaf relative growth rate from a Look Up Table (LUT) 'leaf_table2'

        - Sink activities of the stems of vegetative shoots and fruits are calculated
            based on functional representations of their potential relative growth rates.
            These are the normalized derivatives of Gompertz growth curves fitted to the maximum potential growth
        - Sink activities of leaves are taken from LUT tables
        - Carbon demand = organ mass x organ specific activity;
          organ specific activity = f(GDD) x Degree Days of the current day
        - leaf supply = leaf surface x
          (if for constant_supply = TRUE:) a constant (photo_activity)
          (if for constant_supply = FALSE:) photosynthesis values obtained by radiative model
        - Units:
            length in m
            leaf_area in m^2,
            fruit mass in kg,
            DD is in C Celsius

    """
    print("Activities_bis_functional: computes Demands and Supplies. Assigns supplies to Leftover(optional). Initialize C-allocation, C-excess to 0, at scale", max(g.scales()))
    import math
    import csv

    boost = 2 # factor multiplied by carbon demand. Default for "no-effect" is "1"

    timestep = 86400  # timestep in seconds (1800 s = half an hour; 86400 s = 24 hours)

    dry_mass_to_C_mass = old_div(1,0.475) # dry matter to C mass ratio

    trunk_activity = 0.000031 * boost # g * g^-1 * C^-1
    spec_leaf_surface  = 0.080 # kg/m2

    photo_activity = 0.000648 # (C from CO2 absorption per hour) kg/h*m^2 at Amax
    costant_leaf_activity = 0.030
    # Gompertz parameters for different organs
    a_prol = 2.406
    b_prol = -2.84
    c_prol = -0.001648

    a_fruit = 63.703
    b_fruit = -4.896
    c_fruit = -0.00104

    # RGR is equal to the derivative of the fitted Gompertz: a*b*c*exp(b*(exp(c*x))+c*x)
    RGR_prol = boost * (a_prol * b_prol * c_prol * math.exp(b_prol * (math.exp(c_prol * GDD)) + c_prol * GDD)) / (a_prol * math.exp(b_prol * math.exp(c_prol * GDD)))
    #RGR_epic = (a_epic * b_epic * c_epic * math.exp(b_epic * (math.exp(c_epic * GDD)) + c_epic * GDD)) / (a_epic * math.exp(b_epic * math.exp(c_epic * GDD)))
    RGR_fruit = boost * (a_fruit * b_fruit * c_fruit * math.exp(b_fruit * (math.exp(c_fruit * GDD))+ c_fruit * GDD)) / (a_fruit * math.exp(b_fruit * math.exp(c_fruit * GDD)))

    leaf_table2 = csv.reader (open ("D:\PROJECTS\L-PEACH_adaptation\Maximum Growth Potential Experiment\Rilievi Analisi\LeafRGR.csv" , 'r'))
    Col1_leaf = "ColumnName1"
    Col2_leaf = "ColumnName2"
    leaf_RGR = {Col1_leaf:[], Col2_leaf:[]}
    for row in leaf_table2:
        leaf_RGR[Col1_leaf].append(float(row[0]))
        leaf_RGR[Col2_leaf].append(int(row[1]))

    for j in range(0,len(leaf_RGR[Col2_leaf])):
        if leaf_RGR[Col2_leaf][j] == int(GDD):
            leaf_activity = leaf_RGR[Col1_leaf][j] * boost
    # Initialize C.Demand and Supplies
    for vid in g.vertices(scale=max(g.scales())):
        node = g.node(vid)
        node.C_allocation = 0
        # node.C_supply += node.C_excess               # Use C-excess as supply in next day
        node.C_excess = 0

        # leaf demand
        if leaf_activity > 0:
            node.leaf_demand = leaf_activity * node.leaf_area * spec_leaf_surface *  DD
        else:
            node.leaf_demand = 0
        # fruit demand
        node.fruit_demand = RGR_fruit * node.fruit_dry_weight * DD

        # trunk (old wood)
        if node.leaf_area == 0 and node.label != 'R':
            node.stem_demand = trunk_activity * node.woody_mass * DD
            node.tot_C_supply = 0            #variable used to compute total C_supply of individual leaf over the whole simulation

        # root
        elif node.label == 'R':
            node.stem_demand = RGR_prol * node.woody_mass   * DD
            node.tot_C_supply = 0            #to compute total C_supply of metamer

        # current year shoots
        else:
            #if g.complex(vid).length > 35:               ## on a fixed structure: shoots longer than 35cm follow epicormics RGR curves
            #    RGR_shoot = RGR_epic
            #else:
            #    RGR_shoot = RGR_prol                     ## on a fixed structure: shoots shorter than 35cm follow proleptic RGR curves
            RGR_shoot = RGR_prol
            node.stem_demand = RGR_shoot * node.woody_mass   * DD

            # Carbon Supply
            # If a constant photosynthetic activity is preferred to the use of a radiative model
            if constant_supply:
                node.C_supply = photo_activity * node.leaf_area
            else:
            # Conversion of micromol CO2 s-1 to g organic dry matter (C to dry matter)
            # g organic dry matter = micromol CO2 * CO2 molar weight in grams * g/microg * Carbon to Oxygen mass ratio in CO2 * mass of dry matter per mass of C
            # CO2 molar weight in grams = 44; g/microg = 10exp-6; u.m.a. ratio between C and CO2 = 12/44; dry matter to C mass ratio = 1/0.475; g in Kg = 0.001; timestep in seconds
                node.Photo_unit_surf = old_div((node.PhotResu), (node.LeafAreaResu))
                node.C_supply = (node.Photo_unit_surf * node.leaf_area * 44 * 0.000001 * (old_div(12.,44)) * dry_mass_to_C_mass * 0.001 *  timestep)
    #           node.C_supply = [node.C_excess] + (node.Photo_unit_surf * node.leaf_area * molar_weight_CO2 * umol_to_mol * molCO2_to_molC * C_to_dry_matter * g_to_kg * timestep)
                if i==day:                           #to compute total C_supply of metamer
                    node.tot_C_supply = 0            #to compute total C_supply of metamer
                node.tot_C_supply += node.C_supply   #to compute total C_supply of metamer
        node.C_leftover = node.C_supply
        node.C_demand = node.fruit_demand + node.stem_demand + node.leaf_demand
    #print "(trunk activity is", trunk_activity, "shoot activity is", RGR_shoot, "fruit activity is", RGR_fruit, "leaf activity is", leaf_activity, ")"
    print("activities computed")

def aggregate_props(g, coarse_scale, properties_names):
    """
    Aggregates a property of an MTG up to the 'coarse' scale
    starting from the properties at the finer scale

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'coarse_scale' (int) - Scale to which aggregating the properties
        - 'properties_names' (list) - list of properties to aggregate

    :Returns:
        An MTG with the properties in the 'properties_names' list aggregated at the coarse scale

    """
    print("Aggregates ", properties_names.split(), "to coarse scale")
    properties_names = properties_names.split()
    for property_name in properties_names:

        prop = g.property(property_name)
        GUs = g.vertices(scale=coarse_scale)
        for i in GUs:      #assegna a ogni elemento della scala piu grezza la somma dei suoi componenti alla scala piu fine, per ogni proprieta
            #print property_name, "for vertex", i
            try:
                prop[i] = sum( prop[vid] for vid in g.components(i) )
                #print "aggregated"
            except:
                #print "aggregation not possible on vtx"
                pass

def aggregate_multi(g, coarse_scale):
    """
    Aggregates a given set of properties of an MTG at the "coarse" scale
        starting from the properties at the finer scale

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'coarse_scale' (int) - Scale to which aggregating the properties

    :Returns:
        An MTG with the properties in the 'PropS' list aggregated at the coarse scale

    :Remark:
    - it calls the function 'aggregate_props' over the list of properties 'PropS'
    - the properties in 'PropS' are:
        C_allocation C_demand C_excess C_leftover C_supply
        fruit_demand fruit_dry_weight leaf_area LeafAreaResu leaf_demand leaf_dry_weight length
        init_fruit_dry_weight init_leaf_dry_weight init_woody_mass observation stem_demand woody_mass

    """
    PropS = "C_demand stem_demand fruit_demand leaf_demand C_excess C_supply C_allocation C_leftover fruit_dry_weight leaf_area observation leaf_dry_weight woody_mass LeafAreaResu init_fruit_dry_weight init_leaf_dry_weight init_woody_mass" #volume
    aggregate_props(g, coarse_scale, PropS)

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

def distances(g, selected_scale):
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
    Dists = DistanceMatrix(g.vertices(scale=selected_scale))
    nb_vertices = len(g.vertices(scale=selected_scale))
    nb_vertices_done = 0
    score_init = 0
    for Leaf in g.vertices(scale=selected_scale):
        for Demand in g.vertices(scale=selected_scale):
            # Dists[Leaf,Demand] = compute_geometrical_distance(g,Leaf,Demand, selected_scale)                # substituted to be able to have distance from component to same component = 0
            Dists[Leaf,Demand] = compute_distance(g,Leaf,Demand)

        nb_vertices_done += 1
        score = sum((old_div((nb_vertices_done*100),nb_vertices)) >= np.array([ 10,  20,  30,  40, 50,  60, 70, 80, 90, 100]))
        if score > score_init:
            print("computed",(old_div((nb_vertices_done*100),nb_vertices)), "% ofvertices (", nb_vertices_done,")")
            score_init = score
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
        if g.node(Leaf).leaf_area > 0:
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

def check(g,v=1531):
    node = g.node(v)
    sum_C_allocated = 0
    sum_C_demand = 0
    sum_C_excess = 0
    sum_C_for_growth = 0
    for i in g.components(v):
        sum_C_allocated += g.node(i).C_allocation
        sum_C_demand += g.node(i).C_demand
        sum_C_excess += g.node(i).C_excess
        sum_C_for_growth += g.node(i).C_for_growth
    print("At fine scale:")
    print("sum_C_allocated = ", sum_C_allocated)
    print("sum_C_demand = ", sum_C_demand)
    print("sum_C_excess = ", sum_C_excess)
    print("sum_C_for_growth = ", sum_C_for_growth)
    print("")
    print("At coarse scale:")
    print("C_allocated ", node.C_allocation)
    print("C_demand ", node.C_demand)
    print("C_excess ", node.C_excess)
    print("C_for_growth ", node.C_for_growth)

def growth(g, DD):
    """
    Computes carbon allocated, updated dry weight and relative growth rate of the organ types present
        on an individual plant component (at the finest scale) according to the C_allocated to it.

        - C_allocated is distribuited to the woody, leafy and
            fruit parts proportionally to their C_demands
        - Carbon allocated to a component in excess in respect to the demand (C_excess) is assigned to the supply of the component

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'DD' (real) - Growing Degree Days of the current day (after cutoff outside the range 4.5 - 35 C)

    :Returns:
        an MTG with attributes estimated at the finest scale for the current time step:
        - 'allo_stem' - carbon allocated to the internode of a component,
        - 'allo_leaf' - carbon available to the leaves of a component,
        - 'allo_fruit' - carbon allocated to the fruit of a component,

        - 'wood_RGR' - relative growth rate of the internode of a component,
        - 'leaf_RGR' - relative growth rate of the leaves of a component,
        - 'fruit_RGR' - relative growth rate of the fruit of a component,


        - 'woody_mass' - updated dry weight of the internode of a component,
        - 'leaf_area' - updated dry weight of the leaves of a component,
        - 'fruit_dry_weight' - updated dry weight of the fruit of a component,

        - 'C_supply' - updated carbon supply of the component

    """
    print("Computes growth: C going to wood, leaf and fruit; RGRs; growth of the current day, all at scale", g.max_scale())
    for vid in g.vertices(scale = g.max_scale()):
        # Carbon available for growth = C_for_growth * organ_demand / node.demand
        # assign to the component (e.g. to the fruit) a fraction of the carbon allocated to its node proportional to its own weight
        g.node(vid).allo_stem = g.node(vid).C_for_growth * (old_div(g.node(vid).stem_demand, g.node(vid).C_demand)) # to wood
        #if vid in g.components(1467):
            # print "in vid ", vid
            # print "C_for_growth ", g.node(vid).C_for_growth
            # print "allo_stem ", g.node(vid).allo_stem
        g.node(vid).allo_leaf = g.node(vid).C_for_growth * (old_div(g.node(vid).leaf_demand, g.node(vid).C_demand)) # leaf
        g.node(vid).allo_fruit= g.node(vid).C_for_growth * (old_div(g.node(vid).fruit_demand, g.node(vid).C_demand)) # fruit

        # RGR
        g.node(vid).wood_RGR = (old_div((old_div(g.node(vid).allo_stem, g.node(vid).woody_mass)), DD))*1000
        g.node(vid).leaf_RGR = 0 if g.node(vid).leaf_area == 0 else (old_div((old_div(g.node(vid).allo_leaf, g.node(vid).leaf_dry_weight)), DD))*1000
        g.node(vid).fruit_RGR = 0 if g.node(vid).fruit == 0 else (old_div((old_div(g.node(vid).allo_fruit, g.node(vid).fruit_dry_weight)), DD))*1000

        # Growth
        g.node(vid).woody_mass += g.node(vid).allo_stem
        g.node(vid).leaf_area = g.node(vid).leaf_area #+= g.node(vid).allo_leaf
        g.node(vid).fruit_dry_weight += g.node(vid).allo_fruit
    print("growth computed")


meteo_values = None

def meteo_init():
    """
    Extract Growing Degree Days (GDD) of single days and Accumulated GDD from full bloom from a meteorological file

    :Returns:
         - the variable 'meteo_values' containing couples of values 'cumDD' (Accumulated GDD), 'DD'(GDD) for every DOY.

    """

    global meteo_values
    import csv
    T = open("D:\PROJECTS\L-PEACH_adaptation\Maximum Growth Potential Experiment\Rilievi Analisi\DegreeDays.csv" , 'r')
    csv_T = csv.reader(T)
    meteo_values = dict()
    for row in csv_T:
        if (row[2])!= "cumDD":          #except for heading
            day = int(row[0])
            cumDD = float(row[2]) #int(float(row[2]))
            DDraw = float(row[1]) #int(float(row[1]))
            DD = DDraw - 4.5 #4.5 is cutoff temperatures
            meteo_values[day] = (cumDD, DD)

def meteo(day):
    """
    Extract Growing Degree Days (GDD) of single days and Accumulated GDD from full bloom from a meteorological file

    :Parameters:
        - 'day' - the current timestep

    :Returns:
         - the couple of values 'cumDD' (Accumulated GDD) and 'DD'(GDD) of the current timestep ('day').

    """
    if meteo_values is None: meteo_init()
    return meteo_values[day]

def leafy_shoot_axis_and_root(g):
    """
    Set topological boundaries of a new scale in the place where at the finest scale there is:
    - a root component (labelled 'R')
    - main trunk components whose parent generated a branch
    - vegetative shoots

    :Parameters:
        - 'v' - a vertex

    :Returns:
         - a boolean for the vertex 'v', indicating if that vertex correspond to a boundary.

    :Remark:
        - it is used in the 'insert_scale' function.

    """
    def partition(v):
        if v == g.vertices(scale=g.max_scale())[0] or g.label(v)=='R' or (g.node(v).leaf_area > 0 and g.node(g.parent(v)).leaf_area == 0) or (v in g.Trunk(3) and g.Sons(g.parent(v), EdgeType= '+')):
            return True
        else:
            return False
    return partition

def branches_and_root(g):
    """
    Set topological boundaries of a new scale in the place where at the finest scale there is:
    - a root component (labelled 'R')
    - main branches originating from the trunk

    :Parameters:
        - 'v' - a vertex

    :Returns:
         - a boolean for the vertex 'v', indicating if that vertex correspond to a boundary.

    :Remark:
        - it is used in the 'insert_scale' function.

    """
    def partition(v):
        if v == g.vertices(scale=g.max_scale())[0] or g.label(v)=='R' or (v not in g.Trunk(3) and g.parent(v) in g.Trunk(3) and g.Sons(g.parent(v), EdgeType= '+')):
            return True
        else:
            return False
    return partition

def trunk_branches_shoots(g):
    """
    Set topological boundaries of a new scale in the place where, at the finest scale, there is:
    - a root component (labelled 'R')
    - the trunk starts
    - main branches originating from the trunk
    - a leafy shoot

    :Parameters:
        - 'v' - a vertex

    :Returns:
         - a boolean for the vertex 'v', indicating if that vertex correspond to a boundary.

    :Remark:
        - it is used in the 'insert_scale' function.

    """
    def partition(v):
        if v == g.vertices(scale=g.max_scale())[0] or g.label(v)=='R' or (v not in g.Trunk(3) and g.parent(v) in g.Trunk(3) and g.Sons(g.parent(v), EdgeType= '+')) or (g.node(v).leaf_area > 0 and g.node(g.parent(v)).leaf_area == 0):
            return True
        else:
            return False
    return partition


def metamer(g):
    """
    Set topological boundaries of a new scale in the place where at the finest scale there is:
    - a metamer

    :Parameters:
        - 'v' - a vertex

    :Returns:
         - a boolean for the vertex 'v', indicating if that vertex correspond to a boundary.

    :Remark:
        - it is used in the 'insert_scale' function.

    """
    def partition(v):
        if v == g.vertices(scale=g.max_scale()):
            return True
        else:
            return False
    return partition

def growth_unit(g):
    """
    Set topological boundaries of a new scale in the place where at the finest scale there is:
    - a metamer

    :Parameters:
        - 'v' - a vertex

    :Returns:
         - a boolean for the vertex 'v', indicating if that vertex correspond to a boundary.

    :Remark:
        - it is used in the 'insert_scale' function.
    """
    def partition(v):
        if v == g.vertices(scale=2):
            return True
        else:
            return False
    return partition

def branches_interbranches_and_root(g):
    """
    Set topological boundaries of a new scale in the place where at the finest scale there is:
    - a root component (labelled 'R')
    - main branches originating from the trunk
    - sections of trunk comprised between two branches originating from the trunk

    :Parameters:
        - 'v' - a vertex

    :Returns:
         - a boolean for the vertex 'v', indicating if that vertex correspond to a boundary.

    :Remark:
        - it is used in the 'insert_scale' function.

    """
    def partition(v):
        if v == g.vertices(scale=g.max_scale())[0] or g.label(v)=='R' or (v in g.Trunk(3) and g.Sons(v, EdgeType= '+')) or (v not in g.Trunk(3) and g.parent(v) in g.Trunk(3) and g.Sons(g.parent(v), EdgeType= '+')):
            return True
        else:
            return False
    return partition

def leafy_shoot_and_root(g):
    """
    Set topological boundaries of a new scale in the place where at the finest scale there is:
    - a root component (labelled 'R')
    - leafy shoots

    :Parameters:
        - 'v' - a vertex

    :Returns:
         - a boolean for the vertex 'v', indicating if that vertex correspond to a boundary.

    :Remark:
        - it is used in the 'insert_scale' function.

    """
    def partition(v):
        if v == g.vertices(scale=g.max_scale())[0] or g.label(v)=='R' or (g.node(v).leaf_area > 0 and g.node(g.parent(v)).leaf_area == 0):
            #print v
            return True
        else:
            return False
    return partition

def quotient_leaf_and_ramif(g):
    """
    Set topological boundaries of a new scale in the place where at the finest scale there is:
    - a root component (labelled 'R')
    - a branching connection '+'
    - a leaf

    :Parameters:
        - 'v' - a vertex

    :Returns:
         - a boolean for the vertex 'v', indicating if that vertex correspond to a boundary.

    :Remark:
        - it is used in the 'insert_scale' function.

    """
    def partition(v):
        if v == g.vertices(scale=g.max_scale())[0] or g.label(v)=='R' or g.edge_type(v) != '<':
            return True
        else:
            if g.node(v).leaf_area > 0:   #leaf_state == "growing" :
                return True
            else :
                return False
    return partition

def fruiting_unit_and_root(g):
    """
    Set topological boundaries of a new scale in the place where at the finest scale there is:
    - a root component (labelled 'R')
    - the property 'new_scale' different than zero

    :Parameters:
        - 'v' - a vertex

    :Returns:
         - a boolean for the vertex 'v', indicating if that vertex correspond to a boundary.

    :Remark:
        - it is used in the 'insert_scale' function.

    """
    def partition(v):
        if v == g.vertices(scale=g.max_scale())[0] or g.label(v)=='R' or (g.node(v).new_scale):
            return True
        else:
            return False
    return partition

def replace_scale_fun(g, partition):
    """
    Replaces a scale of representation denoted '2' with a new one.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'partition' - the function containing the rules used to set topological boundaries of the new scale

    :Returns:
        - 'g' (MTG) - a multi scale tree graph with a new scale

    :Remark:
        If the 'partition' is 'fruiting_unit_and_root' it assigns positive values to
        a new attribute 'new_scale' in selected vertices. This is used by the 'fruiting_unit_and_root'
        'partition' to locate the position where fruiting units should start.

    """
    print("Replace scale")
    if partition == fruiting_unit_and_root:
        partition = partition(g)
        for v in g.vertices(scale=max(g.scales())):
            if g.Sons(g.complex(v)):
                if v == g.components(g.complex(v))[0] and g.node(g.Sons(g.complex(v))[0]).leaf_area and not g.node(v).leaf_area:
                    #print v
                    g.node(v).new_scale = 1
    else:
        partition = partition(g)

    XX = g.node(g.roots(scale=2)[0]).XX  # store basal root coordinates
    YY = g.node(g.roots(scale=2)[0]).YY
    ZZ = g.node(g.roots(scale=2)[0]).ZZ
    g.remove_scale(2)
    g = g.insert_scale(inf_scale=g.max_scale(), partition = partition, default_label='LS')
    g.node(g.roots(scale=2)[0]).XX = XX   # re-store basal root coordinates
    g.node(g.roots(scale=2)[0]).YY = YY
    g.node(g.roots(scale=2)[0]).ZZ = ZZ
    print("scale replaced")

def main_init(g, GDD, selected_scale, init_root, replace_scale, partition, coarse_Scale=2):
    """
    Estimates volumes and dry biomasses of individual plant components
    at the finest scale of an MTG in the format produced by the MAppleT model,
    and optionally replace scale '2'

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'GDD' (real) - Accumulated Growing Degree Days from full bloom (after cutoff outside the range 4.5 - 35 C)
        - 'init_root' (boolean) - If 'True', proxy_root function is called.
        - 'replace_scale' - If 'True', scale '2' is sibstituted by calling the 'replace_scale_fun' function.
        - 'partition' - the function containing the rules used to set topological boundaries of the new scale

    :Returns:
        - 'g' (MTG) - a multi scale tree graph with volumes, biomasses and optionally a new scale.

    """
    define_attributes_MappleT(g)
    Volumes_MappleT(g)                                       # Compute Volumes from mtg from MappleT
    aggregate_props(g, coarse_Scale, 'volume length leaf_area')        # Aggregate Volumes and Lengths from M to coarse scale
    bary_base(g, g.max_scale())
    vertices_semilength(g,g.max_scale())
    compute_biomasses(g, GDD, init_root)                     # Compute Biomasses from Volumes and GDD dependent allometric relationships

    if replace_scale == True:
#        pdb.set_trace()
        replace_scale_fun(g, partition)
        for i in g.vertices(scale=g.max_scale()):
            if g.node(i).label=='R': g.node(g.complex(i)).label = 'R'
    bary_base(g, selected_scale)
    vertices_semilength(g,selected_scale)

def radiative_model (g, DayToExtract):
    """
    Computes photosythesis performed by individual leaves starting from direct
    and diffused light, using the RATP radiative model.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'i' - the currently simulated timestep

    :Returns:
        - 'g' - an MTG with 'PhotResu' and 'LeafAreaResu', the properties containing photosynthesis values

    :Remark:
        It uses as input files:
        - 'gridfile' for the discretization of space into voxels
        - 'vfnfile' for the vegetation parameters used in RATP
        - 'skyfile' the discretization of sky for incoming radiation
        - 'name' the meteofile

        Photosynthesis is predicted for 30 minutes timebands and then averaged over the whole day.
    """
    if DayToExtract: # == day:
#		os.chdir('D:/Softwares/tutorial/src/alinea/tutorial')
        print("Radiative model")
        #### Initialize all input files
        grdfile='grid3Da_2004.grd'              #Grid
        #print 2
        ###!!!
        #MTGfile='0.mtg'             #3D plant
        #MTGfile='output_MappleT_low_fruit_load.mtg'             #3D plant
        #MTGfile='output_MappleT_medium_fruit_load.mtg'             #3D plant

        vfnfile='vegetationa_2004.vfn'          #List of plant functioning files
        skyfile='skyvaultsoc.skv'               #Turtle
        #metfile='hh_meteofile_RATP_Caldaro2010_DOY_210_279.mto' #Meteo
        name='D:/Softwares/tutorial/src/alinea/tutorial/meteofiles/hh_meteofile_RATP_Caldaro2010_' #Meteo
        metfile = name + str(DayToExtract) + '.mto'
        #print 3
        OutputfileVoxels='C:\outputVoxels.vtk'  #Output file at voxel scale
        OutputfileLeaves='C:\outputLeaves.vtk'  #Output file at leaf scale
        #print 4
        #g = MTG(MTGfile)

        metamer_number = 0
        for vid in g.vertices(scale=3):
                node = g.node(vid)
                if node.leaf_area !=0:
                        metamer_number += 1
        #print 5
        X = np.zeros(metamer_number)
        Y = np.zeros(metamer_number)
        Z = np.zeros(metamer_number)
        leaf_area  = np.zeros(metamer_number)
        nitrogen = np.zeros(metamer_number)
        entity = np.zeros(metamer_number)
        #print 6

        acc = 0
        for vid in g.vertices(scale=3):
                node = g.node(vid)
                if node.leaf_area !=0:
                        node.id = acc
                        X[acc] = node.XX
                        Y[acc] = node.YY
                        Z[acc] = node.ZZ
                        leaf_area[acc] = node.leaf_area
                        nitrogen[acc] = 2
                        acc += 1
                else :
                        node.id = 99999999
        #print 7
        X = X - 1.2*np.min(X)
        Y = Y - 1.2*np.min(Y)
        #print 8
        #### Initialize the grid from a grid file
        grid = ggr.Grid.read(grdfile)
        #print "X is", X
        #print "Y is", Y
        #print "Z is", Z
        #print "leaf area is", leaf_area
        #print "nitrogen is", nitrogen
        ####fill the grid
        grid,map = ggr.Grid.fill(entity, X, Y, Z, leaf_area, nitrogen, grid)
        #print 9
        ###Write to VGX - Check the Voxels
        ggr.gridToVGX(grid,os.getcwd(),"\gridVGX.vgx")
        #print 10
        #### Vegetation Type
        vegetation = Vegetation.read(vfnfile)
        ##Acces aux parametres physio
        #vegeType = pyratp.vegetation_types
        #vegeType.mu[0] = 0.6
        #print 11
        #### Sky
        sky = Skyvault.read(skyfile)  ###
        #print 12
        #### Meteo
        met = MicroMeteo.read(metfile)
        #print 13
        #### Run RATP
        res = runRATP.DoAll()
        #print 14
        ### get the STAR_sky
        hemi3D= pyratp.hemi_interception
        #print 15
        vegeType = pyratp.vegetation_types
        #vegeType.mu[0] = 0.6
        #print 16
        ##### change the mu to test its impact on STAR
        #mutarget = np.arange(0.5,1,0.02)
        #for i in mutarget:
        #	vegeType.mu = i*np.ones(grid.nent)
        #	resIrrad = runRATP.DoIrradiation()
        #	Tree_STARsky = hemi3D.starsky_canopy
        #	print i
        #	print Tree_STARsky
        #	print '-----------'
        #print 17
        ##[Iteration,day,hour,VoxelId,ShadedPAR,SunlitPAR,ShadedArea,SunlitArea]= res[0].T
        ##colnames = ['Entity','Iteration','day','hour','AirTemperature','VoxelId','ShadedTemp','SunlitTemp','VoxelTemp','STARDirect','STARSky','ShadedPhoto','SunlitPhoto','ShadedTranspi','SunlitTranspi','ShadedArea','SunlitArea','ShadedGs','SunlitGs','VoxelNitrogen']
        colnames = ['Iteration','day','hour','AirTemperature','VoxelId','ShadedTemp','SunlitTemp','VoxelTemp','STARDirect','STARSky','ShadedPhoto','SunlitPhoto','ShadedTranspi','SunlitTranspi','ShadedArea','SunlitArea','ShadedGs','SunlitGs','VoxelNitrogen']
        ##TTTT = res[0].T[1:]
        ##TT = TTTT.T
    #print 18
    ############
    ### Extract a variable at a given time (day and hour)
    ##!!!!!!!!!! cycle to average Photo over the hours of the day
    # Columns: ShadedPhoto=10, SunlitPhoto=11, ShadedArea=14, SunlitArea=15
    indexvar = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
    #DayToExtract =   210 #184
    HourToExtract = 13

    ShaPhot = [0] * len(ExtractLight(map, res[0], DayToExtract, HourToExtract, 10))  #At leaf scale)
    SunPhot = [0] * len(ShaPhot)
    PhotResu = [0] * len(ShaPhot)
    LeafAreaResu = [0] * len(ShaPhot)
    #
    #print 19
    nb_time_bands = 48
    #nb_time_bands = 1 #for 1 hour at noon
    for i in range(0, nb_time_bands):
    #for i in range(24,25): #for 1 hour at noon
            # save semi-hourly values in a temporary list,
            TempShaPhot  = (ExtractLight(map, res[0], DayToExtract, old_div(i,2.0), 10))
            TempSunPhot  = (ExtractLight(map, res[0], DayToExtract, old_div(i,2.0), 11))

            TempSunAreaResu = ExtractLight(map, res[0], DayToExtract, old_div(i,2.0), 15)
            TempShaAreaResu = ExtractLight(map, res[0], DayToExtract, old_div(i,2.0), 14)

            # than loop through the temporary list in order to update the lists PhotResu and LeafAreaResu
            for j in range(len(ShaPhot)):
                    PhotResu[j] += old_div(TempShaPhot[j],nb_time_bands) + old_div(TempSunPhot[j],nb_time_bands)
                    LeafAreaResu[j] += old_div(TempSunAreaResu[j],nb_time_bands) + old_div(TempShaAreaResu [j],nb_time_bands)
    ##!!!!!!!!!
    ###############
    #print 20
    leaf_area = 0
    for vid in g.vertices(scale=3):
            node = g.node(vid)
            leaf_area += node.leaf_area
    #print leaf_area

    ColToExtract = 11  # what for?
    VarResuVoxels = ExtractVoxels(res[0], DayToExtract, HourToExtract,ColToExtract) #At voxel scale
    #print 21
    g.add_property('PhotResu')
    g.add_property('LeafAreaResu')
    for vid in g.vertices(scale = 3) :
            node = g.node(vid)
            if (node.id) != 99999999 :
                    setattr(node,'PhotResu',PhotResu[node.id]) #new
                    setattr(node,'LeafAreaResu',LeafAreaResu[node.id]) #new
                    last_id = vid

    variable_liste = []
    #print 22
    node = g.node(last_id)
    liste = node.properties()
    #print 23
    for key in liste :
        if ((key != 'edge_type') and (key != '_line') and (key != 'index') and (key != 'label') and (key!= 'id')) :
            if isinstance(liste[key],int) :
                type_mtg = 'INT'
            if isinstance(liste[key],float):
                type_mtg = 'REAL'
            if isinstance(liste[key],str):
                type_mtg = 'ALPHA'
            tuple_ = (key,type_mtg)
#	        print tuple_
            variable_liste.append(tuple_)
    #print 24
    l= open('output_MTG',"w")
    l.write(io.write_mtg(g,variable_liste))
    l.close()
    print("radiative model DONE")

def Final_outputs(g, init_GDD, GDD, pre, suff, N, i, day, constant_supply, output_length, beginning_time, time_Dists, selected_scale):
    """
    Computes spatially explicit and implicit information about the whole simulated period
        for the different organ types present on each plant component at the finest scale.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'DD' (real) - Growing Degree Days of the current day (after cutoff outside the range 4.5 - 35 C)
        - 'nn' - path and filename generated by the 'output_filename' function
        - 'N' - number of simulated timesteps
        - 'i' - the currently simulated timestep
        - 'day' - first simulated timestep (first day of simulation if the unit for timestep is days)
        - 'constant_supply' - boolean True if the radiative model is replaced by an arbitrary constant per leaf value of carbon assimilation

    :Returns:
        an MTG with attributes estimated at the scale finest scale.

    :Remark:
        - it calls the '...' functions to generate outputs.
    """
    # Compute
    print("Growth_and_RGR")
    Growth_and_RGR(g, init_GDD, GDD)
    print("tree_axil")
    tree_axil(g)
    print("Spatial_Info")
    Spatial_Info(g, pre, suff, constant_supply, output_length)
    print("Growth_fruit")
    Growth_fruit, Growth_wood = compartmental_budget(g,pre, suff)

    elapsed_time = time_run(g, beginning_time, time_Dists, pre, suff, selected_scale, output_length)
    print("Elapsed time", elapsed_time)

def Growth_and_RGR(g, init_GDD, GDD):
    print("Compute growth_and_RGR for the whole simulated period")
    delta_GDD = GDD - init_GDD
    for vid in g.vertices(scale=max(g.scales())):
        node = g.node(vid)
        node.growth_fruit = node.fruit_dry_weight - node.init_fruit_dry_weight
        node.growth_leaf = node.leaf_dry_weight - node.init_leaf_dry_weight
        node.growth_wood = node.woody_mass - node.init_woody_mass

        node.RGR_woody = old_div((old_div(node.growth_wood, node.init_woody_mass)),delta_GDD)

        if node.fruit_dry_weight:
            node.RGR_fruit = old_div((old_div(node.growth_fruit, node.init_fruit_dry_weight)),delta_GDD)
        else:
            node.RGR_fruit = 0

        if node.init_leaf_dry_weight:
            node.RGR_leaf = old_div((old_div(node.growth_leaf, node.init_leaf_dry_weight)),delta_GDD)
        else:
            node.RGR_leaf = 0
    print("growth and RGR computed")

def time_run(g, beginning_time, time_Dists, pre, suff, selected_scale, output_length):
    """
    Counts and writes in a file the time taken to compute the distance matrix and for the rest of the simulation.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'beginning_time' - time at the beginning of the simulation
        - 'time_Dists' - time spent to compute distance matrix
        - 'pre' - prefix to be added to the output filename
        - 'suff' - suffix to be added to the output filename
        - 'selected_scale' - scale at which the simulation was run
        - 'output_length' - if not "synthetic" the prefix "neighbor" is added to the output filename

    :Returns:
        a file containing the time taken to compute the distance matrix (time_Dists) and for the rest of the simulation (time_non_Dist).

    :Remark:
        - it calls the '...' functions to generate outputs.
        - it takes into account the cases when the distance matrix was previously calculated for another friction parameter or scale but on the same tree structure
    """
    time_non_Dist = time.time()- beginning_time

    nb_components = len( g.vertices( scale = selected_scale ) )

    nb_leaves = 0
    for i in g.vertices( scale = selected_scale ):
        if g.node(i).leaf_area:
            nb_leaves = nb_leaves + 1

    if output_length == "synthetic":
        nn = pre + "\\time_" + suff
    else:
        nn = pre + "\\neighborhood_time_" + suff
    print("time not for distance matrix" + str(time_non_Dist))
    print("nb_components", nb_components)
    print("nb_supplies", nb_leaves)
    n = open (nn,"w")
    n.write("time_for_Dists"+'\t'+"time_non_Dist"+'\t'+ "nb_components"+'\t'+"supplies"+'\n')
    n.write( str(time_Dists) +'\t'+ str(time_non_Dist) +'\t'+ str(nb_components) +'\t'+ str(nb_leaves) +'\n')
    n.close()
    return(time_non_Dist)

def output_filename(MTGfile, day, N, h, scale_name, init_root):
    """
    Generates a path and filename including name of MTG, starting day of simulation,
        number of simulated timesteps and friction parameter.

    :Parameters:
        - 'MTGfile' (MTG) - name of the multi scale tree graph of the plant
        - 'day' - first simulated timestep (first day of simulation if the unit for timestep is days)
        - 'N' - number of simulated timesteps
        - 'h' (real) - friction parameter

    :Returns:
        - 'nn' - path and filename of the output file

    """
    # filepath = "D:\\PROJECTS\\L-Fuji in OpenAlea\\Training MTG\\C_model_Outputs\\Output_mtgs_3cultivars"
    # filepath = "D:\\PROJECTS\\L-Fuji in OpenAlea\\Training MTG\\C_model_Outputs\\Output_mtgs_1cultivars"
    filepath = "D:\\PROJECTS\\L-Fuji in OpenAlea\\Training MTG\\C_model_Outputs\\Output_mtgs_1cultivar_manyh"
    # if selected_scale == 3:
    #     scalename = "M"
    # elif replace_scale == False:
    #     scalename = "GU"
    # elif partition == fruiting_unit_and_root:
    #     scalename = "FU"
    # elif partition == branches_interbranches_and_root:
    #     scalename = "BR"
    if init_root:
        root_declaration = "rootY"
    else:
        root_declaration = "rootN"
    pre, suff = filepath, "_" + MTGfile + "_day" + str(day) + "_for"+ str(N) + "_H" + str(h) + "_" + scale_name + "_" + root_declaration + ".txt"
    return (pre, suff)

def Spatial_Info(g, pre, suff, constant_supply, output_length, selected_scale = 3):
    """
    Generates a  file with spatially explicit information about the simulated mtg, including:
        'vid' - identifier of plant part at fine scale,
        'Complex_id' - identifier of complex of 'vid',
        'organ_type' - type of organ of 'vid',
        'X', 'Y', 'Z' - basal spatial coordinates of 'vid'
        'Photo_unit_surf' - assimilated carbon per unit surface
        'init_fruit_DW' - initial fruit weight on 'vid'
        'growth_fruit' - absolute growth in fruit dry weight present on 'vid'
        'RGR_fruit' - relative growth rate of fruit present on 'vid'

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'nn' - path and filename generated by the 'output_filename' function and passed to the 'output_2' function for generating an output
        - 'constant_supply' (boolean) - passed to the 'output_2' function for generating an output
        - 'selected_scale' - scale at which producing the output

    :Returns:
         - a file of output 'nn'.

    """
    if output_length == "synthetic":
        print("")
        print("Start writing synthetic output file")
        nn = pre + "\\spatial_" + suff
        f = open (nn,"w")
        f.write("Vid"+'\t'+ "Complex_id"+'\t'+"organ_type"+'\t'+ "X"+'\t'+ "Y"+'\t'+ "Z"+'\t'+ "Photo_unit_surf"+'\t'+
                "init_fruit_DW"+'\t'+ "growth_fruit"+'\t'+ "RGR_fruit"+'\t'+ "total_photo"+'\t'+
                "Barycenter_XX"+'\t'+ "Barycenter_YY"+'\t'+ "Barycenter_ZZ"+'\t'+'\n')

        first_metamer = next(g.component_roots_at_scale_iter(g.root, scale=selected_scale))
        list_trunk = g.Trunk(first_metamer, Scale = selected_scale)

        print("Computes properties in neighborhood")
        for i in g.vertices(scale = selected_scale):
            node = g.node(i)
            XX = node.XX
            YY = node.YY
            ZZ = node.ZZ
            node.complex_id = g.complex(i)     # Assign the id of the complex to the complex_id variable of every metamer
            Complex_id = node.complex_id

            Barycenter_XX = g.node(g.complex(i)).barycentre_XX
            Barycenter_YY = g.node(g.complex(i)).barycentre_YY
            Barycenter_ZZ = g.node(g.complex(i)).barycentre_ZZ
            if constant_supply:
                Photo_unit_surf = "none"
            else:
                Photo_unit_surf = node.Photo_unit_surf
            growth_fruit = node.growth_fruit
            RGR_fruit = node.RGR_fruit
            init_fruit_DW = node.init_fruit_dry_weight
            final_fruit_DW = node.fruit_dry_weight
            total_photo = node.tot_C_supply           #to compute total C_supply of metamer

            if node.leaf_area > 0:
                organ_type = "shoot"
            elif i in list_trunk:
                organ_type = "old_wood"
            else:
                organ_type = "trunk"
    #        pdb.set_trace()
            f.write(str(i)+'\t'+ str(Complex_id)+'\t'+ organ_type+'\t'+ str(XX)+'\t'+ str(YY)+"\t"+ str(ZZ)+"\t"+ str(Photo_unit_surf)+"\t"+
                    str(init_fruit_DW)+"\t"+ str(growth_fruit)+"\t"+ str(RGR_fruit)+"\t"+ str(total_photo)+'\t'+
                    str(Barycenter_XX)+'\t'+ str(Barycenter_YY)+'\t'+ str(Barycenter_ZZ)+'\t'+'\n')
        print("spatial output recorded" + nn)

    elif output_length == "complete":
        print("")
        print("Start writing complete output file, with properties of neighborhood")
        nn = pre + "\\neighborhood_spatial_" + suff
        f = open (nn,"w")
        f.write("Vid"+'\t'+ "Complex_id"+'\t'+"organ_type"+'\t'+ "X"+'\t'+ "Y"+'\t'+ "Z"+'\t'+ "Photo_unit_surf"+'\t'+ "init_fruit_DW"+'\t'+ "growth_fruit"+'\t'+ "RGR_fruit"+'\t'+ "total_photo"+'\t'+
                "Photo_neigh005" +'\t'+ "Fruit_neigh005" +'\t'+ "Photo_neigh015" +'\t'+ "Fruit_neigh015" +'\t'+
                "Photo_neigh025" +'\t'+ "Fruit_neigh025" +'\t'+ "Photo_neigh035" +'\t'+ "Fruit_neigh035" +'\t'+
                "Photo_neigh045" +'\t'+ "Fruit_neigh045" +'\t'+ "Photo_neigh055" +'\t'+ "Fruit_neigh055" +'\t'+
                "Photo_neigh065" +'\t'+ "Fruit_neigh065" +'\t'+
                "Photo_neigh075" +'\t'+ "Fruit_neigh075" +'\t'+
                "Photo_neigh095" +'\t'+ "Fruit_neigh095" +'\t'+
                "Photo_neigh115" +'\t'+ "Fruit_neigh115" +'\t'+
                "Photo_neigh135" +'\t'+ "Fruit_neigh135" +'\t'+
                "Axil_neigh005" +'\t'+ "Axil_neigh015" +'\t'+
                "Axil_neigh025" +'\t'+ "Axil_neigh035" +'\t'+
                "Axil_neigh045" +'\t'+ "Axil_neigh055" +'\n')

        first_metamer = next(g.component_roots_at_scale_iter(g.root, scale=selected_scale))
        list_trunk = g.Trunk(first_metamer, Scale = selected_scale)

        print("Computes properties in neighborhood")
        for i in g.vertices(scale = selected_scale):
            node = g.node(i)
            XX = node.XX
            YY = node.YY
            ZZ = node.ZZ
            node.complex_id = g.complex(i)     # Assign the id of the complex to the complex_id variable of every metamer
            Complex_id = node.complex_id
            if constant_supply:
                Photo_unit_surf = "none"
            else:
                Photo_unit_surf = node.Photo_unit_surf
            growth_fruit = node.growth_fruit
            RGR_fruit = node.RGR_fruit
            init_fruit_DW = node.init_fruit_dry_weight
            final_fruit_DW = node.fruit_dry_weight
            total_photo = node.tot_C_supply           #to compute total C_supply of metamer

            if node.leaf_area > 0:
                organ_type = "shoot"
            elif i in list_trunk:
                organ_type = "old_wood"
            else:
                organ_type = "trunk"
            #print "vertex is", i
            neighborhood(g, i)
            Axil_neigh005 = node.axil_neigh005
            Axil_neigh015 = node.axil_neigh015
            Axil_neigh025 = node.axil_neigh025
            Axil_neigh035 = node.axil_neigh035
            Axil_neigh045 = node.axil_neigh045
            Axil_neigh055 = node.axil_neigh055

            Photo_neigh005 = node.photo_neigh005
            Fruit_neigh005 = node.fruit_neigh005
            Photo_neigh015 = node.photo_neigh015
            Fruit_neigh015 = node.fruit_neigh015
            Photo_neigh025 = node.photo_neigh025
            Fruit_neigh025 = node.fruit_neigh025
            Photo_neigh035 = node.photo_neigh035
            Fruit_neigh035 = node.fruit_neigh035
            Photo_neigh045 = node.photo_neigh045
            Fruit_neigh045 = node.fruit_neigh045
            Photo_neigh055 = node.photo_neigh055
            Fruit_neigh055 = node.fruit_neigh055
            Photo_neigh065 = node.photo_neigh065
            Fruit_neigh065 = node.fruit_neigh065
            Photo_neigh075 = node.photo_neigh075
            Fruit_neigh075 = node.fruit_neigh075
            #Photo_neigh085 = node.photo_neigh085
            #Fruit_neigh085 = node.fruit_neigh085
            Photo_neigh095 = node.photo_neigh095
            Fruit_neigh095 = node.fruit_neigh095
            #Photo_neigh105 = node.photo_neigh105
            #Fruit_neigh105 = node.fruit_neigh105
            Photo_neigh115 = node.photo_neigh115
            Fruit_neigh115 = node.fruit_neigh115
            #Photo_neigh125 = node.photo_neigh125
            #Fruit_neigh125 = node.fruit_neigh125
            Photo_neigh135 = node.photo_neigh135
            Fruit_neigh135 = node.fruit_neigh135
            #Photo_neigh145 = node.photo_neigh145
            #Fruit_neigh145 = node.fruit_neigh145
            #Photo_neigh155 = node.photo_neigh155
            #Fruit_neigh155 = node.fruit_neigh155
            #Photo_neigh165 = node.photo_neigh165
            #Fruit_neigh165 = node.fruit_neigh165
            #Photo_neigh175 = node.photo_neigh175
            #Fruit_neigh175 = node.fruit_neigh175
    #        pdb.set_trace()
            f.write(str(i)+'\t'+ str(Complex_id)+'\t'+ organ_type+'\t'+ str(XX)+'\t'+ str(YY)+"\t"+ str(ZZ)+"\t"+ str(Photo_unit_surf)+"\t"+ str(init_fruit_DW)+"\t"+ str(growth_fruit)+"\t"+ str(RGR_fruit)+"\t"+ str(total_photo)+'\t'+
                    str(Photo_neigh005)+'\t'+ str(Fruit_neigh005) +'\t'+ str(Photo_neigh015)+'\t'+ str(Fruit_neigh015) +'\t'+
                    str(Photo_neigh025)+'\t'+ str(Fruit_neigh025) +'\t'+ str(Photo_neigh035)+'\t'+ str(Fruit_neigh035) +'\t'+
                    str(Photo_neigh045)+'\t'+ str(Fruit_neigh045) +'\t'+ str(Photo_neigh055)+'\t'+ str(Fruit_neigh055) +'\t'+
                    str(Photo_neigh065)+'\t'+ str(Fruit_neigh065) +'\t'+
                    str(Photo_neigh075)+'\t'+ str(Fruit_neigh075) +'\t'+
                    str(Photo_neigh095)+'\t'+ str(Fruit_neigh095) +'\t'+
                    str(Photo_neigh115)+'\t'+ str(Fruit_neigh115) +'\t'+
                    str(Photo_neigh135)+'\t'+ str(Fruit_neigh135) +'\t'+
                    str(Axil_neigh005)+'\t'+ str(Axil_neigh015) +'\t'+
                    str(Axil_neigh025)+'\t'+ str(Axil_neigh035) +'\t'+
                    str(Axil_neigh045)+'\t'+ str(Axil_neigh055) +'\n')
                    #str(Photo_neigh005)+'\t'+ str(Fruit_neigh005) +'\t'+ str(Photo_neigh015)+'\t'+ str(Fruit_neigh015) +'\t'+
                    #str(Photo_neigh025)+'\t'+ str(Fruit_neigh025) +'\t'+ str(Photo_neigh035)+'\t'+ str(Fruit_neigh035) +'\t'+
                    #str(Photo_neigh045)+'\t'+ str(Fruit_neigh045) +'\t'+ str(Photo_neigh055)+'\t'+ str(Fruit_neigh055) +'\t'+
                    #str(Photo_neigh065)+'\t'+ str(Fruit_neigh065) +'\t'+ str(Photo_neigh075)+'\t'+ str(Fruit_neigh075) +'\t'+
                    #str(Photo_neigh085)+'\t'+ str(Fruit_neigh085) +'\t'+ str(Photo_neigh095)+'\t'+ str(Fruit_neigh095) +'\t'+
                    #str(Photo_neigh105)+'\t'+ str(Fruit_neigh105) +'\t'+ str(Photo_neigh115)+'\t'+ str(Fruit_neigh115) +'\t'+
                    #str(Photo_neigh125)+'\t'+ str(Fruit_neigh125) +'\t'+ str(Photo_neigh135)+'\t'+ str(Fruit_neigh135) +'\t'+
                    #str(Photo_neigh145)+'\t'+ str(Fruit_neigh145) +'\t'+ str(Photo_neigh155)+'\t'+ str(Fruit_neigh155) +'\t'+
                    #str(Photo_neigh165)+'\t'+ str(Fruit_neigh165) +'\t'+ str(Photo_neigh175)+'\t'+ str(Fruit_neigh175) +"\n")
        print("spatial output recorded" + nn)

def define_scale_specific_parameters(partition):
    """
    Define selected_scale, replace_scale and scale_name parameters based on the chosen "partition".

    :Parameters:
        - 'MTGfile' (MTG) - name of the multi scale tree graph of the plant
        - 'day' - first simulated timestep (first day of simulation if the unit for timestep is days)
        - 'N' - number of simulated timesteps
        - 'h' (real) - friction parameter

    :Returns:
        - 'nn' - path and filename of the output file

    """
    if partition == metamer:
        selected_scale = 3
        replace_scale = False
        scale_name = "M"

    elif partition == growth_unit:
        selected_scale = 2
        replace_scale = False
        scale_name = "GU"

    elif partition == leafy_shoot_axis_and_root:
        selected_scale = 2
        replace_scale = True
        scale_name = "SA"

    elif partition == branches_and_root:
        selected_scale = 2
        replace_scale = True
        scale_name = "BR1"

    elif partition == branches_interbranches_and_root:
        selected_scale = 2
        replace_scale = True
        scale_name = "BR2"

    elif partition == trunk_branches_shoots:
        selected_scale = 2
        replace_scale = True
        scale_name = "TBS"

    elif partition == fruiting_unit_and_root:
        selected_scale = 2
        replace_scale = True
        scale_name = "FU"
    return(selected_scale, replace_scale, scale_name)

def iterations(Dists, time_Dists, scale_name_prec, MTGfile_prec, g, MTGfile, day, N, partition, h=0.5, constant_supply = False, init_root = False, output_length="complete"):
    """
    Iterates the computation of C_demands, C_supplies, C_allocation and growth
        over N iterations, starting from "day", using distance parameter "h"

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'MTGfile' (MTG) - name of the multi scale tree graph of the plant
        - 'day' - the current timestep
        - 'N' - number of simulated timesteps
        - 'partition' - the function containing the rules used to set topological boundaries of the new scale
        - 'selected_scale' (int) - topological scale of the plant components between which carbon allocation is computed
        - 'h' (real) - friction parameter
        - 'constant_supply' (boolean) - If True, photosynthesis is calculated with a radiative model,
            otherwise it is set to a constant per unit surface value
        - 'init_root' (boolean) - If 'True', a root is initialized by calling the proxy_root function
        - 'replace_scale' - If 'True', scale '2' is substituted by calling the 'replace_scale_fun' function

    :Returns:
        - 'g' (MTG) - a multi scale tree graph with initial and updated values of organs dry weight after the simulation period.
        - 'nn' - an output file that can be analyzed, e.g. in the R environment

    """
#    pdb.set_trace()
    beginning_time = time.time()
    for i in range (day, day+N):
        print("")
#        pdb.set_trace()
        if constant_supply == False:
            radiative_model(g, DayToExtract = i) #call RATP
        GDD, DD = meteo(i)                # obtain Cumulated and current day Degree Days from meteorological file
        print("Iteration: on DOY ", i, "accumulated Degree Days from bloom are ", GDD, "; Degree Days of the current day are ", DD)

        if i == day:
#            pdb.set_trace()
            selected_scale, replace_scale, scale_name = define_scale_specific_parameters(partition)
            main_init(g, GDD, selected_scale, init_root, replace_scale, partition)

            if not Dists or MTGfile != MTGfile_prec or scale_name != scale_name_prec:
                scale_name_prec = scale_name
                MTGfile_prec = MTGfile
                beg_time_dist = time.time()

                Dists = distances(g, selected_scale)

                time_Dists = time.time() - beg_time_dist
                beginning_time = beginning_time + time_Dists
                print(str(time_Dists) + " to compute distance matrix")
            init_GDD = GDD - DD
        activities_bis_functional(g, GDD, DD, i, day, constant_supply)  # define sink strengths of the current day  #usually (g, 1000, 0.8, 0.0101)

        if selected_scale != max(g.scales()): # aggregate mtg values when working at coarse scales
            aggregate_multi(g, coarse_scale=selected_scale)

        C_allocation_SIMWAL2ses(g, h, selected_scale, Dists) # Carbon allocation

        if selected_scale != max(g.scales()):
            clean_properties(g, selected_scale)
            disaggregate(g, selected_scale, equality=False)
        daily_diagnostics(g)
        growth(g, DD)                                            # Add C_allocated to biomasses
        if i == (day + N -1):
            pre, suff = output_filename(MTGfile, day, N, h, scale_name, init_root)
            Final_outputs(g, init_GDD, GDD, pre, suff, N, i, day, constant_supply, output_length, beginning_time, time_Dists, selected_scale)
    return(Dists, time_Dists, scale_name_prec, MTGfile_prec)


def neighborhood(g, i):
    print("Compute neighborhood for ", i)
    # Neighborhood properties
    # Whole tree approach
    neighborhood005 = 0.05
    neighborhood015 = 0.15
    neighborhood025 = 0.25
    neighborhood035 = 0.35
    neighborhood045 = 0.45
    neighborhood055 = 0.55
    neighborhood065 = 0.65
    neighborhood075 = 0.75
#    neighborhood085 = 0.85
    neighborhood095 = 0.95
#    neighborhood105 = 0.105
    neighborhood115 = 1.15
#    neighborhood125 = 0.125
    neighborhood135 = 1.35
#    neighborhood145 = 0.145
#    neighborhood155 = 1.55
#    neighborhood165 = 0.165
#    neighborhood175 = 0.175
    node = g.node(i)
#    Photo_neigh005, Photo_neigh015, Photo_neigh025, Photo_neigh035, Photo_neigh045, Photo_neigh055, Photo_neigh065 = 0,0,0,0,0,0,0
#    Photo_neigh075, Photo_neigh085, Photo_neigh095, Photo_neigh105, Photo_neigh115, Photo_neigh125, Photo_neigh135 = 0,0,0,0,0,0,0
#    Photo_neigh145, Photo_neigh155, Photo_neigh165, Photo_neigh175 = 0,0,0,0,0,0,0

#    Fruit_neigh005, Fruit_neigh015, Fruit_neigh025, Fruit_neigh035, Fruit_neigh045, Fruit_neigh055, Fruit_neigh065 = 0,0,0,0,0,0,0
#    Fruit_neigh075, Fruit_neigh085, Fruit_neigh095, Fruit_neigh105, Fruit_neigh115, Fruit_neigh125, Fruit_neigh135 = 0,0,0,0,0,0,0
#    Fruit_neigh145, Fruit_neigh155, Fruit_neigh165, Fruit_neigh175 = 0,0,0,0,0,0,0

    Photo_neigh005, Photo_neigh015, Photo_neigh025, Photo_neigh035, Photo_neigh045, Photo_neigh055 = 0,0,0,0,0,0
    Photo_neigh065, Photo_neigh075, Photo_neigh095, Photo_neigh115, Photo_neigh135= 0,0,0,0,0

    Fruit_neigh005, Fruit_neigh015, Fruit_neigh025, Fruit_neigh035, Fruit_neigh045, Fruit_neigh055 = 0,0,0,0,0,0
    Fruit_neigh065, Fruit_neigh075, Fruit_neigh095, Fruit_neigh115, Fruit_neigh135 = 0,0,0,0,0

    Axil_neigh005, Axil_neigh015, Axil_neigh025, Axil_neigh035, Axil_neigh045, Axil_neigh055 = 0,0,0,0,0,0

    if node.fruit:
        for j in g.vertices(scale = g.max_scale()):
            dist = compute_distance(g,i,j)

            if dist <= neighborhood005:
                Photo_neigh005 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh005 += 1
                if g.node(j).leaf_area:
                    Axil_neigh005 += 1

            if dist <= neighborhood015:
                Photo_neigh015 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh015 += 1
                if g.node(j).leaf_area:
                    Axil_neigh015 += 1

            if dist <= neighborhood025:
                Photo_neigh025 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh025 += 1
                if g.node(j).leaf_area:
                    Axil_neigh025 += 1

            if dist <= neighborhood035:
                Photo_neigh035 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh035 += 1
                if g.node(j).leaf_area:
                    Axil_neigh035 += 1

            if dist <= neighborhood045:
                Photo_neigh045 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh045 += 1
                if g.node(j).leaf_area:
                    Axil_neigh045 += 1

            if dist <= neighborhood055:
                Photo_neigh055 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh055 += 1
                if g.node(j).leaf_area:
                    Axil_neigh055 += 1

            if dist <= neighborhood065:
                Photo_neigh065 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh065 += 1

            if dist <= neighborhood075:
                Photo_neigh075 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh075 += 1

            #if dist <= neighborhood085:
            #    Photo_neigh085 += g.node(j).tot_C_supply
            #    if g.node(j).fruit:
            #        Fruit_neigh085 += 1

            if dist <= neighborhood095:
                Photo_neigh095 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh095 += 1

            #if dist <= neighborhood105:
            #    Photo_neigh105 += g.node(j).tot_C_supply
            #    if g.node(j).fruit:
            #        Fruit_neigh105 += 1

            if dist <= neighborhood115:
                Photo_neigh115 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh115 += 1

            #if dist <= neighborhood125:
            #    Photo_neigh125 += g.node(j).tot_C_supply
            #    if g.node(j).fruit:
            #        Fruit_neigh125 += 1

            if dist <= neighborhood135:
                Photo_neigh135 += g.node(j).tot_C_supply
                if g.node(j).fruit:
                    Fruit_neigh135 += 1

            #if dist <= neighborhood145:
            #    Photo_neigh145 += g.node(j).tot_C_supply
            #    if g.node(j).fruit:
            #        Fruit_neigh145 += 1

            #if dist <= neighborhood155:
            #    Photo_neigh155 += g.node(j).tot_C_supply
            #    if g.node(j).fruit:
            #        Fruit_neigh155 += 1

            #if dist <= neighborhood165:
            #    Photo_neigh165 += g.node(j).tot_C_supply
            #    if g.node(j).fruit:
            #        Fruit_neigh165 += 1
            #
            #if dist <= neighborhood175:
            #    Photo_neigh175 += g.node(j).tot_C_supply
            #    if g.node(j).fruit:
            #        Fruit_neigh175 += 1
        #print "vertex " + str(i)
        #print "in 5cm " + str(Fruit_neigh005)+ " fruits and " + str(Photo_neigh005) + " accumulated photosynthesis"
        #print "in 25cm " + str(Fruit_neigh025) + " fruits and " + str(Photo_neigh025) + " accumulated photosynthesis"
        #print "in 75cm " + str(Fruit_neigh075) + " fruits and " + str(Photo_neigh075) + " accumulated photosynthesis"
        #print "in 105cm " + str(Fruit_neigh105) + " fruits and " + str(Photo_neigh105) + " accumulated photosynthesis"
        #print "in 175cm " + str(Fruit_neigh175) + " fruits and " + str(Photo_neigh175) + " accumulated photosynthesis"
    node.photo_neigh005 = Photo_neigh005
    node.fruit_neigh005 = Fruit_neigh005
    node.axil_neigh005 = Axil_neigh005

    node.photo_neigh015 = Photo_neigh015
    node.fruit_neigh015 = Fruit_neigh015
    node.axil_neigh015 = Axil_neigh015

    node.photo_neigh025 = Photo_neigh025
    node.fruit_neigh025 = Fruit_neigh025
    node.axil_neigh025 = Axil_neigh025

    node.photo_neigh035 = Photo_neigh035
    node.fruit_neigh035 = Fruit_neigh035
    node.axil_neigh035 = Axil_neigh035

    node.photo_neigh045 = Photo_neigh045
    node.fruit_neigh045 = Fruit_neigh045
    node.axil_neigh045 = Axil_neigh045

    node.photo_neigh055 = Photo_neigh055
    node.fruit_neigh055 = Fruit_neigh055
    node.axil_neigh055 = Axil_neigh055

    node.photo_neigh065 = Photo_neigh065
    node.fruit_neigh065 = Fruit_neigh065
    node.photo_neigh075 = Photo_neigh075
    node.fruit_neigh075 = Fruit_neigh075
    #node.photo_neigh085 = Photo_neigh085
    #node.fruit_neigh085 = Fruit_neigh085
    node.photo_neigh095 = Photo_neigh095
    node.fruit_neigh095 = Fruit_neigh095
    #node.photo_neigh105 = Photo_neigh105
    #node.fruit_neigh105 = Fruit_neigh105
    node.photo_neigh115 = Photo_neigh115
    node.fruit_neigh115 = Fruit_neigh115
    #node.photo_neigh125 = Photo_neigh125
    #node.fruit_neigh125 = Fruit_neigh125
    node.photo_neigh135 = Photo_neigh135
    node.fruit_neigh135 = Fruit_neigh135
    #node.photo_neigh145 = Photo_neigh145
    #node.fruit_neigh145 = Fruit_neigh145
    #node.photo_neigh155 = Photo_neigh155
    #node.fruit_neigh155 = Fruit_neigh155
    #node.photo_neigh165 = Photo_neigh165
    #node.fruit_neigh165 = Fruit_neigh165
    #node.photo_neigh175 = Photo_neigh175
    #node.fruit_neigh175 = Fruit_neigh175


def compartmental_budget(g,pre, suff):
    """
    Computes growth at compartmental scale for the whole plant

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant

    :Returns:
        - Growth_fruit - the total growth achieved by the fruit compartment
        - Growth_wood - the total growth achieved by the wood compartment

    """
    print("Computes compartmental budget for the whole simulation: ")

    nn = pre + "\\compBudget_" + suff
    budget = open (nn,"w")
    budget.write("Growth_old_wood"+'\t'+ "Growth_shoot"+'\t'+"Growth_fruit"+'\t'+ "Growth_leaves"+'\t'+"Growth_root"+'\t'+
                 "Nb_fruits"+'\t'+ "Nb_shoots" +'\t'+ "Total_shoot_length" +'\t'+ "Total_shoot_volume" +'\t'+
                 "Total_old_wood_biomass"+'\t'+ "Total_shoot_biomass" +'\t'+ "Total_fruit_biomass" +'\t'+ "Total_leafy_biomass" +'\t'+ "Total_root_biomass"+'\n')

    Growth_wood = 0
    Growth_old_wood = 0
    Growth_shoot = 0
    Growth_fruit = 0
    Growth_leaves = 0
    Growth_root = 0
    Nb_fruits = 0
    Nb_shoots = 0
    Total_woody_biomass = 0
    Total_root_biomass = 0
    Total_fruit_biomass = 0
    Total_leafy_biomass = 0
    Total_shoot_length = 0
    Total_shoot_volume = 0
    Total_shoot_biomass = 0

    extremities = g.Extremities(g.vertices(scale=g.max_scale())[0])
    print("extremities" + str(extremities))

    for i in g.vertices(scale=3):
        node = g.node(i)

        if i == g.root:
            Total_root_biomass  += node.woody_mass
            Growth_root += node.growth_wood
            print("I am in the ROOT!")

        if (node.fruit > 0):
            Nb_fruits += 1
            Growth_fruit  += node.growth_fruit
            Total_fruit_biomass += node.init_fruit_dry_weight

        if node.leaf_area > 0:
            Growth_shoot += node.growth_wood
            Total_shoot_biomass += node.woody_mass
            Growth_leaves += node.growth_leaf
            Total_leafy_biomass += node.leaf_dry_weight

            Total_shoot_length += node.length
            Total_shoot_volume += node.volume

            if i in extremities:
                Nb_shoots += 1

        Total_woody_biomass += node.woody_mass
        Growth_wood += node.growth_wood

    Growth_old_wood = Growth_wood - (Growth_shoot + Growth_root)
    Total_old_wood_biomass = Total_woody_biomass - (Total_shoot_biomass + Total_root_biomass)

    budget.write(str(Growth_old_wood) +'\t'+ str(Growth_shoot) +'\t'+ str(Growth_fruit) +'\t'+ str(Growth_leaves) + '\t'+ str(Growth_root) +'\t'+
                 str(Nb_fruits) +'\t'+  str(Nb_shoots) +'\t'+ str(Total_shoot_length) +'\t'+ str(Total_shoot_volume) +'\t'+
                 str(Total_old_wood_biomass) +'\t'+  str(Total_shoot_biomass) +'\t'+ str(Total_fruit_biomass) +'\t'+  str(Total_leafy_biomass) +'\t'+ str(Total_root_biomass) +'\n')
    budget.close()

    print("(total carbon allocated to fruit is " + str(Growth_fruit), ")")
    print("(total carbon allocated to old wood  is " + str(Growth_old_wood), ")")
    print("(total carbon allocated to shoots  is " + str(Growth_shoot), ")")
    print("(total carbon allocated to root  is " + str(Growth_root), ")")
    print("compartmental budget computed")
    return(Growth_fruit, Growth_old_wood)

def daily_diagnostics(g):
    """
    Provides aggregated information about current day C-allocation computed at max_scale

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant

    :Returns:
    """
    print("Computes diagnostics for the current day: ")

    N_tot_components = len( g.vertices(scale = g.max_scale()) )
    total_C_demand = 0
    total_C_supply = 0
    total_C_allocation = 0
    N_excess = 0
    for i in g.vertices(scale = g.max_scale()):
        node = g.node(i)
        total_C_demand += node.C_demand
        total_C_supply += node.C_supply
        total_C_allocation += node.C_allocation

        if node.C_excess > 0:
            N_excess += 1

    # print "  - total carbon supply " + str(total_C_supply)
    # print "  - total carbon demand " + str(total_C_demand)
    # print "  - carbon allocated " + str(total_C_allocation)
    # print "  - carbon demand satisfied for " + str(N_excess) + " of " + str(N_tot_components) + " components"

    tot_C_Supply = "tot_C_supply = " + str(total_C_supply)
    tot_C_Demand = "tot_C_demand = " + str(total_C_demand)
    tot_C_Allocated = "tot_C_allocated = " + str(total_C_allocation)
    nb_C_D_satisfied = "nb_C_demand_satisfied" + str(N_excess)
    nb_C_Demand = "nb_C_demand" + str(N_tot_components)

    return(tot_C_Supply,tot_C_Demand,tot_C_Allocated,nb_C_D_satisfied,nb_C_Demand)


def tree_axil(g):
    """
    Compute axialization (ratio between leaf dry weight and shoot dry weight) for each leafy shoot of an mtg.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant

    :Returns:
        - 'g' (MTG) - mtg with the "axil" property containing the axialization of the current leafy shoot

    """
    print("Compute axilization in this tree")
    vtxs = list()
    for i in g.vertices( scale=g.max_scale() ):
        if g.node(i).leaf_area:
            vtxs.append(i)

    while vtxs:
        vtxs = list(vtxs)
        i = vtxs[0]
        Axil, vtxs_met_i = shoot_axil(g,i)

        for i in vtxs_met_i:
            g.node(i).axil = Axil
        # vtxs = set(vtxs) ^ set(vtxs_met_i)

        vtxs = set(vtxs) ^ (set(vtxs_met_i).intersection(vtxs))
        # print (len(vtxs))
    print("axilization computed")

def shoot_axil(g,i):
    """
    Computes axillization of the shoot related to a given vertex at max scale.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'i' - vertex
    :Returns:
        - Axil - the axilization value of the shoot to which the vertex belongs
        - vtxs_met_i - the vertices composing the shoot
    """
    # print "Computing axilization of vertex ", i
    shoot_down = i
    shoot_up = i
    vtxs_met_i = list()
    vtxs_met_i.append(i)

    Leaf_dry_weight = 0
    Shoot_biomass = 0
    Leaf_dry_weight += g.node(i).leaf_dry_weight
    Shoot_biomass += g.node(i).woody_mass
    while g.node( g.Successor(shoot_up) ):
        #print shoot_up
        shoot_up = g.Successor(shoot_up)
        Leaf_dry_weight += g.node(shoot_up).leaf_dry_weight
        Shoot_biomass += g.node(shoot_up).woody_mass
        vtxs_met_i.append(shoot_up)

    test_leaf = True
    while g.node(g.Predecessor(shoot_down)):
        #print shoot_down
        if g.node(g.Predecessor(shoot_down)).leaf_dry_weight:
            shoot_down = g.Predecessor(shoot_down)
            Leaf_dry_weight += g.node(shoot_down).leaf_dry_weight
            Shoot_biomass += g.node(shoot_down).woody_mass
            vtxs_met_i.append(shoot_down)
        else:
            shoot_down = g.Predecessor(shoot_down)
            test_leaf= False

    Axil = old_div(Leaf_dry_weight, Shoot_biomass)
    # print Axil, vtxs_met_i
    return(Axil, vtxs_met_i)

def Z_adjustment(g):
    """
    Translates ZZ coordinates of an mtg in order to have the minimum ZZ value equal to zero.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
    :Returns:
        - minimum_Z - the minimum ZZ coordinate value before adjustment
    """
    firts = True
    for i in g.vertices(scale=3):
        if firts:
            minimum_Z = g.node(i).ZZ
            firts = False
        if g.node(i).ZZ < minimum_Z:
            minimum_Z = g.node(i).ZZ
    if minimum_Z != 0:
        print("perform Z-adjustment")
        for i in g.vertices(scale=3):
            g.node(i).ZZ = g.node(i).ZZ - minimum_Z
    return (minimum_Z)


def batch_sim(h_values = [0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,8,16],  # 0.25,0.5,1,2,4,8,16
              mtg_names = ['Med_cul05_1996_6.mtg'], output_length ="synthetic"):
    """
    Launch a batch of simulations with all combinations of tree structures and friction parameters in input,
    and for all scales defined in the code.

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
    :Remarks:
        - distance matrices are not re-computed if only the scale or the friction parameter are changed.
    """
    #'high_1996_7.mtg','High_cul05_1996_6.mtg','High_cul10_1996_6.mtg','High_cul45_1996_6.mtg','Med_cul05_1996_6.mtg','Med_cul10_1996_6.mtg','Med_cul45_1996_6.mtg'
    # 'normal_fuji_1997_7.mtg', 'Med_cul05_1996_6.mtg','Med_cul10_1996_6.mtg'
    Dists, time_Dists, scale_name_prec, MTGfile_prec = None, None, None, None
    for mtg_name in mtg_names: # metamer,growth_unit, trunk_branches_shoots, branches_and_root, fruiting_unit_and_root
        for topo_scale in [growth_unit, metamer, trunk_branches_shoots, branches_and_root, fruiting_unit_and_root]: # metamer growth_unit, leafy_shoot_axis_and_root, branches_interbranches_and_root,fruiting_unit_and_root, branches_and_root, trunk_branches_shoots
            for j in h_values:
                print("")
                print("Input filename:",mtg_name, topo_scale, "h value", j)
                print("So we start: ciao!")
                # if __name__ == '__main__':
                MTGfile = mtg_name
                g = MTG(mtg_name)
                Z_adjustment(g)
                Dists, time_Dists, scale_name_prec, MTGfile_prec = iterations(Dists, time_Dists, scale_name_prec, MTGfile_prec,
                                                                                  g, MTGfile, day = 182, N=1, partition=topo_scale, h=j,
                                                                                  constant_supply = False, init_root = True, output_length=output_length)



# Old functions

def output_2(g, nn, constant_supply, selected_scale = 3):
    """
    Generates an output file from the simulation containing neighborhood properties, including:
        'vid' - identifier of plant part at fine scale,
        'Complex_id' - identifier of complex of 'vid',
        'organ_type' - type of organ of 'vid',
        'X', 'Y', 'Z' - basal spatial coordinates of 'vid'
        'Photo_unit_surf' - assimilated carbon per unit surface
        'init_fruit_DW' - initial fruit weight on 'vid'
        'growth_fruit' - absolute growth in fruit dry weight present on 'vid'
        'RGR_fruit' - relative growth rate of fruit present on 'vid'

    :Parameters:
        - 'g' (MTG) - multi scale tree graph of the plant
        - 'nn' - path and filename generated by the 'output_filename' function and passed to the 'output_2' function for generating an output
        - 'constant_supply' (boolean) - passed to the 'output_2' function for generating an output
        - 'selected_scale' - scale at which producing the output

    :Returns:
         - a file of output 'nn'.

    """
    print("")
    print("Start writing file with properties of neighborhood")
    f = open (nn,"w")
    f.write("Vid"+'\t'+ "Complex_id"+'\t'+"organ_type"+'\t'+ "X"+'\t'+ "Y"+'\t'+ "Z"+'\t'+ "Photo_unit_surf"+'\t'+ "init_fruit_DW"+'\t'+ "growth_fruit"+'\t'+ "RGR_fruit"+'\t'+ "total_photo"+'\t'+
            "Photo_neigh005" +'\t'+ "Fruit_neigh005" +'\t'+ "Photo_neigh015" +'\t'+ "Fruit_neigh015" +'\t'+
            "Photo_neigh025" +'\t'+ "Fruit_neigh025" +'\t'+ "Photo_neigh035" +'\t'+ "Fruit_neigh035" +'\t'+
            "Photo_neigh045" +'\t'+ "Fruit_neigh045" +'\t'+ "Photo_neigh055" +'\t'+ "Fruit_neigh055" +'\t'+
            "Photo_neigh065" +'\t'+ "Fruit_neigh065" +'\t'+
            "Photo_neigh075" +'\t'+ "Fruit_neigh075" +'\t'+
            "Photo_neigh095" +'\t'+ "Fruit_neigh095" +'\t'+
            "Photo_neigh115" +'\t'+ "Fruit_neigh115" +'\t'+
            "Photo_neigh135" +'\t'+ "Fruit_neigh135" +'\n')

            #"Photo_neigh005" +'\t'+ "Fruit_neigh005" +'\t'+ "Photo_neigh015" +'\t'+ "Fruit_neigh015" +'\t'+
            #"Photo_neigh025" +'\t'+ "Fruit_neigh025" +'\t'+ "Photo_neigh035" +'\t'+ "Fruit_neigh035" +'\t'+
            #"Photo_neigh045" +'\t'+ "Fruit_neigh045" +'\t'+ "Photo_neigh055" +'\t'+ "Fruit_neigh055" +'\t'+
            #"Photo_neigh065" +'\t'+ "Fruit_neigh065" +'\t'+ "Photo_neigh075" +'\t'+ "Fruit_neigh075" +'\t'+
            #"Photo_neigh085" +'\t'+ "Fruit_neigh085" +'\t'+ "Photo_neigh095" +'\t'+ "Fruit_neigh095" +'\t'+
            #"Photo_neigh105" +'\t'+ "Fruit_neigh105" +'\t'+ "Photo_neigh115" +'\t'+ "Fruit_neigh115" +'\t'+
            #"Photo_neigh125" +'\t'+ "Fruit_neigh125" +'\t'+ "Photo_neigh135" +'\t'+ "Fruit_neigh135" +'\t'+
            #"Photo_neigh145" +'\t'+ "Fruit_neigh145" +'\t'+ "Photo_neigh155" +'\t'+ "Fruit_neigh155" +'\t'+
            #"Photo_neigh165" +'\t'+ "Fruit_neigh165" +'\t'+ "Photo_neigh175" +'\t'+ "Fruit_neigh175" +'\t'+'\n')

    first_metamer = next(g.component_roots_at_scale_iter(g.root, scale=selected_scale))
    list_trunk = g.Trunk(first_metamer, Scale = selected_scale)

    print("Computes properties in neighborhood")
    for i in g.vertices(scale = selected_scale):
        node = g.node(i)
        XX = node.XX
        YY = node.YY
        ZZ = node.ZZ
        node.complex_id = g.complex(i)     # Assign the id of the complex to the complex_id variable of every metamer
        Complex_id = node.complex_id
        if constant_supply:
            Photo_unit_surf = "none"
        else:
            Photo_unit_surf = node.Photo_unit_surf
        growth_fruit = node.growth_fruit
        RGR_fruit = node.RGR_fruit
        init_fruit_DW = node.init_fruit_dry_weight
        final_fruit_DW = node.fruit_dry_weight
        total_photo = node.tot_C_supply           #to compute total C_supply of metamer

        if node.leaf_area > 0:
            organ_type = "shoot"
        elif i in list_trunk:
            organ_type = "old_wood"
        else:
            organ_type = "trunk"
        #print "vertex is", i
        neighborhood(g, i)
        Photo_neigh005 = node.photo_neigh005
        Fruit_neigh005 = node.fruit_neigh005
        Photo_neigh015 = node.photo_neigh015
        Fruit_neigh015 = node.fruit_neigh015
        Photo_neigh025 = node.photo_neigh025
        Fruit_neigh025 = node.fruit_neigh025
        Photo_neigh035 = node.photo_neigh035
        Fruit_neigh035 = node.fruit_neigh035
        Photo_neigh045 = node.photo_neigh045
        Fruit_neigh045 = node.fruit_neigh045
        Photo_neigh055 = node.photo_neigh055
        Fruit_neigh055 = node.fruit_neigh055
        Photo_neigh065 = node.photo_neigh065
        Fruit_neigh065 = node.fruit_neigh065
        Photo_neigh075 = node.photo_neigh075
        Fruit_neigh075 = node.fruit_neigh075
        #Photo_neigh085 = node.photo_neigh085
        #Fruit_neigh085 = node.fruit_neigh085
        Photo_neigh095 = node.photo_neigh095
        Fruit_neigh095 = node.fruit_neigh095
        #Photo_neigh105 = node.photo_neigh105
        #Fruit_neigh105 = node.fruit_neigh105
        Photo_neigh115 = node.photo_neigh115
        Fruit_neigh115 = node.fruit_neigh115
        #Photo_neigh125 = node.photo_neigh125
        #Fruit_neigh125 = node.fruit_neigh125
        Photo_neigh135 = node.photo_neigh135
        Fruit_neigh135 = node.fruit_neigh135
        #Photo_neigh145 = node.photo_neigh145
        #Fruit_neigh145 = node.fruit_neigh145
        #Photo_neigh155 = node.photo_neigh155
        #Fruit_neigh155 = node.fruit_neigh155
        #Photo_neigh165 = node.photo_neigh165
        #Fruit_neigh165 = node.fruit_neigh165
        #Photo_neigh175 = node.photo_neigh175
        #Fruit_neigh175 = node.fruit_neigh175
#        pdb.set_trace()
        f.write(str(i)+'\t'+ str(Complex_id)+'\t'+ organ_type+'\t'+ str(XX)+'\t'+ str(YY)+"\t"+ str(ZZ)+"\t"+ str(Photo_unit_surf)+"\t"+ str(init_fruit_DW)+"\t"+ str(growth_fruit)+"\t"+ str(RGR_fruit)+"\t"+ str(total_photo)+'\t'+
                str(Photo_neigh005)+'\t'+ str(Fruit_neigh005) +'\t'+ str(Photo_neigh015)+'\t'+ str(Fruit_neigh015) +'\t'+
                str(Photo_neigh025)+'\t'+ str(Fruit_neigh025) +'\t'+ str(Photo_neigh035)+'\t'+ str(Fruit_neigh035) +'\t'+
                str(Photo_neigh045)+'\t'+ str(Fruit_neigh045) +'\t'+ str(Photo_neigh055)+'\t'+ str(Fruit_neigh055) +'\t'+
                str(Photo_neigh065)+'\t'+ str(Fruit_neigh065) +'\t'+
                str(Photo_neigh075)+'\t'+ str(Fruit_neigh075) +'\t'+
                str(Photo_neigh095)+'\t'+ str(Fruit_neigh095) +'\t'+
                str(Photo_neigh115)+'\t'+ str(Fruit_neigh115) +'\t'+
                str(Photo_neigh135)+'\t'+ str(Fruit_neigh135) +'\n')
                #str(Photo_neigh005)+'\t'+ str(Fruit_neigh005) +'\t'+ str(Photo_neigh015)+'\t'+ str(Fruit_neigh015) +'\t'+
                #str(Photo_neigh025)+'\t'+ str(Fruit_neigh025) +'\t'+ str(Photo_neigh035)+'\t'+ str(Fruit_neigh035) +'\t'+
                #str(Photo_neigh045)+'\t'+ str(Fruit_neigh045) +'\t'+ str(Photo_neigh055)+'\t'+ str(Fruit_neigh055) +'\t'+
                #str(Photo_neigh065)+'\t'+ str(Fruit_neigh065) +'\t'+ str(Photo_neigh075)+'\t'+ str(Fruit_neigh075) +'\t'+
                #str(Photo_neigh085)+'\t'+ str(Fruit_neigh085) +'\t'+ str(Photo_neigh095)+'\t'+ str(Fruit_neigh095) +'\t'+
                #str(Photo_neigh105)+'\t'+ str(Fruit_neigh105) +'\t'+ str(Photo_neigh115)+'\t'+ str(Fruit_neigh115) +'\t'+
                #str(Photo_neigh125)+'\t'+ str(Fruit_neigh125) +'\t'+ str(Photo_neigh135)+'\t'+ str(Fruit_neigh135) +'\t'+
                #str(Photo_neigh145)+'\t'+ str(Fruit_neigh145) +'\t'+ str(Photo_neigh155)+'\t'+ str(Fruit_neigh155) +'\t'+
                #str(Photo_neigh165)+'\t'+ str(Fruit_neigh165) +'\t'+ str(Photo_neigh175)+'\t'+ str(Fruit_neigh175) +"\n")
    print("spatial output recorded" + nn)
