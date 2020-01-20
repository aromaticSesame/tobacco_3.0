import pymatgen as pm
from ase.geometry.cell import cellpar_to_cell
import os
from datetime import datetime
import numpy as np
import warnings

tol = 1E-2 #tolerance for distances
scale = 10 #scale lattice constants by this factor
cgd_filename = 'RCSRnets-2019-06-01.cgd' #http://rcsr.anu.edu.au/systre

vnames = [
	'V','Er','Ti','Ce','S',
	'H','He','Li','Be','B',
	'C','N','O','F','Ne',
	'Na','Mg','Al','Si','P',
	'Cl','Ar','K','Ca','Sc',
	'Cr','Mn','Fe','Co','Ni'] #names of vertices
edge_center_name = 'Lr' #placeholder edge name

if edge_center_name in vnames:
	raise ValueError('Edge center name must not be in vnames',edge_center_name)

#initialize lists
topologies_all = [] #all topologies
groups_all = [] #all spacegroups
cellpars_all = [] #all [a,b,c,alpha,beta,gamma]
vertices_all = [] #all [x,y,z] fractional positions of vertices
edges_center_all = [] #all [x,y,z] fractional positions of edge centers
edges_head_all = [] #all [x,y,z] fractional positions of edge heads
cn_all = [] #all vertex coordination numbers, coded as dictionaries

#Make sure .cgd file is present
if not os.path.exists(cgd_filename):
	raise ValueError('Missing RCSR .cgd data file', cgd_filename)

#Forbidden names on Windows
forbidden_names = ['con','prn','aux','nul']

#Read info from .cgd file
with open(cgd_filename,'r') as r:
	for line in r:
		line = line.strip()

		#Initialize values for new topology
		if 'crystal' in line.lower():
			three_dim = True
			vertices = []
			edges_center = []
			edges_head = []
			cn = {}
			vertices_count = 0

		#Get the topology name
		elif 'name' in line.lower():
			topology_val = line.lower().split('name')[-1].replace('*','_star').replace('-','').strip()
			if topology_val in forbidden_names:
				topology_val += '0'

		#Get the spacegroup
		elif 'group' in line.lower():

			#Do not alter capitalization of spacegroups
			group_val = line.split('GROUP')[-1].split('group')[-1].strip()

			#Use updated group name of Cmca
			if group_val == 'Cmca':
				group_val = 'Cmce'

		#Get the lattice constants (with scale*(a,b,c))
		elif 'cell' in line.lower():
			cell_val = line.lower().split('cell')[-1]
			cell_val = [float(i) for i in cell_val.split()]
			cell_val[0] = cell_val[0]*scale
			cell_val[1] = cell_val[1]*scale
			cell_val[2] = cell_val[2]*scale

		#Get the vertices (make sure it's 3D) and get CNs
		elif 'node' in line.lower() or 'atom' in line.lower():
			vert_val = line.lower().split('node')[-1].split('atom')[-1].strip()
			vert_val = [i for i in vert_val.split()]

			#Make sure there is a coordination number and [x,y,z]
			if len(vert_val) != 5:
				three_dim = False
				continue

			vertices.append([float(vert_val[2]),float(vert_val[3]),float(vert_val[4])])
			vertices_count += 1

			#Make sure there are enough names for the vertices
			if vertices_count > len(vnames):
				raise ValueError('More verties than vnames for '+topology)

			#Make a coordination number dictionary
			cn[vnames[vertices_count-1]] = int(vert_val[1])

		#Get edge centers
		elif 'edge_center' in line.lower():
			edge_center_val = line.lower().split('edge_center')[-1].strip()
			edge_center_val = [float(i) for i in edge_center_val.split()]
			edges_center.append(edge_center_val)

			#Make sure there are [x,y,z] coordinates
			if len(edge_center_val) != 3:
				three_dim = False
				continue

		#Get edge endpoints
		elif 'edge' in line.lower():
			edge_val = line.lower().split('edge')[-1].strip()
			edge_val = [float(i) for i in edge_val.split()]
			edges_head.append(edge_val[0:3])

		#Store results for topology
		elif line.lower() == 'end':

			#Skip 2D topologies
			if not three_dim:
				continue

			#Skip weirdly formatted cgd entries
			if len(cn) != len(vertices):
				warnings.warn('Error: skipping '+topology_val+' because it is not formatted properly in .cgd file',Warning)
				continue
			elif len(edges_head) != len(edges_center):
				warnings.warn('Error: skipping '+topology_val+' because it is not formatted properly in .cgd file',Warning)
				continue				

			topologies_all.append(topology_val)
			groups_all.append(group_val)
			cellpars_all.append(cell_val)
			vertices_all.append(vertices)
			edges_center_all.append(edges_center)
			edges_head_all.append(edges_head)
			cn_all.append(cn)

		#Ignore NC nets (assumed to be at bottom of .cgd file)
		elif 'nc nets' in line.lower():
			break

#Make folders to store topology CIFs
if not os.path.exists('templates_database'):
	os.mkdir('templates_database')
if not os.path.exists('templates_errors'):
	os.mkdir('templates_errors')

#Cycle through all topologies and make CIFs
for i in range(0,len(topologies_all)):

	#Flag for skipping CIF generation
	bad = False

	#Get all .cgd info for given topology, i
	topology = topologies_all[i]
	group = groups_all[i]
	cellpars = cellpars_all[i]
	vertices = vertices_all[i]
	edges_center = edges_center_all[i]
	edges_head = edges_head_all[i]
	cn_vec = cn_all[i]

	#Make list of vertex and edge center symbols
	sym_vertices = []
	for j in range(len(vertices)):
		sym_vertices.append(vnames[j])
	sym_collection = sym_vertices+[edge_center_name]*len(edges_center)

	#Get lattice vectors (using ASE function because it's easy)
	lattice_vectors = cellpar_to_cell(cellpars)

	#Get vertex and edge positions
	basis_collection = np.array(vertices+edges_center)

	#Make pymatgen structure
	pm_structure = pm.Structure.from_spacegroup(group,lattice_vectors,sym_collection,basis_collection)
	pm_structure.merge_sites(mode='delete')	

	#Calculate distance between edge centers and edge ends
	dummy_edges = pm.Structure(lattice_vectors,[edge_center_name]*len(edges_head),edges_head)
	dummy_centers = pm.Structure(lattice_vectors,[edge_center_name]*len(edges_center),edges_center)
	d_tests = []
	for j, dummy_center in enumerate(dummy_centers):
		d_test = 2*dummy_center.distance(dummy_edges[j])
		d_tests.append(d_test)

	#Defined bond distance as 2*distance between edge centers and edge ends (i.e. vertices)
	bond_dists = [] 
	if np.abs(np.max(d_tests)-np.min(d_tests)) < 2*tol:
		bond_dists.append(np.average(d_tests))
	else:
		bond_dists.extend(d_tests)

	#Make lattice constants > bond dist
	if np.max(bond_dists) < scale:
		extend = tol+scale
	else:
		extend = tol+np.max(bond_dists)
	n_supercells = [np.ceil(extend/cellpars[0]),np.ceil(extend/cellpars[1]),np.ceil(extend/cellpars[2])]
	if n_supercells != [1,1,1]:
		pm_structure.make_supercell(n_supercells)

	#Get atoms of edge centers and vertices
	vertices_indices = [atom_idx for atom_idx, atom in enumerate(pm_structure) if atom.species_string != edge_center_name]
	edge_center_indices = [atom_idx for atom_idx, atom in enumerate(pm_structure) if atom.species_string == edge_center_name]

	#Make text for top of CIF
	top_text = 'data_'+topology+'\n'+'_audit_creation_date              '+datetime.today().strftime('%Y-%m-%d')+'\n'+"_audit_creation_method            'Pymatgen'\n"+"_symmetry_space_group_name_H-M    'P1'\n"+'_symmetry_Int_Tables_number       1\n'
	cellpar_text = 'loop_\n_symmetry_equiv_pos_as_xyz\n  x,y,z\n'+'_cell_length_a                    '+str(np.round(pm_structure.lattice.abc[0],4))+'\n'+'_cell_length_b                    '+str(np.round(pm_structure.lattice.abc[1],4))+'\n'+'_cell_length_c                    '+str(np.round(pm_structure.lattice.abc[2],4))+'\n'+'_cell_angle_alpha                 '+str(np.round(pm_structure.lattice.angles[0],4))+'\n'+'_cell_angle_beta                 '+str(np.round(pm_structure.lattice.angles[1],4))+'\n'+'_cell_angle_gamma                 '+str(np.round(pm_structure.lattice.angles[2],4))+'\n'
	pos_text = 'loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n'

	#Make (minimum image) distance matrix
	dist_mat = pm_structure.distance_matrix

	#Initialization
	bonded_pairs = [] #list for the indices of bonded vertices
	bonded_edge_centers = [] #list for the indices of bonded edge centers
	img_list = [] #list of image displacements
	d_list = [] #list of bond distances

	#Cycle through every vertex to find its bonded atoms
	for j, vertex_idx in enumerate(vertices_indices):

		#Initialization
		vertex_atom = pm_structure[vertex_idx] #Site object
		cn = cn_vec[vertex_atom.species_string] #int
		pm_structure[vertex_idx].index = j #store the index, excluding edge centers

		#Make string containing fract position for atom j
		pos_text += vertex_atom.species_string+str(j+1)+'     '+vertex_atom.species_string+'     '+str(np.round(vertex_atom.frac_coords[0],4))+'   '+str(np.round(vertex_atom.frac_coords[1],4))+'   '+str(np.round(vertex_atom.frac_coords[2],4))+'\n'
		
		#Find all edge centers connected to vertex j
		edge_overlap_indices = []
		for bond_dist in bond_dists:
			edges_shell_temp = pm_structure.get_neighbors_in_shell(pm_structure[vertex_idx].coords,bond_dist/2,tol,include_index=True)
			edges_shell = [k for k in edges_shell_temp if k[0].species_string == edge_center_name]
			for edge_shell in edges_shell:
				if edge_shell[2] not in edge_overlap_indices:
					edge_overlap_indices.append(edge_shell[2])

		#Check if an edge center must be counted again due to PBCs
		edge_overlap_to_add = []
		for edge_overlap_idx in edge_overlap_indices:
			edge_overlap_atom = pm_structure[edge_overlap_idx]

			d_vec = vertex_atom.coords-edge_overlap_atom.coords
			dummy_atom = pm.Structure(pm_structure.lattice.matrix,[edge_center_name],[edge_overlap_atom.coords+d_vec*2],coords_are_cartesian=True)
			
			if edge_overlap_atom.is_periodic_image(dummy_atom[0],tolerance=tol):
				edge_overlap_to_add.append(edge_overlap_idx)
		edge_overlap_indices += edge_overlap_to_add

		#Make sure the right number of edges are detected
		if len(edge_overlap_indices) != cn:
			warnings.warn('Error: '+topology+'. Incorrect number of edges',Warning)
			pm_structure.to(filename=os.path.join('templates_errors',topology+'.cif'))
			bad = True
			break

		#Find all vertices bound to vertex j
		vertex_overlap_indices = []
		bonded_set = []
		for bond_dist in bond_dists:

			#Get all vertices in bond_dist to vertex j
			vertices_shell_temp = pm_structure.get_neighbors_in_shell(pm_structure[vertex_idx].coords,bond_dist,2*tol,include_index=True,include_image=True)	
			vertices_shell = [k for k in vertices_shell_temp if k[0].species_string != edge_center_name]
			vertices_shell_indices = [k[2] for k in vertices_shell]

			#For each edge connected to vertex j
			for edge_overlap_idx in edge_overlap_indices:

				#Get all vertices in bond_dist/2 to edge center that are also within bond_dist to vertex j
				vertices_shell2_temp = pm_structure.get_neighbors_in_shell(pm_structure[edge_overlap_idx].coords,bond_dist/2,tol,include_index=True)
				vertices_shell2_indices = [k[2] for k in vertices_shell2_temp if k[0].species_string != edge_center_name and k[2] in vertices_shell_indices]

				#For every bonded vertex, get properties
				for vertex2_idx in vertices_shell2_indices:

					#Get properties of the bonded vertex, but don't double-count entries
					locs = np.where(vertices_shell_indices==vertex2_idx)[0].tolist()
					for loc in locs:
						bonded_vertex = vertices_shell[loc]
						img_temp = bonded_vertex[3].tolist()
						img = [int(ii) for ii in img_temp]
						bonded_set_temp = [vertex_idx,bonded_vertex[2],img]
						if bonded_set_temp not in bonded_set:
							bonded_set.append(bonded_set_temp)
							vertex_overlap_indices.append(bonded_vertex[2])
							img_list.append(img)
							d_list.append(pm_structure.get_distance(vertex_idx,bonded_vertex[2],jimage=img))
							break
						else:
							continue

		#Check coordination number
		if len(vertex_overlap_indices) != cn:
			warnings.warn('Error: '+topology+'. Incorrect number of bonded vertices',Warning)
			pm_structure.to(filename=os.path.join('templates_errors',topology+'.cif'))
			bad = True
			break

		#Add set of bonded pair of indices to list
		for vertex_overlap_idx in vertex_overlap_indices:
			bonded_pairs.append([vertex_idx,vertex_overlap_idx])
		
	if bad:
		continue

	#Make the bonding text for the CIF
	bond_text = 'loop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n_geom_bond_distance\n_geom_bond_site_symmetry_2\n_ccdc_geom_bond_type\n'

	#For every bonded pair, get bonding/symmetry info
	done_dot_indices = [] #completed bond pairs with . symmetry
	for j, bonded_pair in enumerate(bonded_pairs):

		#Get distance/image properties
		atom1 = pm_structure[bonded_pair[0]] #Vertex1
		atom2 = pm_structure[bonded_pair[1]] #Vertex2
		output_indices = [atom1.index+1,atom2.index+1] #indices to write in CIF
		img = img_list[j]
		d = d_list[j]
		reverse = False
		if output_indices[0] > output_indices[1]:
			reverse = True
			img = [-1*k for k in img]

		#Make symmetry text
		if img == [0,0,0]:
			symmetry_sym = '.'
		else:
			symmetry_sym = '1_'+str(img[0]+5)+str(img[1]+5)+str(img[2]+5)			

		#Complete bond text string
		if symmetry_sym == '.' and (output_indices in done_dot_indices or [output_indices[1],output_indices[0]] in done_dot_indices):
			continue
		if reverse:
			bond_text += atom2.species_string+str(output_indices[1])+'     '+atom1.species_string+str(output_indices[0])+'    '+str(np.round(d,3))+'   '+symmetry_sym+'     S\n'
		else:
			bond_text += atom1.species_string+str(output_indices[0])+'     '+atom2.species_string+str(output_indices[1])+'    '+str(np.round(d,3))+'   '+symmetry_sym+'     S\n'
		if symmetry_sym == '.':
			done_dot_indices.append(output_indices)
			
	#Write the topology CIF
	with open(os.path.join('templates_database',topology+'.cif'),'w') as w:
		w.write(top_text+cellpar_text+pos_text+bond_text)

	print('Success: '+topology)