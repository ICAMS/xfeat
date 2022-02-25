extern void atom_set_up(){
    /* Read atomic structure and set atoms types to which boundary 
    conditions will be applied later. Write out new IMD input file with types.
    Definitions:
    type 0: standard atoms w/o restrictions
    type 2: interface atoms between free and restricted layers,
            move freely during relaxation and serve as BC for inner XFEM boundary nodes
    type 4: restricted atoms that are coupled to XFEM distortions, do not move
            during relaxation
            
    reads: relaxed_perfect_crystal.imd
    creates: relaxed_perfect_crystal_with_atom_type.imd
    
    Attributes
    ----------
    atom_id
    coords
    masses
    natom
    n_type2
    n_type4
    aID
    ind_type4
    
    */
	ifstream inputfile;
	inputfile.open(temp_dir + "/relaxed_perfect_crystal.imd");

	if (inputfile == NULL) {
		cout << "File 'relaxed_perfect_crystal.imd' does not exist!" << endl;
		exit(EXIT_FAILURE);
	}

	ofstream myfile;
	myfile.open(temp_dir + "/relaxed_perfect_crystal_with_atom_type.imd");

	int number, type;
	double mass, x, y, z;

	int count = 0;

	std::string line;
	char *token;

	double ref_atom[3];
	double plane_dist[3];
	double min_coord[3];
	int finished;

    cout.precision(12);

	while (std::getline(inputfile, line)) {

		if (line[0] == '#' && line[1] == 'F')
			myfile << "#F A 1 1 1 3 0 0" << endl;

		if (line[0] == '#' && line[1] == 'C')
			myfile << "#C number type mass x y z" << endl;

		if (line[0] == '#' && line[1] == 'X') {

			myfile << line << endl;

			int ntokens = 1;

			token = strtok(&line[0], " ");

			while ((token = strtok(NULL, " ")) != NULL) {

				ntokens = ntokens + 1;

				switch (ntokens) {
				case 2:
					Lx = strtod(token, NULL);
					break;
				}
			}
		}

		if (line[0] == '#' && line[1] == 'Y') {

			myfile << line << endl;

			int ntokens = 1;

			token = strtok(&line[0], " ");

			while ((token = strtok(NULL, " ")) != NULL) {

				ntokens = ntokens + 1;

				switch (ntokens) {
				case 3:
					Ly = strtod(token, NULL);
					break;
				}
			}
		}

		if (line[0] == '#' && line[1] == 'Z') {

			myfile << line << endl;

			int ntokens = 1;

			token = strtok(&line[0], " ");

			while ((token = strtok(NULL, " ")) != NULL) {

				ntokens = ntokens + 1;

				switch (ntokens) {
				case 4:
					Lz = strtod(token, NULL);
					break;
				}
			}
		}

		if (line[0] == '#' && line[1] == 'E')
			myfile << line << endl;

		if (line[0] != '#') {

			token = strtok(&line[0], " ");
			number = atoi(token);

			int ntokens = 1;

			while ((token = strtok(NULL, " ")) != NULL) {

				ntokens = ntokens + 1;

				switch (ntokens) {
				case 2:
					type = atoi(token);
					break;
				case 3:
					mass = strtod(token, NULL);
					break;
				case 4:
					x = strtod(token, NULL);
					break;
				case 5:
					y = strtod(token, NULL);
					break;
				case 6:
					z = strtod(token, NULL);
					break;
				}

			}

			if (natom == 0) {
				ref_atom[0] = x;
				ref_atom[1] = y;
				ref_atom[2] = z;
				dist[0] = Lx;
				dist[1] = Ly;
				dist[2] = Lz;
				plane_dist[0] = Lx;
				plane_dist[1] = Ly;
				plane_dist[2] = Lz;
				min_coord[0] = x;
				min_coord[1] = y;
				min_coord[2] = z;
				bv = 0.0;
			}
			else {
				double diffy = fabs(ref_atom[1] - y);
				if (fabs(x-ref_atom[0]) < eps) {
					if (fabs(z-ref_atom[2]) < eps) {
						if (diffy < dist[1]) {
							dist[1] = diffy;
						}
					}
				}
				if (diffy < plane_dist[1] && diffy > eps) {
					plane_dist[1] = diffy;
				}
				double diffx = fabs(ref_atom[0] - x);
				if (fabs(y-ref_atom[1]) < eps) {
					if (fabs(z-ref_atom[2]) < eps) {
						if (diffx < dist[0]) {
							dist[0] = diffx;
						}
					}
				}
				if (diffx < plane_dist[0] && diffx > eps) {
					plane_dist[0] = diffx;
				}
				double diffz = fabs(ref_atom[2] - z);
				if (fabs(x-ref_atom[0]) < eps) {
					if (fabs(y-ref_atom[1]) < eps) {
						if (diffz < dist[2]) {
							dist[2] = diffz;
						}
					}
				}
				if (diffz < plane_dist[2] && diffz > eps) {
					plane_dist[2] = diffz;
				}
				if (min_coord[0] > x) min_coord[0] = x;
				if (min_coord[1] > y) min_coord[1] = y;
				if (min_coord[2] > z) min_coord[2] = z;
        		bv = dist[2];
			}

			atom_id[natom] = number;
			coords[number][0] = x;
			coords[number][1] = y;
			coords[number][2] = z;
			masses[natom] = mass;

			natom++;

		} //looping over lines which don't have #

	} // End of input while loop on atoms

	if (minwidth_x > (Lx/2 - plane_dist[0])) {
		cout << " 'minwidth_x' (" << minwidth_x << ") too large for 'Lx' (" << Lx << ")" << endl;
		exit(-1);
	}
	if (minwidth_y > (Ly/2 - plane_dist[1])) {
		cout << " 'minwidth_y' (" << minwidth_y << ") too large for 'Ly' (" << Ly << ")" << endl;
		exit(-1);
	}

	ref_coord_low[0] = min_coord[0];
	ref_coord_low[1] = min_coord[1];
	ref_coord_low[2] = min_coord[2];

	finished = 0;
	while (!finished) {
		if (ref_coord_low[0]+dist[0] >= minwidth_x) {
			finished = 1;
		}
		ref_coord_low[0] += dist[0];
	}

	ref_coord_high[0] = ref_coord_low[0];
	finished = 0;
	while (!finished) {
		if (ref_coord_high[0]+dist[0] >= Lx-minwidth_x) {
			finished = 1;
		}
		if (!finished) {
			ref_coord_high[0] += dist[0];
			num_space_x++;
		}
	}
	finished = 0;
	while (!finished) {
		if (ref_coord_low[1]+dist[1] >= minwidth_y) {
			finished = 1;
		}
		ref_coord_low[1] += dist[1];
	}

	ref_coord_high[1] = ref_coord_low[1];
	finished = 0;
	while (!finished) {
		if (ref_coord_high[1]+dist[1] >= Ly-minwidth_y) {
			finished = 1;
		}
		if (!finished) {
			ref_coord_high[1] += dist[1];
			num_space_y++;
		}
		if (finished && (num_space_y % 2 == 0)) {
			num_space_y--;
			ref_coord_high[1] -= dist[1];
		}
	}
	
    // assign atom types
    n_type4 = 0;
	n_type2 = 0;
	for (count = 0; count < natom; count++) {
	    x = coords[atom_id[count]][0];
	    y = coords[atom_id[count]][1];
	    z = coords[atom_id[count]][2];
		if (((fabs(x-ref_coord_low[0]) < 0.2*plane_dist[0]) && (y >= ref_coord_low[1]-0.2*plane_dist[1]) && (y <= ref_coord_high[1]+0.2*plane_dist[1])) ||
		    ((fabs(y-ref_coord_low[1]) < 0.2*plane_dist[1]) && (x >= ref_coord_low[0]-0.2*plane_dist[0]) && (x <= ref_coord_high[0]+0.2*plane_dist[0])) ||
		    ((fabs(x-ref_coord_high[0]) < 0.2*plane_dist[0]) && (y >= ref_coord_low[1]-0.2*plane_dist[1]) && (y <= ref_coord_high[1]+0.2*plane_dist[1])) ||
		    ((fabs(y-ref_coord_high[1]) < 0.2*plane_dist[1]) && (x >= ref_coord_low[0]-0.2*plane_dist[0]) && (x <= ref_coord_high[0]+0.2*plane_dist[0]))) {
			type = 2;
            aID[n_type2][0] = atom_id[count];
            aID[n_type2][1] = x;
            aID[n_type2][2] = y;
            aID[n_type2][3] = z;
            ++n_type2;
		}
		else if (x > ref_coord_low[0] && y > ref_coord_low[1] &&
		    x < ref_coord_high[0] && y < ref_coord_high[1]) {
                        type = 0;
                }
		else {
			type = 4;
            ind_type4[n_type4] = atom_id[count];
            ++ n_type4;
		}
		myfile << atom_id[count] << "  " << type << "  " << masses[count] << "  " << x << "  " << y << "  " << z << endl;
	}

	inputfile.close();
	myfile.close();

	sys_width[0] = ref_coord_high[0]-ref_coord_low[0];
	sys_width[1] = ref_coord_high[1]-ref_coord_low[1];

	double ref_coord_low_shift[2];
	ref_coord_low_shift[0] = ref_coord_low[0] - 0.5*Lx;
	ref_coord_low_shift[1] = ref_coord_low[1] - 0.5*Ly;
	double ref_coord_high_shift[2];
	ref_coord_high_shift[0] = ref_coord_high[0] - 0.5*Lx;
	ref_coord_high_shift[1] = ref_coord_high[1] - 0.5*Ly;
    shift[0] = 0.5*sys_width[0]-ref_coord_high_shift[0];
    shift[1] = -0.5*sys_width[1]-ref_coord_low_shift[1];
    //myfile_shift << shift[0] << "  " << shift[1] << endl;
	
	sys_reduced[0] = (0.5*sys_width[0]-0.05)*2;
	sys_reduced[1] = (0.5*sys_width[1]-0.05)*2;
	double fem_elems[2];
	fem_elems[0] = (ceil(fem_size/dist[0])*dist[0]-sys_reduced[0])/dist[0];
	fem_elems[1] = (ceil(fem_size/dist[1])*dist[1]-sys_reduced[1])/dist[1];

	fem_elems_round[0] = floor(fem_elems[0]);
	if (fem_elems_round[0] % 2 == 1) {
		fem_elems_round[0]++;
	}
	fem_elems_round[1] = floor(fem_elems[1]);
	if (fem_elems_round[1] % 2 == 1) {
		fem_elems_round[1]++;
	}
	fem_big_box[0] = 0.5*(((double)fem_elems_round[0])*dist[0]+sys_reduced[0]);
	fem_big_box[1] = 0.5*(((double)fem_elems_round[1])*dist[1]+sys_reduced[1]);
	fem_big_box[2] = Lz;

    cout << "FEM Inner BOX DIMENTIONS(X,Y) : " << -(0.5*sys_width[0]-0.05)<<"," << -(0.5*sys_width[1]-0.05)  <<"   "<< 0.5*sys_width[0]-0.05<< "," << 0.5*sys_width[1]-0.05  << endl;
	cout << "FEM BIG BOX (X,Y,Z): " << -(fem_big_box[0]) << "," << -(fem_big_box[1]) <<"   "<< fem_big_box[0] << "," << fem_big_box[1] << "   " << fem_big_box[2] << endl;
	cout << "NUMBER OF MESHES OF HALF EDGES IN FEM (X,Y): " << fem_elems_round[0]/2 << "  " << fem_elems_round[1]/2 << endl;
	cout << "NUMBER OF GAPS BETWEEN ATOMS IN MD MODEL AND THE FEM INNER BOX MESH NUMBER (X,Y): " << num_space_x << "  " << num_space_y << endl;
	cout << "dist is (" << dist[0] << "," << dist[1] << "," << dist[2] << ") ; plane_dist is (" << plane_dist[0] << "," << plane_dist[1] << "," << plane_dist[2] << ")" << endl;
	
}
