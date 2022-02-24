void atom_element() {
    /* reads atomic structure and assigns type 4 atoms in overlapping FEM region to
    elements, and type 2 atoms to nodes. Type 4 atoms are always moved with the FE
    solution, type 2 atoms define the boundary conditions for inner nodes
    
    THIS SHOULD BE DONE WHEN ATOM TYPES ARE DEFINED!!!
    
    Attributes
    ----------
    a_ele_ID
    interaction_atom_element: element number, unique atom number in IMD file
    n_type2: number of type 2 atoms
    aID: unique atom ID and position of type 2 atoms
    
    */

	ifstream inputfile_atom;
	inputfile_atom.open("temp/relaxed_perfect_crystal_with_atom_type.imd");

	if (inputfile_atom == NULL) {
		cout << "File 'relaxed_perfect_crystal_with_atom_type.imd' does not exist!" << endl;
		exit(EXIT_FAILURE);
	}

	double vx, vy, vz, Epot, eam_rho;
	int number, type;
	double mass, x, y, z;
	double fx, fy, fz;

	int count_type;

	count_type = 0;

	double xloc[9], yloc[9], zloc[9];
	double xmin, ymin, zmin, xmax, ymax, zmax;
	double xat, yat, zat;

	double refx, refy, refz; // reference values of first atom
	double diffvalx, diffvaly, diffvalz; // temporary difference values
	int atom_counter; // counts the number of atoms in input file

	bool found;

	ATOM_COUNT = 0;
	n_type2 = 0;

	std::string line;
	char *token;

	atom_counter = 0;
	diffx = Lx;
	diffy = Ly;
	diffz = Lz;
	lxmin = Lx;
	lxmax = 0.0;
	lymin = Ly;
	lymax = 0.0;
    
	while (std::getline(inputfile_atom, line)) {

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

			if (type == 2) {
				if (x < lxmin) {
					lxmin = x;
				}
				if (y < lymin) {
					lymin = y;
				}
				if (x > lxmax) {
					lxmax = x;
				}
				if (y > lymax) {
					lymax = y;
				}
				x = x - Lx * 0.5 + shift_x;
				y = y - Ly * 0.5 + shift_y;
				aID[n_type2][0] = number;
				aID[n_type2][1] = x;
				aID[n_type2][2] = y;
				aID[n_type2][3] = z;
				++n_type2;
				x = x + Lx * 0.5 - shift_x;
				y = y + Ly * 0.5 - shift_y;
			}
			if (atom_counter == 0) {
				refx = x;
				refy = y;
				refz = z;
			}
			else {
				diffvalx = fabs(refx-x);
				if (diffvalx > EPS && diffvalx < diffx) {
					diffx = diffvalx;
				}
				diffvaly = fabs(refy-y);
				if (diffvaly > EPS && diffvaly < diffy) {
					diffy = diffvaly;
				}
				if (diffvalx + diffvaly < EPS) {
					diffvalz = fabs(refz-z);
					if (diffvalz < diffz) {
						diffz = diffvalz;
					}
				}
			}
			atom_counter++;

			x = x - Lx * 0.5 + shift_x;  // Current
			y = y - Ly * 0.5 + shift_y;

			if (type == 4) { //atoms outside the hollow box
				found = false;
				for (int i = 0; i < NEL; i++) {

					xmin = ymin = zmin = 1e16;
					xmax = ymax = zmax = -1e16;

					for (int l = 1; l < 9; l++) {
						xloc[l] = Xglob[LOTOGO[8 * i + (l - 1)]-1];
						yloc[l] = Yglob[LOTOGO[8 * i + (l - 1)]-1];
						zloc[l] = Zglob[LOTOGO[8 * i + (l - 1)]-1];
					}

					for (int j = 1; j < 9; j++) {

						if ((xloc[j] <= xmin) && (yloc[j] <= ymin)
								&& (zloc[j] <= zmin)) {
							xmin = xloc[j];
							ymin = yloc[j];
							zmin = zloc[j];
						}

						if ((xloc[j] >= xmax) && (yloc[j] >= ymax)
								&& (zloc[j] >= zmax)) {
							xmax = xloc[j];
							ymax = yloc[j];
							zmax = zloc[j];
						}
					} // end of loop over the nodes of an element

					if ((xmin >= 1e15) || (ymin >= 1e15) || (xmin >= 1e15)) {
						cout << "Error in finding xmin, ymin and zmin!" << endl;
						exit(EXIT_FAILURE);
					}

					if ((xmax <= -1e15) || (ymax <= -1e15) || (xmax <= -1e15)) {
						cout << "Error in finding xmax, ymax and zmax!" << endl;
						exit(EXIT_FAILURE);
					}

					xat = x;
					yat = y;
					zat = z;

					if ((xat >= xmin) && (yat >= ymin) && (zat >= zmin)
							&& (xat <= xmax) && (yat <= ymax)
							&& (zat <= zmax)) {
						a_ele_ID[ATOM_COUNT][0] = number;
						a_ele_ID[ATOM_COUNT][1] = x;
						a_ele_ID[ATOM_COUNT][2] = y;
						a_ele_ID[ATOM_COUNT][3] = z;
						interaction_atom_element[ATOM_COUNT][0] = i;  // number of element
						interaction_atom_element[ATOM_COUNT][1] = number;  // atom number in IMD file
						//myfile << number << "    " << i << endl;
						++ATOM_COUNT;
						//cout << ATOM_COUNT << endl;

						found = true;
						break;
					}

				} // end of loop over elements

				if (found == false) {
					cout << "Atom number: " << number << endl;
					cout << xat << " " << yat << "  " << endl;
					x = x + Lx * 0.5 - shift_x;  // Current
					y = y + Ly * 0.5 - shift_y;
					cout << x << " " << y << "  " << z << "  " << endl;
					cout << "Error in finding points inside an element!"
							<< endl;
					exit(EXIT_FAILURE);
				}
			} //looping over type 4 atoms

		} //looping over lines which don't have #

	} // End of input while loop on atoms

	cout << "(DIFFX,DIFFY,DIFFZ) IS (" << diffx << "," << diffy << "," << diffz << ")" << endl;
	cout << "(LXMIN,LXMAX) IS (" << lxmin << "," << lxmax << ")" << endl;
	cout << "(LYMIN,LYMAX) IS (" << lymin << "," << lymax << ")" << endl;
	cout << "Number of type 4 atoms found in the mesh: " << ATOM_COUNT << endl;

}

void displacement_interpolation() {
    /* 
    Calculate atomic displacements based on interpolated nodal solution in elements
    in overlapping region. To be applied to type 4 atoms.
    
    Uses a_ele_ID
    
    Attributes
    ----------
    atom_disp
    */

	double SHI, NEW, PHI;

	double xloc[9], yloc[9], zloc[9];
	double xmin, ymin, zmin, xmax, ymax, zmax;
	double xat, yat, zat;

	int ele, atom, inode;
	double lx, ly, lz;
	double dx, dy, dz;

	bool found;
	int n[9];
	double DP[25], SF[9];
	double dax, day, daz;

	double arctan;
	double w;

	int check_count_up, check_count_low, check_count;

	check_count = 0;
	check_count_up = 0;
	check_count_low = 0;

	for (int l = 0; l < ATOM_COUNT; l++) {
		ele = interaction_atom_element[l][0];
		atom = interaction_atom_element[l][1];

		xmin = ymin = zmin = 1e16;
		xmax = ymax = zmax = -1e16;

		for (int l = 1; l < 9; l++) {
			xloc[l] = Xglob[LOTOGO[8 * ele + (l - 1)]-1];
			yloc[l] = Yglob[LOTOGO[8 * ele + (l - 1)]-1];
			zloc[l] = Zglob[LOTOGO[8 * ele + (l - 1)]-1];
		}

		for (int j = 1; j < 9; j++) {

			if ((xloc[j] <= xmin) && (yloc[j] <= ymin) && (zloc[j] <= zmin)) {
				xmin = xloc[j];
				ymin = yloc[j];
				zmin = zloc[j];
			}

			if ((xloc[j] >= xmax) && (yloc[j] >= ymax) && (zloc[j] >= zmax)) {
				xmax = xloc[j];
				ymax = yloc[j];
				zmax = zloc[j];
			}
		} // end of loop over the nodes of an element

		if ((xmin >= 1e15) || (ymin >= 1e15) || (xmin >= 1e15)) {
			cout << "Displacement error in finding xmin, ymin and zmin!"
					<< endl;
			exit(EXIT_FAILURE);
		}

		if ((xmax <= -1e15) || (ymax <= -1e15) || (xmax <= -1e15)) {
			cout << "Displacement error in finding xmax, ymax and zmax!"
					<< endl;
			exit(EXIT_FAILURE);
		}

		bool found1, found2, found3, found4, found5, found6, found7, found8;

		found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 =
				false;

		for (int i = 1; i < 9; i++) {

			if ((xloc[i] == xmin) && (yloc[i] == ymin) && (zloc[i] == zmax)) {
				n[1] = i;
				found1 = true;
			}

			if ((xloc[i] == xmin) && (yloc[i] == ymin) && (zloc[i] == zmin)) {
				n[2] = i;
				found2 = true;
			}

			if ((xloc[i] == xmin) && (yloc[i] == ymax) && (zloc[i] == zmin)) {
				n[3] = i;
				found3 = true;
			}

			if ((xloc[i] == xmin) && (yloc[i] == ymax) && (zloc[i] == zmax)) {
				n[4] = i;
				found4 = true;
			}

			if ((xloc[i] == xmax) && (yloc[i] == ymin) && (zloc[i] == zmax)) {
				n[5] = i;
				found5 = true;
			}

			if ((xloc[i] == xmax) && (yloc[i] == ymin) && (zloc[i] == zmin)) {
				n[6] = i;
				found6 = true;
			}

			if ((xloc[i] == xmax) && (yloc[i] == ymax) && (zloc[i] == zmin)) {
				n[7] = i;
				found7 = true;
			}

			if ((xloc[i] == xmax) && (yloc[i] == ymax) && (zloc[i] == zmax)) {
				n[8] = i;
				found8 = true;
			}
		}

		if (!found1 || !found2 || !found3 || !found4 || !found5 || !found6
				|| !found7 || !found8) {
			cout << "Error in node arrangement!" << endl;
			exit(EXIT_FAILURE);
		}

		lx = fabs(xmax - xmin);
		ly = fabs(ymax - ymin);
		lz = fabs(zmax - zmin);

		found = false;

		for (int j = 0; j < ATOM_COUNT; j++) {
			if (a_ele_ID[j][0] == atom) {
				xat = a_ele_ID[j][1];
				yat = a_ele_ID[j][2];
				zat = a_ele_ID[j][3];
				found = true;
				check_count = check_count + 1;
				break;
			}
		}

		if (found == false) {
			cout << "Atom not found!" << endl;
			exit(EXIT_FAILURE);
		}

		dx = fabs(xat - xmin);
		dy = fabs(yat - ymin);
		dz = fabs(zat - zmin);

		if ((dx > lx) || (dy > ly) || (dz > lz)) {
			cout << "Error in calculating dx, dy and dz!" << endl;
			exit(EXIT_FAILURE);
		}

		SHI = ((2.0 * xat) - (xmin + xmax)) / lx;
		NEW = ((2.0 * yat) - (ymin + ymax)) / ly;
		PHI = ((2.0 * zat) - (zmin + zmax)) / lz;

		if ((SHI > 1.0) || (SHI < -1.0) || (NEW > 1.0) || (NEW < -1.0)
				|| (PHI > 1.0) || (PHI < -1.0)) {
			cout << SHI << "  " << NEW << "  " << PHI << endl;
			cout << xat << " " << yat << "  " << zat << endl;
			cout << xmin << " " << ymin << "  " << zmin << "   " << xmax << " "
					<< ymax << "   " << zmax << endl;
			cout << lx << " " << ly << "  " << lz << "  " << dx << " " << dy
					<< "  " << dz << endl;
			cout << "Error in calculating SHI, NEW and PHI!" << endl;
			exit(EXIT_FAILURE);
		}

		for (int i = 1; i < 9; i++) {
		    inode = LOTOGO[8 * ele + n[i] - 1] - 1;
			DP[3*i-2] = node_dis[3*inode + (1 - 1)];
			DP[3*i-1] = node_dis[3*inode + (2 - 1)];
			DP[3*i  ] = node_dis[3*inode + (3 - 1)];
		}

		SF[1] = 0.125 * (1.0 - SHI) * (1.0 - NEW) * (1.0 + PHI);
		SF[2] = 0.125 * (1.0 - SHI) * (1.0 - NEW) * (1.0 - PHI);
		SF[3] = 0.125 * (1.0 - SHI) * (1.0 + NEW) * (1.0 - PHI);
		SF[4] = 0.125 * (1.0 - SHI) * (1.0 + NEW) * (1.0 + PHI);
		SF[5] = 0.125 * (1.0 + SHI) * (1.0 - NEW) * (1.0 + PHI);
		SF[6] = 0.125 * (1.0 + SHI) * (1.0 - NEW) * (1.0 - PHI);
		SF[7] = 0.125 * (1.0 + SHI) * (1.0 + NEW) * (1.0 - PHI);
		SF[8] = 0.125 * (1.0 + SHI) * (1.0 + NEW) * (1.0 + PHI);

		dax = 0.0;
		day = 0.0;
		daz = 0.0;

		for (int i = 1; i < 9; i++) {
			dax = dax + SF[i] * DP[(3 * i - 2)];
			day = day + SF[i] * DP[(3 * i - 1)];
			daz = daz + SF[i] * DP[3 * i];
		}

		atom_disp[atom][0] = dax;
		atom_disp[atom][1] = day;
		atom_disp[atom][2] = daz;
		
		if ((xat == 0) && (yat == 0)) {
			xat = 0.25;
			arctan = atan2(double(yat), double(xat));
		} else
			arctan = atan2(double(yat), double(xat));

		w = (bv / (2.0 * PI)) * arctan;

		if ((Xglob[LOTOGO[8 * ele + (1 - 1)]-1] <= 0.0)
				&& (Xglob[LOTOGO[8 * ele + (2 - 1)]-1] <= 0.0)
				&& (Xglob[LOTOGO[8 * ele + (3 - 1)]-1] <= 0.0)
				&& (Xglob[LOTOGO[8 * ele + (4 - 1)]-1] <= 0.0)
				&& (Xglob[LOTOGO[8 * ele + (5 - 1)]-1] <= 0.0)
				&& (Xglob[LOTOGO[8 * ele + (6 - 1)]-1] <= 0.0)
				&& (Xglob[LOTOGO[8 * ele + (7 - 1)]-1] <= 0.0)
				&& (Xglob[LOTOGO[8 * ele + (8 - 1)]-1] <= 0.0)) {
			bool condition1, condition2;

			condition1 = false;
			condition2 = false;

			if ((Yglob[LOTOGO[8 * ele + (1 - 1)]-1] < 0.0)
					|| (Yglob[LOTOGO[8 * ele + (2 - 1)]-1] < 0.0)
					|| (Yglob[LOTOGO[8 * ele + (3 - 1)]-1] < 0.0)
					|| (Yglob[LOTOGO[8 * ele + (4 - 1)]-1] < 0.0)
					|| (Yglob[LOTOGO[8 * ele + (5 - 1)]-1] < 0.0)
					|| (Yglob[LOTOGO[8 * ele + (6 - 1)]-1] < 0.0)
					|| (Yglob[LOTOGO[8 * ele + (7 - 1)]-1] < 0.0)
					|| (Yglob[LOTOGO[8 * ele + (8 - 1)]-1] < 0.0))
				condition1 = true;

			if ((Yglob[LOTOGO[8 * ele + (1 - 1)]-1] > 0.0)
					|| (Yglob[LOTOGO[8 * ele + (2 - 1)]-1] > 0.0)
					|| (Yglob[LOTOGO[8 * ele + (3 - 1)]-1] > 0.0)
					|| (Yglob[LOTOGO[8 * ele + (4 - 1)]-1] > 0.0)
					|| (Yglob[LOTOGO[8 * ele + (5 - 1)]-1] > 0.0)
					|| (Yglob[LOTOGO[8 * ele + (6 - 1)]-1] > 0.0)
					|| (Yglob[LOTOGO[8 * ele + (7 - 1)]-1] > 0.0)
					|| (Yglob[LOTOGO[8 * ele + (8 - 1)]-1] > 0.0))
				condition2 = true;

			if (condition1 && condition2) {
				if (yat >= 0) {
					atom_disp[atom][2] = w;		//DISP[LOTOGO[ele][n[3]]][1];

				} else {
					atom_disp[atom][2] = -w;		//DISP[LOTOGO[ele][n[1]]][1];
				}
			}
		} // end of loop on xfem elements

	} // end of for loop on 'l' atoms

	cout << "Number of times atom found: " << check_count << endl;
}

void init_atom_config() {
    /* Produce IMD input file with atomic shifts from FEM solution in boundary region
    applied to type 4 atoms. All other atom types are shifted according to Volterra 
    solution for screw dislocation.
    
    */
	ifstream inputfile_atom;
	inputfile_atom.open("temp/relaxed_perfect_crystal_with_atom_type.imd");

	if (inputfile_atom == NULL) {
		cout << "File 'relaxed_perfect_crystal_with_atom_type.imd' does not exist!" << endl;
		exit(EXIT_FAILURE);
	}

	ofstream myfile;
	myfile.open("temp/atomistic_dislocation_with_fem_solution.imd");

	double vx, vy, vz, Epot, eam_rho;
	int number, type;
	double mass, x, y, z;
	double fx, fy, fz;

	double arctan;
	double w;

	int type_atom;

	type_atom = 0;

	std::string line;
	char *token;

	while (std::getline(inputfile_atom, line)) {

		if (line[0] == '#' && line[1] == 'F')
			myfile << "#F A 1 1 1 3 0 0" << endl;

		if (line[0] == '#' && line[1] == 'C')
			myfile << "#C number type mass x y z" << endl;

		if (line[0] == '#' && line[1] == 'X')
			myfile << line << endl;

		if (line[0] == '#' && line[1] == 'Y')
			myfile << line << endl;

		if (line[0] == '#' && line[1] == 'Z')
			myfile << line << endl;

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

			x = x - Lx * 0.5 + shift_x;  // Current
			y = y - Ly * 0.5 + shift_y;

			if (type == 4) {

				z += atom_disp[number][2];
				y += atom_disp[number][1];
				x += atom_disp[number][0];

				if (atom_disp[number][2] == 0.0) {
					cout << "Error in atom displacement calculations!" << endl;
					exit(EXIT_FAILURE);
				}

			} else {

				if ((x == 0) && (y == 0)) {
					x = 0.25;
					arctan = atan2(double(y), double(x));
				} else
					arctan = atan2(double(y), double(x));

				w = (bv / (2.0 * PI)) * arctan;

				z += w;

			}

			x = x + Lx * 0.5 - shift_x;  // Current
			y = y + Ly * 0.5 - shift_y;

			myfile << number << "  " << type << "  " << mass << "  " << x
					<< "  " << y << "  " << z << endl;

		} //looping over lines which don't have #

	} // End of input while loop on atoms

	//cout << "Total type 2 atoms inside the system box : " << type_atom << endl;

}

void atom_configuration() {
    /*  Apply nodal positons to type 4 atoms and write updated atomic
        configuration. Other atom types remain untouched.
        
        Reads atomistic configuration from file: 
            relaxed_atomistic_dislocation_structure.00000.ss
        Produces file:
            atomistic_dislocation_with_fem_solution.imd
    
    */
	ifstream inputfile_atom;
	inputfile_atom.open("temp/relaxed_atomistic_dislocation_structure.00000.ss");

	if (inputfile_atom == NULL) {
		cout << "File relaxed_atomistic_dislocation_structure.00000.ss does not exist!" << endl;
		exit(EXIT_FAILURE);
	}

	ofstream myfile;
	myfile.open("temp/atomistic_dislocation_with_fem_solution.imd");

	double vx, vy, vz, Epot, eam_rho;
	int number, type;
	double mass, x, y, z;
	double fx, fy, fz;

	std::string line;
	char *token;

	while (std::getline(inputfile_atom, line)) {

		if (line[0] == '#' && line[1] == 'F')
			myfile << "#F A 1 1 1 3 0 0" << endl;

		if (line[0] == '#' && line[1] == 'C')
			myfile << "#C number type mass x y z" << endl;

		if (line[0] == '#' && line[1] == 'X')
			myfile << line << endl;

		if (line[0] == '#' && line[1] == 'Y')
			myfile << line << endl;

		if (line[0] == '#' && line[1] == 'Z')
			myfile << line << endl;

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

			if (type == 4) {

				z = cor_type4[number][2] + atom_disp[number][2];
				y = cor_type4[number][1] + atom_disp[number][1];
				x = cor_type4[number][0] + atom_disp[number][0];

				if (atom_disp[number][2] == 0.0) {
					cout << number << " " << atom_disp[number][2] << "  " << x
							<< "  " << y << "  " << z << endl;
					cout << "Error in atom displacement calculations!" << endl;
					exit(EXIT_FAILURE);
				}
			} // end of input while loop

			myfile << number << "  " << type << "  " << mass << "  " << x
					<< "  " << y << "  " << z << endl;

		} //looping over lines which don't have #

	} // End of input while loop on atoms

}

void nodal_displacement() {
    /* Get displacement of relaxed atoms and apply as BC to XFEM model
        Reads relaxed atomistic configuration from file: 
            relaxed_atomistic_dislocation_structure.00000.ss
        Atomic displacements are calculated relative to:
            relaxed_perfect_crystal_with_atom_type.imd
    
    */

	ifstream inputfile_atom_disp;
	inputfile_atom_disp.open("temp/relaxed_atomistic_dislocation_structure.00000.ss");

	if (inputfile_atom_disp == NULL) {
		cout << "File 'relaxed_atomistic_dislocation_structure' does not exist!" << endl;
		exit(EXIT_FAILURE);
	}

	ifstream inputfile_atom;
	inputfile_atom.open("temp/relaxed_perfect_crystal_with_atom_type.imd");
	/* THIS INFORMATION SHOULD BE STORED WHEN GENERATING THE FILE:
	data_disp_un */
	if (inputfile_atom == NULL) {
		cout << "File 'relaxed_perfect_crystal_with_atom_type.imd' does not exist!" << endl;
		exit(EXIT_FAILURE);
	}

	double vx, vy, vz, Epot, eam_rho;
	int number, type;
	double mass, x, y, z;
	double fx, fy, fz;

	std::string line;
	char *token;

	while (std::getline(inputfile_atom_disp, line)) {

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

			data_disp[number][0] = x;
			data_disp[number][1] = y;
			data_disp[number][2] = z;

		} //looping over lines which don't have #

	} //End of input while loop on atoms

	std::string line_string;
	char *token_string;

	while (std::getline(inputfile_atom, line_string)) {

		if (line_string[0] != '#') {

			token_string = strtok(&line_string[0], " ");
			number = atoi(token_string);

			int ntokens = 1;

			while ((token_string = strtok(NULL, " ")) != NULL) {

				ntokens = ntokens + 1;

				switch (ntokens) {
				case 2:
					type = atoi(token_string);
					break;
				case 3:
					mass = strtod(token_string, NULL);
					break;
				case 4:
					x = strtod(token_string, NULL);
					break;
				case 5:
					y = strtod(token_string, NULL);
					break;
				case 6:
					z = strtod(token_string, NULL);
					break;
				}
			}

			data_disp_un[number][0] = x;
			data_disp_un[number][1] = y;
			data_disp_un[number][2] = z;

		} //looping over lines which don't have #

	} //End of input while loop on atoms

	double dx, dy, dz;

    int ind, inode;
	for (int i = 0; i < NCONSNODE; i++) {
	    ind = interaction_atom_node[i][1];
		dz = data_disp[ind][2] - data_disp_un[ind][2];
		dy = data_disp[ind][1] - data_disp_un[ind][1];
		dx = data_disp[ind][0] - data_disp_un[ind][0];

		if (fabs(dz) > (Lz * 0.5)) {
			if (dz > 0.0)
				dz = dz - Lz;
			else
				dz = dz + Lz;
		}
		
		inode = interaction_atom_node[i][0];
		ENFRDISPglob[3*inode    ] = dx;
		ENFRDISPglob[3*inode + 1] = dy;
	    ENFRDISPglob[3*inode + 2] = dz;
	}

}

void atom_node() {
    /* Assign atoms of type 2 to nodes on inner boundary nodes
    
    Attributes
    ----------
    interaction_atom_node
    */

	double min, min_tol, dis;
	double x1, x2, z1, y1, y2, z2;

	for (int i = 0; i < NCONSNODE; i++) {
		min_tol = 1.0;

		x1 = nID[i][1];
		y1 = nID[i][2];
		z1 = nID[i][3];

		for (int j = 0; j < n_type2; j++) {

			if (z1 > (Lz - 0.05*diffz)) //change this!!!!!!!!!
				z1 = z1 - Lz;

			x2 = aID[j][1];
			y2 = aID[j][2];
			z2 = aID[j][3];

			dis = pow(
					pow((x1 - x2), 2.0) + pow((y1 - y2), 2.0)
							+ pow((z1 - z2), 2.0), 0.5);

			min = dis;

			if (min < min_tol) {
				min_tol = min;
				interaction_atom_node[i][0] = nID[i][0];
				interaction_atom_node[i][1] = aID[j][0];
			}

		} // end of j

		if (min_tol >= 1.0) {
			cout << "Error in finding atom node pair! " << min << endl;
			exit(EXIT_FAILURE);
		}

	} // end of i
}

void create_volterra_dis() {
    /* create BC for Volterra dislocation on constraint nodes
    
    Parameters:

    */
    int i, inode;
    double yc;
    
    // outer boundary
    for (i = 0; i < NENFD; i++) {
		inode = ncstr_o[i];
		ENFRDISPglob[3*inode    ] = 0.0;
		ENFRDISPglob[3*inode + 1] = 0.0;
		ENFRDISPglob[3*inode + 2] = bv * atan2(Yglob[inode-1], Xglob[inode-1])/(2.0*PI);
	}	
	//inner boundary
	for (i = 0; i < NCONSNODE; i++) {
	    inode = ncstr_i[i];
		ENFRDISPglob[3*inode    ] = 0.0;
		ENFRDISPglob[3*inode + 1] = 0.0;
	    ENFRDISPglob[3*inode + 2] = bv * atan2(Yglob[inode-1], Xglob[inode-1])/(2.0*PI);
	}
}

