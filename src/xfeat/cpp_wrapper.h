#ifndef HEADER_XFEAT_
#define HEADER_XFEAT_

#define NDLOC 25  // Local dof K is 8 x 3 + 1
#define NODEN 9 //number of nodes per element + 1
#define ND 4  // 3D + 1
#define NDS 7 // Number of stress components 6 + 1
#define NELEMAX 40000  // size of number of elements
#define NNODEMAX 40000   // size of number of nodes
#define NELEIN 5000  // size of number of nodes/elements/atoms constrained
#define ATOMTOTAL 20000
#define EPS 0.00000000001 // general precision value for the program

// definitions for variables shared with Cython part
// constants and model parameters
double PI;
std::string temp_dir;
double fem_size;
double EE3D[NDS][NDS];

// atomistic quantities
double bv;
double shift[2];

// XFEM quantities
int NNODE, NLAG, NEL, NDF, xdof;
int nelem_full_big_box_x;
int nelem_full_big_box_y;
int nelem_full_big_box_z;
double Lx, Ly, Lz;
double sigma[3];
double at_disp[ATOMTOTAL][3];
std::vector<int> Iglob;
std::vector<int> Jglob;
std::vector<int> JJglob;
std::vector<int> LOTOGO;
std::vector<double> NEWAglob;
std::vector<double> Dglob;
std::vector<double> Fglob;
std::vector<double> ENFRDISPglob;
std::vector<double> Xglob;
std::vector<double> Yglob;
std::vector<double> Zglob;
std::vector<double> node_dis;

// shared functions
void atom_set_up();
void create_mesh();
void init_matrix();
void update_fglob();
void atom_node();
void atom_element();
void displacement_interpolation();
void atom_configuration();
void init_screw_dis();
void nodal_displacement();
void create_volterra_dis();
void apply_e23_outer(double e23);
void STRESS(double XLOC[NODEN], double YLOC[NODEN], double ZLOC[NODEN], int IEL);

// other global variables
int NCF, NENFD, NCONT;
int n_type2; // number of type 2 atoms
int n_type4; // number of type 4 atoms
int NCONSNODE; // number of nodes on the inner surface
int num_space_x = 0;
int num_space_y = 0;
int ctrmat = 0;
int natom = 0;
int fem_elems_round[2];
int interaction_atom_element[ATOMTOTAL][2];
int interaction_atom_node[ATOMTOTAL][2];
int atom_id[ATOMTOTAL];
int ind_type4[ATOMTOTAL];  // store numbers of type 4 atoms


double Det;
double shift_x, shift_y;
double diffx, diffy, diffz; // space between planes in each dimension (diffz is approx. the Burgers vector?)
double lxmin, lxmax, lymin, lymax; // limit coordinate values for type 2 atoms
double minwidth_x = 11.0;
double minwidth_y = 11.0;
double eps = 1.e-7;
double dist[3];
double atom_disp[ATOMTOTAL][3];
double coords[ATOMTOTAL][3];
double masses[ATOMTOTAL];
double nID[NELEIN][4];
double aID[ATOMTOTAL][4];
double a_ele_ID[ATOMTOTAL][4];
double ref_coord_high[3];
double ref_coord_low[3];
double sys_width[2];
double sys_reduced[2];
double fem_big_box[3];
double Ajac[ND][ND];   //jacobian matrix
double AjacInv[ND][ND]; //jacobian inverse matrix
double AKLOC[NDLOC][NDLOC]; //Local stiffness matrix
double ENFLOC[NDLOC];
double STR[NDS];
double DISLOCATIONF[NDLOC];

std::vector<int> NODE;
std::vector<int> ContraintNodes;  // nodes with Lagrange constraints along slip plane
std::vector<int> ncstr_o;  // list of constraint nodes on outer boundary
std::vector<int> ncstr_i;  // list of constraint nodes on inner boundary
std::vector<double> EEglob;
std::vector<double> FJointglob;
std::vector<double> Aglob;
std::vector<std::map<int, double> > MAT;

// other functions
int EqNum(int * JJglob, int Njoint);
void STIFF(double XLOC[NODEN], double YLOC[NODEN], double ZLOC[NODEN]);
void Assem(int Iel);
void JacInv();
void skylineTOcompressedarray();

#endif  /*  HEADER_XFEAT_  */
