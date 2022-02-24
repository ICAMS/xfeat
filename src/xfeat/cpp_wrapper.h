#ifndef HEADER_XFEAT_
#define HEADER_XFEAT_

#define NDLOC 25  // Local dof K is 8 x 3 + 1
#define NODEN 9 //number of nodes per element + 1
#define ND 4  // 3D + 1
#define NDS 7 // Number of stress components 6 + 1
#define NELEMAX 40000  // size of number of elements
#define NNODEMAX 40000   // size of number of nodes
#define NELEIN 5000  // size of number of nodes/elements/atoms constrained
#define ATOMTOTAL 40000
#define EPS 0.00000000001 // general precision value for the program

double PI;

// definitions for Atom_Set_up
double fem_size;
double Lx, Ly, Lz;
int natom = 0;
int *atom_id;
double *coords;
double *masses;
double minwidth_x = 11.0;
double minwidth_y = 11.0;
double eps = 1.e-7;

double shift[2];
double dist[3];
double bv; //burgers vector

void get_n_atoms();
void atom_set_up();

//definitions for FEM_Mesh
//int nodecount;
int nelem_full_big_box_x;
int nelem_full_big_box_y;
int nelem_full_big_box_z;

double ref_coord_high[3];
double ref_coord_low[3];
int num_space_x = 0;
int num_space_y = 0;
double sys_width[2];
int fem_elems_round[2];
double sys_reduced[2];
double fem_big_box[3];
//int size_Aglob;
int xdof;

int NDF;  // no. of degrees of freedom
int NNODE, NEL, NCF, NENFD, NCONT;
int NLAG;
int ctrmat = 0;

double Ajac[ND][ND];   //jacobian matrix
double AjacInv[ND][ND]; //jacobian inverse matrix
double Det;
double EE3D[NDS][NDS];
double sigma[3];

double AKLOC[NDLOC][NDLOC]; //Local stiffness matrix
double ENFLOC[NDLOC];
double STR[NDS];
double DISLOCATIONF[NDLOC];

std::vector<std::map<int, double> > MAT;

std::vector<double> Xglob;
std::vector<double> Yglob;
std::vector<double> Zglob;
std::vector<double> Aglob;
std::vector<double> Dglob;
std::vector<double> NEWAglob;
std::vector<double> Fglob;
std::vector<double> EEglob;
std::vector<double> FJointglob;
std::vector<double> ENFRDISPglob;
std::vector<double> node_dis;

std::vector<int> Iglob;
std::vector<int> Jglob;
std::vector<int> JJglob;
std::vector<int> LOTOGO;
std::vector<int> NODE;
std::vector<int> ContraintNodes;  // nodes with Lagrange constraints along slip plane
std::vector<int> ncstr_o;  // list of constraint nodes on outer boundary
std::vector<int> ncstr_i;  // list of constraint nodes on inner boundary

int EqNum(int * JJglob, int Njoint);

void create_mesh();
void init_matrix();
void iter_matrix();
void STIFF(double XLOC[NODEN], double YLOC[NODEN], double ZLOC[NODEN]);
void Assem(int Iel);
void STRESS(double XLOC[NODEN], double YLOC[NODEN], double ZLOC[NODEN],
		int IEL);
void JacInv();
void skylineTOcompressedarray();
void apply_e23_outer(double e23);

//definitions for Init_Dislocation.cpp
int n_type2;
//int INELE; // number of elements on the inner surface
int NCONSNODE; // number of nodes on the inner surface
int ATOM_COUNT; // number of type 2 atoms on the outside of the inner surface
double cor_type4[ATOMTOTAL][3];  // store coordinates of type 4 atoms

//double XA[ATOMTOTAL], YA[ATOMTOTAL], ZA[ATOMTOTAL];

double shift_x, shift_y;

double diffx, diffy, diffz; // space between planes in each dimension (diffz is approx. the Burgers vector?)
double lxmin, lxmax, lymin, lymax; // limit coordinate values for type 2 atoms

double atom_disp[ATOMTOTAL][3];
double nID[NELEIN][4];
double aID[ATOMTOTAL][4];
double a_ele_ID[ATOMTOTAL][4];

int interaction_atom_element[ATOMTOTAL][2];
int interaction_atom_node[ATOMTOTAL][2];

void atom_element();
void displacement_interpolation();
void atom_configuration();
void init_atom_config();
void nodal_displacement();
void create_volterra_dis();
void atom_node();

// definitions for Nodal_Calc.cpp
double data_disp[ATOMTOTAL][3];
double data_disp_un[ATOMTOTAL][3];

// definitions for Iterative_Step.cpp
void displacement_interpolation();
void atom_configuration();


#endif  /*  HEADER_XFEAT_  */