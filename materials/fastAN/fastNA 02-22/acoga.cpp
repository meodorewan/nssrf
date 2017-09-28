// author	: Do Duc Dong, ITI.VNU
//            Dang Cao Cuong, UET.VNU
// date		: 22/2/2014


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <queue>

using namespace std;

typedef pair<int , int> ii;
typedef vector<ii> vii;

const int MAX_NODE    = 10000;

char ver_name[100]			= "An Ant Colony Optimization for Graph Alignment Problem";
char input1_file[100]		= "celeg.ga";
char input2_file[100]		= "dmela.ga";
char output_file[100]		= "";
char output2_file[100]		= "";
char blat_input_file[100]   = "";
char log_file[100]			= "";

ofstream flog;

clock_t startProgram,			// do thoi gian chay chuong trinh
		endProgram,
		startTry,				// do thoi gian mot lan chay
		endTry;

int max_tries,                  // so lan chay
	n_try,
    Nc,                         // so vong lap
	iteration,
	seed;

int max_time;					// thoi gian chay toi da

int n_ants;                     // so kien
double rho;						// tham so bay hoi

double alpha;

int	ls_flag;

int kAlg;
// 0. -pheromone
// 1. +pheromone

struct ant_struct{
  int       tour[MAX_NODE];
  int       mark1[MAX_NODE], mark2[MAX_NODE];
  double	tour_length;
  int       E12;
  int 		nodes_LCCS, edges_LCCS;
  int 	    nn, nodeList[MAX_NODE];
  double	maxFit;
  int       maxRelation;

  int		iteration;
  int		time;
};

ant_struct	best_so_far_ant[100],	    // ket qua tot nhat cua tung lan chay, toi da 100 lan chay
			ibest_ant,					//
			current_ant;				// kien dang duoc xay dung

int		*pheromone_mark[MAX_NODE];
double	*pheromone[MAX_NODE];
double	trail_max, trail_min, trail_mid;

int     *matrixGraph2[MAX_NODE];

struct graph{
  vector<vii> adj;
  int		  n;
  int		  m;
};

graph g1, g2;
double *BlatSimilar[MAX_NODE];            // thong tin heuristic dua tren blat

double randomdouble(){
	int x=rand();
	double y=(double)x /(double) RAND_MAX;
	return y;
}

void parse_commandline(int argc, char *argv[]) {

	 max_tries				= -1;
     Nc						= -1;
	 max_time				= -1;
     n_ants					= -1;
	 alpha					= -1;
	 rho					= -1;
	 ls_flag				= -1;
	 kAlg   				= -1;

     int i;

     for (i=1; i<argc; i=i+2) {

         // -i1 : ten file du lieu vao
         if (strcmp(argv[i],"-i1")==0) {
              strcpy(input1_file,"");
		      strcat(input1_file,argv[i+1]);
         }

         // -i2 : ten file du lieu vao
         if (strcmp(argv[i],"-i2")==0) {
              strcpy(input2_file,"");
		      strcat(input2_file,argv[i+1]);
         }

         // -k : loai thuat toan
         if (strcmp(argv[i],"-k")==0) {
              kAlg = atoi(argv[i+1]);
         }

         // -s : so ngau nhien ban dau
         if (strcmp(argv[i],"-s")==0) {
              seed = atoi(argv[i+1]);
         }
		 else seed = (unsigned)time( NULL );

         // -c : so vong lap
         if (strcmp(argv[i],"-c")==0) {
              Nc = atoi(argv[i+1]);
         }

		 // -t : thoi gian chay
         if (strcmp(argv[i],"-t")==0) {
              max_time = atoi(argv[i+1]);
         }

         // -r : so lan chay
         if (strcmp(argv[i],"-r")==0) {
              max_tries = atoi(argv[i+1]);
         }

         // -a : alpha
         if (strcmp(argv[i],"-a")==0) {
              alpha = atof(argv[i+1]);
         }

         // -e : he so bay hoi
         if (strcmp(argv[i],"-e")==0) {
              rho = atof(argv[i+1]);
         }

		 // -m : so kien
         if (strcmp(argv[i],"-m")==0) {
              n_ants = atoi(argv[i+1]);
         }

		 // -l : local search
         if (strcmp(argv[i],"-l")==0) {
              ls_flag = atoi(argv[i+1]);
         }
     }

     // ten file mac dinh them .log
     strcpy(log_file,input1_file);
     strcat(log_file,"_");
     strcat(log_file,input2_file);
     strcat(log_file,".log");

     // ten file mac dinh them .out
     strcpy(output_file,input1_file);
     strcat(output_file,"_");
     strcat(output_file,input2_file);
     strcat(output_file,".out");

	 strcpy(output2_file,input2_file);
     strcat(output2_file,"_");
     strcat(output2_file,input1_file);
     strcat(output2_file,".out");

	 strcpy(blat_input_file,input1_file);
	 strcat(blat_input_file,"_");
	 strcat(blat_input_file,input2_file);
	 strcat(blat_input_file,".pin");
}

void read_InputFile(char input_file[], graph &g) {
    int i, j;
	ifstream fin;
    fin.open(input_file);

    g.n = 0;
	g.m = 0;
    g.adj.assign(g.n + 1 , vii());

    int u, v, w;
	string s1, s2, snode1, snode2;

	while (fin>>s1>>s2) {

		if (s1=="node") {
			g.n++;
			fin>>s1>>s2;
			fin>>s1>>s2;
			fin>>s1>>s2;
			fin>>s1;
		}

		if (s1=="edge") {
			g.m++;
			fin>>snode1>>u;
			fin>>snode2>>v;
			fin>>s1;

			w = 1;
			g.adj[u].push_back(ii(v,w));
			g.adj[v].push_back(ii(u,w));
		}

	}

    fin.close();
}

void read_BlatInputFile(char input_file[]) {
    int i, j;
	double w;
	ifstream fin;
    fin.open(input_file);

    while ((fin>>i>>j>>w)>0) {
		BlatSimilar[i][j] = w;
    }

    fin.close();
}

void set_parameters() {

	 if (max_tries==-1)		max_tries	= 1;			// so lan chay
     if (Nc==-1)			Nc			= 15;			// so loi giai xay dung
	 if (max_time==-1)		max_time	= 0;	   		// thoi gian chay
	 if (alpha==-1)			alpha		= 0.3;
	 if (rho==-1)			rho			= 0.8;
	 if (ls_flag==-1)		ls_flag		= 1;
	 if (kAlg==-1)		    kAlg		= 1;
	 if (n_ants==-1)		n_ants		= (kAlg == 1)? 5 : 1;


	 trail_max = 1.0;
 	 trail_min = trail_max / double (g1.n+g2.n);
	 trail_mid = (trail_min);
}

void print_parameters_log(){
	 flog<<ver_name<<endl;
	 flog<<"Graph 1                        : "<<input1_file<<endl;
	 flog<<"Graph 2                        : "<<input2_file<<endl;
	 flog<<"Output File                    : "<<output_file<<endl;
	 flog<<"(V1,E1)                        : "<<"("<<g1.n<<","<<g1.m<<")"<<endl;
	 flog<<"(V2,E2)                        : "<<"("<<g2.n<<","<<g2.m<<")"<<endl;
	 flog<<endl;
     flog<<"Parameter-settings"<<endl;
     flog<<"max-tries                      : "<<max_tries<<endl;
	 flog<<"max-time                       : "<<max_time<<"s"<<endl;
     flog<<"number of iteration            : "<<Nc<<endl;
     flog<<"num-ants                       : "<<n_ants<<endl;
	 flog<<"rho                            : "<<rho<<endl;
	 flog<<"trail_max/trail_min            : "<<trail_max/trail_min<<endl;
	 flog<<"trail_mid/trail_min            : "<<trail_mid/trail_min<<endl;

     flog<<"alpha                          : "<<alpha<<"("<<(1-alpha)<<"*blat)"<<endl;
	 flog<<"algorithm                      : "<<kAlg<<"("<<((kAlg == 1)?"+pheromone;":"-pheromone;")<<")"<<endl;
	 flog<<"local search                   : "<<ls_flag<<endl;
}

void print_parameters(){
	 cout<<ver_name<<endl;
	 cout<<"Graph 1                        : "<<input1_file<<endl;
	 cout<<"Graph 2                        : "<<input2_file<<endl;
	 cout<<"Output File                    : "<<output_file<<endl;
	 cout<<"(V1,E1)                        : "<<"("<<g1.n<<","<<g1.m<<")"<<endl;
	 cout<<"(V2,E2)                        : "<<"("<<g2.n<<","<<g2.m<<")"<<endl;
	 cout<<endl;
     cout<<"Parameter-settings"<<endl;
     cout<<"max-tries                      : "<<max_tries<<endl;
	 cout<<"max-time                       : "<<max_time<<"s"<<endl;
     cout<<"number of iteration            : "<<Nc<<endl;
     cout<<"num-ants                       : "<<n_ants<<endl;
	 cout<<"rho                            : "<<rho<<endl;
	 cout<<"trail_max/trail_min            : "<<trail_max/trail_min<<endl;
	 cout<<"trail_mid/trail_min            : "<<trail_mid/trail_min<<endl;

     cout<<"alpha                          : "<<alpha<<"("<<(1-alpha)<<"*blat)"<<endl;
	 cout<<"algorithm                      : "<<kAlg<<"("<<((kAlg == 1)?"+pheromone;":"-pheromone;")<<")"<<endl;
	 cout<<"local search                   : "<<ls_flag<<endl;

}

void preprocess_phase() {
	int i, j;

    for (i=0; i<g1.n; i++) {
        pheromone_mark[i] = new int[g2.n];
        pheromone[i] = new double[g2.n];
		BlatSimilar[i] = new double[g2.n];
    }

	//cout<<"memory -> Ok"<<endl;

	for (i=0; i<g1.n; i++)
        for (j=0; j<g2.n; j++)
			BlatSimilar[i][j] = 0.0;

	if (alpha < 1.0) {
		// su dung blat
		read_BlatInputFile(blat_input_file);

	}

    for (int i=0; i<g2.n; i++) {
         matrixGraph2[i] = new int[g2.n];
         for (int j=0; j<g2.n; j++) matrixGraph2[i][j] = 0;
    }

    for (int i=0; i<g2.n; i++)
        for (int j=0; j<g2.adj[i].size(); j++) {
            matrixGraph2[i][g2.adj[i][j].first] = 1;
            matrixGraph2[g2.adj[i][j].first][i] = 1;
        }

}

void init_program(int argc, char *argv[]){

     parse_commandline(argc, argv);

     cout<<"reading file 1: "<<input1_file;
     read_InputFile(input1_file,g1);
     cout<<" -> Ok"<<endl;

     cout<<"reading file 2: "<<input2_file;
     read_InputFile(input2_file,g2);
     cout<<" -> Ok"<<endl;

	 set_parameters();

     flog.open(log_file);

	 print_parameters();
	 print_parameters_log();

}

int termination_condition(int n_try)
/*
      FUNCTION:       checks whether termination condition is met
      INPUT:          none
      OUTPUT:         0 if condition is not met, number neq 0 otherwise
      (SIDE)EFFECTS:  none
*/
{
  return (
             (iteration * n_ants >= Nc)
           ||((max_time > 0) && ((clock() - startTry) > max_time * CLOCKS_PER_SEC))
		  );
}


void init_pheromone_trails( double init_trail ) {
	int i, j, k;

	for (i=0; i<g1.n; i++)
		for (j=0; j<g2.n; j++)
			pheromone[i][j] = init_trail;

}

void init_pheromone_mark ( int init_mark )
{
  int i, j, k;

  for (i=0; i<g1.n; i++)
		for (j=0; j<g2.n; j++)
			pheromone_mark[i][j] = init_mark;

}

void init_try(int n_try){
/*
      FUNCTION: initilialize variables appropriately when starting a trial
      INPUT:    trial number
      OUTPUT:   none
      COMMENTS: none
*/
	startTry = clock(); // bat dau tinh thoi gian cho mot lan chay

	init_pheromone_trails( trail_max );

	best_so_far_ant[n_try].tour_length = 0;

	cout<<endl<<"start try "<<n_try+1<<" ..."<<endl;
	flog<<endl<<"start try "<<n_try+1<<" ..."<<endl;
}

int next_node_to_align() {

    current_ant.maxRelation = -1;

    for (int i=0; i<g1.n; i++)
    if (current_ant.mark1[i]==0) {

        int numRelation = 0;
        for (int j=0; j<g1.adj[i].size(); j++)
            numRelation += current_ant.mark1[g1.adj[i][j].first];

        if (numRelation>current_ant.maxRelation) {
           current_ant.nn = 1;
           current_ant.nodeList[0] = i;
           current_ant.maxRelation = numRelation;
        }
        else if (numRelation==current_ant.maxRelation) {
           current_ant.nodeList[current_ant.nn] = i;
           current_ant.nn++;
        }
    }

    double total = 0;
    for (int j=0; j<current_ant.nn; j++) total += g1.adj[current_ant.nodeList[j]].size();
    double r = randomdouble() * total, r0 = 0;
    for (int j=0; j<current_ant.nn; j++) {
        r0 += g1.adj[current_ant.nodeList[j]].size();
        if (r0>=r) {
           return current_ant.nodeList[j];
        }
    }

    return (current_ant.nodeList[rand() % current_ant.nn]);
}

int find_node_to_match(int i) {

	current_ant.maxFit = -1;

    for (int j=0; j<g2.n; j++)
    if (current_ant.mark2[j]==0) {

        int numRelation = 0;
        for (int k=0; k<g1.adj[i].size(); k++) {
            int u = g1.adj[i][k].first;
            int v = current_ant.tour[u];
            if (current_ant.mark1[u] == 1) numRelation += matrixGraph2[j][v];
        }

        double tmp =  pheromone[i][j] * ( alpha * numRelation + (1-alpha) * BlatSimilar[i][j] ) ;

		if ( tmp > current_ant.maxFit ) {
		   current_ant.maxRelation = numRelation;
		   current_ant.maxFit = tmp;
		   current_ant.nn = 1;
           current_ant.nodeList[0] = j;
        }
    }

    return (current_ant.nodeList[rand() % current_ant.nn]);
}

void deep_first_search(int i, int sv) {
	current_ant.mark1[i] = sv;
	current_ant.nodes_LCCS++;
	for (int j=0; j<g1.adj[i].size(); j++) {
		int nj = g1.adj[i][j].first;
		int u = current_ant.tour[i];
		int v = current_ant.tour[nj];
		if ((current_ant.mark1[nj] == 0) && (matrixGraph2[u][v]))
			deep_first_search(nj,sv);
	}
}

void calLCCS() {
	memset(current_ant.mark1,0,sizeof(current_ant.mark1));

	int sv = 0, vung;
	int max_nodes_LCCS = 0;
	for (int i=0; i<g1.n; i++)
		if (current_ant.mark1[i] == 0) {
			sv++;
			current_ant.nodes_LCCS = 0;
			deep_first_search(i,sv);

			if (max_nodes_LCCS < current_ant.nodes_LCCS) {
				max_nodes_LCCS = current_ant.nodes_LCCS;
				vung = sv;
			}
		}
	current_ant.nodes_LCCS = max_nodes_LCCS;
	current_ant.edges_LCCS = 0;
	for (int i=0; i<g1.n; i++)
		if (current_ant.mark1[i] == vung)
			for (int j=0; j<g1.adj[i].size(); j++)
			if (current_ant.mark1[g1.adj[i][j].first] == vung) {

				int u = current_ant.tour[i];
				int v = current_ant.tour[g1.adj[i][j].first];
				current_ant.edges_LCCS += matrixGraph2[u][v];
			}
	current_ant.edges_LCCS /= 2;
}

void ant_k_construct_solutions() {
     memset(current_ant.mark1,0,sizeof(current_ant.mark1));
     memset(current_ant.mark2,0,sizeof(current_ant.mark2));

     current_ant.tour_length = 0;
     current_ant.E12 = 0;

     int i,j,k;
     for (k=0; k<g1.n; k++) {

         i = next_node_to_align();
         j = find_node_to_match(i);

         current_ant.mark1[i] = 1;
         current_ant.mark2[j] = 1;

         current_ant.tour[i] = j;

         current_ant.E12 += current_ant.maxRelation;
         current_ant.tour_length += (alpha * current_ant.maxRelation  + (1-alpha) * BlatSimilar[i][j]);

         //cout<<i<<" - "<<j<<" -> "<<current_ant.maxRelation<<endl;
     }

	 calLCCS();
}

bool local_search(ant_struct &ls_ant, int nkeep) {

    ant_struct save_ant = ls_ant;
    current_ant = ls_ant;

	double khop[MAX_NODE]={0};
	int ds[MAX_NODE];

	for (int i=0; i<g1.n; i++) {
            for (int k=0; k<g1.adj[i].size(); k++)
                khop[i] += matrixGraph2[current_ant.tour[i]][current_ant.tour[g1.adj[i][k].first]];

			khop[i] = alpha * khop[i]  + (1-alpha) * BlatSimilar[i][current_ant.tour[i]];
			ds[i] = i;

    }

	for (int x=0; x<nkeep; x++)
		for (int y=x+1; y<g1.n; y++)
			if (khop[ds[x]] < khop[ds[y]]) {
				swap(ds[x],ds[y]);
			}

    for (int x=nkeep; x<g1.n; x++)
        {
			int i = ds[x];
            int nedge = 0;

            for (int k=0; k<g1.adj[i].size(); k++)
                if (current_ant.mark1[g1.adj[i][k].first] > 0)
                    nedge += matrixGraph2[current_ant.tour[i]][current_ant.tour[g1.adj[i][k].first]];

            current_ant.tour_length -= (alpha * nedge  + (1-alpha) * BlatSimilar[i][current_ant.tour[i]]);
            current_ant.E12 -= nedge;

            current_ant.mark1[i] = 0;
            current_ant.mark2[ls_ant.tour[i]] = 0;

        }

	int nRebuild = g1.n - nkeep;
    for (int k=0; k<nRebuild; k++) {

         int i = next_node_to_align();
         int j = find_node_to_match(i);

         current_ant.mark1[i] = 1;
         current_ant.mark2[j] = 1;

         current_ant.tour[i] = j;

         current_ant.E12 += current_ant.maxRelation;
         current_ant.tour_length += (alpha * current_ant.maxRelation + (1-alpha) * BlatSimilar[i][j]);

     }

    if (current_ant.tour_length > save_ant.tour_length) {
        ls_ant = current_ant;
		calLCCS();
        return true;
    }
    else {
        ls_ant = save_ant;
        return false;
    }
}

void construct_solutions(int n_try){

	int ak;

	ibest_ant.tour_length = 0;
	init_pheromone_mark(0);

	// tung kien xay dung loi giai
	for (ak=0; ak<n_ants; ak++) {

		ant_k_construct_solutions();

        cout<<iteration * n_ants + ak	<<" : "	<<current_ant.tour_length<<"("
			<<current_ant.E12<<","
			<<current_ant.nodes_LCCS<<","
			<<current_ant.edges_LCCS<<")";

		flog<<iteration * n_ants + ak	<<" : "	<<current_ant.tour_length<<"("
            <<current_ant.E12<<","
			<<current_ant.nodes_LCCS<<","
			<<current_ant.edges_LCCS<<")";

		if (ls_flag == 1) {
			while (local_search(current_ant,int(g1.n * 0.01))) {
				// local search for all ant
				cout<<" -> "<<current_ant.tour_length<<"("
				<<current_ant.E12<<","
				<<current_ant.nodes_LCCS<<","
				<<current_ant.edges_LCCS<<")";

				flog<<" -> "<<current_ant.tour_length<<"("
				<<current_ant.E12<<","
				<<current_ant.nodes_LCCS<<","
				<<current_ant.edges_LCCS<<")";
			}
		}

		cout<<endl;
		flog<<endl;

		// cap nhat loi giai tot nhat, vong lap, thoi gian
		if (current_ant.tour_length > ibest_ant.tour_length) {
			ibest_ant = current_ant;
		}
	}

}

void TLAS_pheromone_trail_update(int n_try){

	int i, j;

	// cap nhat cac canh thuoc ibest_ant
	// cap nhat theo -> trail_max
	for (i=0; i<g1.n; i++) {
        j = ibest_ant.tour[i];
        //j = best_so_far_ant[n_try].tour[i];

		pheromone[i][j] = rho*pheromone[i][j] + (1-rho)*trail_max;
		pheromone_mark[i][j] = 2;
	}

	// cap nhat cac canh KHONG thuoc ibest_ant
	// cap nhat theo -> trail_min
	for (i=0; i<g1.n; i++)
		for (j=0; j<g2.n; j++) {
				if (pheromone_mark[i][j]==1) {

					pheromone[i][j] = rho*pheromone[i][j] + (1-rho)*trail_mid;

				}
				else if (pheromone_mark[i][j]==0){
					pheromone[i][j] = rho*pheromone[i][j] + (1-rho)*trail_min;
				}
        }
}

void search_control_and_statistics(int n_try){
/*
      FUNCTION:       occasionally compute some statistics and check whether algorithm
                      is converged
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated
*/
	// cap nhat loi giai tot nhat trong lan chay n_try
	if (ibest_ant.tour_length > best_so_far_ant[n_try].tour_length) {
		best_so_far_ant[n_try] = ibest_ant;
		best_so_far_ant[n_try].iteration = iteration;
		best_so_far_ant[n_try].time = clock() - startTry;
		cout<<"Iteration "<<setw(5)<<iteration
						  <<" time "<<setw(7)<<setprecision(3)<<fixed<<(double)(best_so_far_ant[n_try].time)/CLOCKS_PER_SEC
						  <<" -----> Best : "<<best_so_far_ant[n_try].tour_length
						  <<"("
						  <<best_so_far_ant[n_try].E12<<","
						  <<best_so_far_ant[n_try].nodes_LCCS<<","
						  <<best_so_far_ant[n_try].edges_LCCS<<")"
						  <<endl;


		flog<<"Iteration "<<setw(5)<<iteration
						  <<" time "<<setw(7)<<setprecision(3)<<fixed<<(double)(best_so_far_ant[n_try].time)/CLOCKS_PER_SEC
						  <<" -----> Best : "<<best_so_far_ant[n_try].tour_length
						  <<"("
						  <<best_so_far_ant[n_try].E12<<","
						  <<best_so_far_ant[n_try].nodes_LCCS<<","
						  <<best_so_far_ant[n_try].edges_LCCS<<")"
						  <<endl;
	}

}

void exit_try(int n_try){
	cout<<"Best Solution in try "<<setw(3)<<n_try+1<<" is "<<best_so_far_ant[n_try].tour_length
										<<"("
										<<best_so_far_ant[n_try].E12<<","
										<<best_so_far_ant[n_try].nodes_LCCS<<","
										<<best_so_far_ant[n_try].edges_LCCS<<")"
			                            <<" at iterations "<<best_so_far_ant[n_try].iteration
										<<" time "<<best_so_far_ant[n_try].time / (double)CLOCKS_PER_SEC<<endl;

	flog<<"Best Solution in try "<<setw(3)<<n_try+1<<" is "<<best_so_far_ant[n_try].tour_length
										<<"("
										<<best_so_far_ant[n_try].E12<<","
										<<best_so_far_ant[n_try].nodes_LCCS<<","
										<<best_so_far_ant[n_try].edges_LCCS<<")"
			                            <<" at iterations "<<best_so_far_ant[n_try].iteration
										<<" time "<<best_so_far_ant[n_try].time / (double)CLOCKS_PER_SEC<<endl;

}

void exit_program(){
	int id_best;
	double best_tour_length, worst_tour_length, avg_tour_length;

	int i, j;

	best_tour_length  = best_so_far_ant[0].tour_length; id_best = 0;
	worst_tour_length = best_so_far_ant[0].tour_length;
	avg_tour_length   = best_so_far_ant[0].tour_length;

	for (i=1; i<max_tries; i++) {
		if (best_so_far_ant[i].tour_length > best_tour_length) {
			best_tour_length = best_so_far_ant[i].tour_length;
			id_best = i;
		}

		if (best_so_far_ant[i].tour_length < worst_tour_length) {
			worst_tour_length = best_so_far_ant[i].tour_length;
		}

		avg_tour_length += best_so_far_ant[i].tour_length;
	}

	avg_tour_length = avg_tour_length / (double) max_tries;

	cout<<"Best : "<<best_tour_length<<"("<<best_so_far_ant[id_best].E12<<","<<best_so_far_ant[id_best].nodes_LCCS<<","<<best_so_far_ant[id_best].edges_LCCS<<")"
	<<" Worst : "<<worst_tour_length<<" Avg : "<< avg_tour_length<<endl;
	flog<<"Best : "<<best_tour_length<<"("<<best_so_far_ant[id_best].E12<<","<<best_so_far_ant[id_best].nodes_LCCS<<","<<best_so_far_ant[id_best].edges_LCCS<<")"
	<<" Worst : "<<worst_tour_length<<" Avg : "<< avg_tour_length<<endl;

	// ghi vao file ket qua cac thong tin
	// 1. cac tham so
	// 2. ket qua tung lan chay
	// 3. thong ke
	// 4. loi giai tot nhat


	ofstream fout;

    fout.open(output_file);
	fout<<best_so_far_ant[id_best].E12<<endl;
	for (int i=0; i<g1.n; i++)
	    fout<<i<<" "<<best_so_far_ant[id_best].tour[i]<<endl;
	fout.close();


    fout.open(output2_file);
	fout<<best_so_far_ant[id_best].E12<<endl;
	for (int i=0; i<g1.n; i++)
	    fout<<best_so_far_ant[id_best].tour[i]<<" "<<i<<endl;
	fout.close();

	endProgram = clock();

	cout<<"time total : "<<(double)(endProgram - startProgram)/CLOCKS_PER_SEC<<"s"<<endl;
	flog<<"time total : "<<(double)(endProgram - startProgram)/CLOCKS_PER_SEC<<"s"<<endl;
	flog.close();

}


int main(int argc, char *argv[])
{
    startProgram = clock(); // bat dau tinh thoi gian chay chuong trinh

    init_program(argc, argv);

    if (g1.n > g2.n) {
       cout<<"ERROR: graph1 > graph2"<<endl;
       return 1;
    }

    cout<<"preprocess_phase";
    preprocess_phase();

    cout<<" -> Ok at time "<<(double)(clock() - startProgram)/CLOCKS_PER_SEC<<"s"<<endl;

	srand( seed );

	for ( n_try = 0 ; n_try < max_tries ; n_try++ ) {

		init_try(n_try);

		iteration=0;

		while ( !termination_condition(n_try) ) {

			construct_solutions(n_try);

			search_control_and_statistics(n_try);

			if (kAlg == 1) TLAS_pheromone_trail_update(n_try);

			iteration++;
		}
		exit_try(n_try);
    }
    exit_program();

    //system("PAUSE");
    return EXIT_SUCCESS;
}

