// author	: Do Duc Dong, ITI.VNU
//            Dang Cao Cuong, UET.VNU
// date		: 31/1/2014


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

const int MAX_NODE    = 20000;

char input1_file[100]		= "g1.ga";
char input2_file[100]		= "g2.ga";
char output_file[100]		= "g1.ga_g2.ga.out";
char ver_name[100]			= "ACO for Graph Alignment Problem";
char log_file[100]			= "g1.ga_g2.ga.log";
ofstream flog;

clock_t startProgram,			// do thoi gian chay chuong trinh
		endProgram,
		startTry,				// do thoi gian mot lan chay
		endTry;

int max_tries,                  // so lan chay
	n_try,
    Nc,                         // so vong lap
	max_reset,					// so lan reset vet mui
	iteration,
	seed,
    reset_count,                // so lan reset lai vet mui
    sl;                         // kiem soat hoi tu

int max_time;					// thoi gian chay toi da

int n_ants;                     // so kien

double alpha,
       beta,
       rho;

int	ls_flag;

int kAlg; 
// 0. -graphlet, -pheromone, -blat
// 1. +graphlet, -pheromone, -blat
// 2. -graphlet, +pheromone, -blat
// 3. +graphlet, +pheromone, -blat

struct ant_struct{
  int       tour[MAX_NODE];
  int       mark1[MAX_NODE], mark2[MAX_NODE];
  int		tour_length; 
  int 		nodes_LCCS, edges_LCCS;
  int 	    nn, nodeList[MAX_NODE];  
  double	fit[MAX_NODE];
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
  int         orbit[MAX_NODE][73];
  double      orbitLog1[MAX_NODE][73], orbitLog2[MAX_NODE][73];
  int		  n;
  int		  m;
};

graph g1, g2;
double *similar[MAX_NODE];        // thong tin heuristic dua tren graphlet

double randomdouble();
void gen_permutation(int p[], int psize);

void parse_commandline(int argc, char *argv[]) {

	 max_tries				= -1;
     Nc						= -1;
	 max_time				= -1;
	 max_reset				= -1;
     n_ants					= -1;
	 alpha					= -1;
	 beta					= -1;
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

		 // -s : so lan reset pheremone
         if (strcmp(argv[i],"-s")==0) {
              max_reset = atoi(argv[i+1]);
         }		 
         
         // -k : loai thuat toan
         if (strcmp(argv[i],"-k")==0) {
              kAlg = atoi(argv[i+1]);
         }

         // -seed : so ngau nhien ban dau
         if (strcmp(argv[i],"-seed")==0) {
              seed = atoi(argv[i+1]);
         }
		 else seed = (unsigned)time( NULL );

         // -c : so vong lap
         if (strcmp(argv[i],"-c")==0) {
			  max_reset = 0;
              Nc = atoi(argv[i+1]);
         }

		 // -t : thoi gian chay
         if (strcmp(argv[i],"-t")==0) {
			  max_reset = 0;
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

         // -b : beta
         if (strcmp(argv[i],"-b")==0) {
              beta = atof(argv[i+1]);
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
}

void read_InputFile(char input_file[], graph &g) {
    int i, j;
	ifstream fin;
    fin.open(input_file);

    fin>>g.n>>g.m;

    g.adj.assign(g.n + 1 , vii());

    int u, v, w;

    for (int i=0; i<g.m; i++) {
        fin>>u>>v; w = 1;

        g.adj[u].push_back(ii(v,w));
        g.adj[v].push_back(ii(u,w));
    }

    for (int i=0; i<g.n; i++)
        for (j=0; j<73; j++) {
            fin>>g.orbit[i][j];
            g.orbitLog1[i][j] = log(g.orbit[i][j]+1);
            g.orbitLog2[i][j] = log(g.orbit[i][j]+2);
        }

    fin.close();
}

void set_parameters() {

	 if (max_tries==-1)		max_tries	= 10;			// so lan chay
     if (Nc==-1)			Nc			= 300;			// so loi giai xay dung
	 if (max_time==-1)		max_time	= 0;	   		// thoi gian chay
	 if (max_reset==-1)		max_reset	= 10;	   		// so lan khoi tao lai vet mui
     if (n_ants==-1)		n_ants		= 30;
	 if (alpha==-1)			alpha		= 1.0;
	 if (beta==-1)			beta		= 1.0;
	 if (rho==-1)			rho			= 0.5;
	 if (ls_flag==-1)		ls_flag		= 4;
	 if (kAlg==-1)		    kAlg		= 1;


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
	 flog<<"number of reset pheremone      : "<<max_reset<<endl;
     flog<<"num-ants                       : "<<n_ants<<endl;
     flog<<"alpha                          : "<<alpha<<endl;
     flog<<"beta                           : "<<beta<<endl;
     flog<<"rho                            : "<<rho<<endl;
	 flog<<"local search                   : "<<ls_flag<<endl;
	 flog<<"algorithm                      : "<<kAlg<<endl;
	 flog<<"trail_max/trail_min            : "<<trail_max/trail_min<<endl;
	 flog<<"trail_mid/trail_min            : "<<trail_mid/trail_min<<endl;
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
	 cout<<"number of reset pheremone      : "<<max_reset<<endl;
     cout<<"num-ants                       : "<<n_ants<<endl;
     cout<<"alpha                          : "<<alpha<<endl;
     cout<<"beta                           : "<<beta<<endl;
     cout<<"rho                            : "<<rho<<endl;
	 cout<<"local search                   : "<<ls_flag<<endl;
	 cout<<"algorithm                      : "<<kAlg<<endl;
	 cout<<"trail_max/trail_min            : "<<trail_max/trail_min<<endl;
	 cout<<"trail_mid/trail_min            : "<<trail_mid/trail_min<<endl;
}

double calSimilar(int i, int j) {
       
    // khong su dung graphlet
    if ((kAlg & 1) == 0) return 1;
    
	int k;
	double t=0;
	for (k=0; k<73; k++) {
		if (g1.orbit[i][k] > g2.orbit[j][k]) 
           t += (g1.orbitLog1[i][k] - g2.orbitLog1[j][k])/g1.orbitLog2[i][k];
        else 
           t += (g2.orbitLog1[j][k] - g1.orbitLog1[i][k])/g2.orbitLog2[j][k];  
    }

    t = t / 73;
	return (1-t);
}

void preprocess_phase() {
	int i, j;

    for (i=0; i<g1.n; i++) {
        pheromone_mark[i] = new int[g2.n];
        pheromone[i] = new double[g2.n];
        similar[i] = new double[g2.n];
    }
			
	for (i=0; i<g1.n; i++) {
        for (j=0; j<g2.n; j++) {
			similar[i][j] = calSimilar(i,j);
        }
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
           &&((clock() - startTry) > max_time * CLOCKS_PER_SEC)
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
	reset_count = 0;

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
    
    return (current_ant.nodeList[rand() % current_ant.nn]);
}

int find_node_to_match(int i) {
    
    current_ant.maxRelation = -1;
                                 
    for (int j=0; j<g2.n; j++) 
    if (current_ant.mark2[j]==0) {
                                 
        int numRelation = 0;
        for (int k=0; k<g1.adj[i].size(); k++) {
            int u = g1.adj[i][k].first;
            int v = current_ant.tour[u];
            if (current_ant.mark1[u] == 1) numRelation += matrixGraph2[j][v];
        }
		
		if (current_ant.maxRelation < numRelation) {
			current_ant.maxRelation = numRelation;
			current_ant.fit[0] =  pheromone[i][j] * similar[i][j];
			current_ant.nodeList[0] = j;
			current_ant.nn = 1;
		}
		else if (current_ant.maxRelation == numRelation) {
			current_ant.fit[current_ant.nn] =  pheromone[i][j] * similar[i][j];
			current_ant.nodeList[current_ant.nn] = j;
			current_ant.nn ++;
		}
    }
    
    int nn = 2;
    for (int x=0; x<min(current_ant.nn,nn); x++)
        for (int y=x+1; y<current_ant.nn; y++)
            if (current_ant.fit[x]<current_ant.fit[y]) {
               swap(current_ant.fit[x],current_ant.fit[y]);
               swap(current_ant.nodeList[x],current_ant.nodeList[y]);
            }
            
    // voi xac suat q0 thi chon gia tri lon nhat
    double q0=0.99;
    if (randomdouble()<=q0) 
        return current_ant.nodeList[0];    
        
    current_ant.nn = min(current_ant.nn,nn);
    
    double total = 0;
    for (int j=0; j<current_ant.nn; j++) total += current_ant.fit[j];
    double r = randomdouble() * total, r0 = 0;
    for (int j=0; j<current_ant.nn; j++) {
        r0 += current_ant.fit[j];
        if (r0>=r) {
           return current_ant.nodeList[j];
        }
    }
    
    return (current_ant.nodeList[current_ant.nn-1]);
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
     
     int i,j,k;
     for (k=0; k<g1.n; k++) {
         
         i = next_node_to_align();
         j = find_node_to_match(i);
         
         current_ant.mark1[i] = 1;
         current_ant.mark2[j] = 1;
         
         current_ant.tour[i] = j;
         current_ant.tour_length += current_ant.maxRelation;
         
         //cout<<i<<" - "<<j<<" -> "<<current_ant.maxRelation<<endl;
     }
     
	 calLCCS();	
}

bool local_search(ant_struct &ls_ant) {
    
    ant_struct save_ant = ls_ant;
    current_ant = ls_ant;
    
	int khop[MAX_NODE]={0}, ds[MAX_NODE];
	
	for (int i=0; i<g1.n; i++) {
            for (int k=0; k<g1.adj[i].size(); k++) 
                khop[i] += matrixGraph2[current_ant.tour[i]][current_ant.tour[g1.adj[i][k].first]];
			ds[i] = i;
    }
	
	int nfit = int(g1.n * 0.01);
	for (int x=0; x<nfit; x++)
		for (int y=x+1; y<g1.n; y++)
			if (khop[ds[x]] < khop[ds[y]]) {
				swap(ds[x],ds[y]);
			}
		
    for (int x=nfit; x<g1.n; x++) 
        {
			int i = ds[x];
            
            for (int k=0; k<g1.adj[i].size(); k++) 
                if (current_ant.mark1[g1.adj[i][k].first] > 0)
                   current_ant.tour_length -= matrixGraph2[current_ant.tour[i]][current_ant.tour[g1.adj[i][k].first]];
            
            current_ant.mark1[i] = 0;
            current_ant.mark2[ls_ant.tour[i]] = 0;
           
        }
    
	int nRebuild = g1.n - nfit;
    for (int k=0; k<nRebuild; k++) {
         
         int i = next_node_to_align();
         int j = find_node_to_match(i);
         
         current_ant.mark1[i] = 1;
         current_ant.mark2[j] = 1;
         
         current_ant.tour[i] = j;
         current_ant.tour_length += current_ant.maxRelation;
    
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

        cout<<ak<<" : "	<<current_ant.tour_length<<"("
						<<current_ant.nodes_LCCS<<","    
						<<current_ant.edges_LCCS<<") -> ";      
		if (ls_flag==4) {
			while (local_search(current_ant)) {
				// local search for all ant
			}
		}

		if (ls_flag==3) {
			while (local_search(current_ant)) {
				// local search for all ant
			}
		}
		
		cout<<current_ant.tour_length<<"("
			<<current_ant.nodes_LCCS<<","    
			<<current_ant.edges_LCCS<<")"<<endl;

		// cap nhat loi giai tot nhat, vong lap, thoi gian
		if (current_ant.tour_length > ibest_ant.tour_length) {
			ibest_ant = current_ant;
		}
	}

	cout<<"Iteration "<<iteration<<" : "<<ibest_ant.tour_length<<"("
										<<ibest_ant.nodes_LCCS<<","    
										<<ibest_ant.edges_LCCS<<") -> ";
	if (ls_flag==2) {
		while (local_search(ibest_ant)) {
			// local search for ibest
		}
	}

	if (ls_flag==1) {
		local_search(ibest_ant);
	}
	cout<<ibest_ant.tour_length<<"("
		<<ibest_ant.nodes_LCCS<<","    
		<<ibest_ant.edges_LCCS<<")"<<endl;
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
						  <<" -> Best : "<<best_so_far_ant[n_try].tour_length
						  <<"("<<best_so_far_ant[n_try].nodes_LCCS<<","    
						  <<best_so_far_ant[n_try].edges_LCCS<<")"
						  <<endl;


		flog<<"Iteration "<<setw(5)<<iteration
						  <<" time "<<setw(7)<<setprecision(3)<<fixed<<(double)(best_so_far_ant[n_try].time)/CLOCKS_PER_SEC
						  <<" -> Best : "<<best_so_far_ant[n_try].tour_length
						  <<"("<<best_so_far_ant[n_try].nodes_LCCS<<","    
						  <<best_so_far_ant[n_try].edges_LCCS<<")"
						  <<endl;
		sl=1;
	}
	else if (ibest_ant.tour_length == best_so_far_ant[n_try].tour_length) {
		sl++;
	}

	if (sl==30) {
		cout<<"Init trail at "<<iteration<<" time "<<(double)(clock() - startTry)/CLOCKS_PER_SEC<<endl;
		flog<<"Init trail at "<<iteration<<" time "<<(double)(clock() - startTry)/CLOCKS_PER_SEC<<endl;
		init_pheromone_trails(trail_max);
		init_pheromone_mark(0);
		TLAS_pheromone_trail_update(n_try);
		sl=0;
		reset_count++;
	}

}

void exit_try(int n_try){
	cout<<"Best Solution in try "<<setw(3)<<n_try+1<<" is "<<best_so_far_ant[n_try].tour_length
			                            <<" at iterations "<<best_so_far_ant[n_try].iteration
										<<" time "<<best_so_far_ant[n_try].time / (double)CLOCKS_PER_SEC<<endl;

	flog<<"Best Solution in try "<<setw(3)<<n_try+1<<" is "<<best_so_far_ant[n_try].tour_length
			                            <<" at iterations "<<best_so_far_ant[n_try].iteration
										<<" time "<<best_so_far_ant[n_try].time / (double)CLOCKS_PER_SEC<<endl;

}

void exit_program(){
	int best_tour_length, id_best, worst_tour_length;
	double avg_tour_length;

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

	cout<<"Best : "<<best_tour_length<<" Worst : "<<worst_tour_length<<" Avg : "<< avg_tour_length<<endl;
	flog<<"Best : "<<best_tour_length<<" Worst : "<<worst_tour_length<<" Avg : "<< avg_tour_length<<endl;

	// ghi vao file ket qua cac thong tin
	// 1. cac tham so
	// 2. ket qua tung lan chay
	// 3. thong ke
	// 4. loi giai tot nhat


	ofstream fout;
    fout.open(output_file);

	fout<<best_tour_length<<endl;
	
	for (int i=0; i<g1.n; i++)
	    fout<<i<<" "<<best_so_far_ant[id_best].tour[i]<<endl;
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

			if ((kAlg & 2) > 0) TLAS_pheromone_trail_update(n_try);

			iteration++;
		}
		exit_try(n_try);
    }
    exit_program();

    //system("PAUSE");
    return EXIT_SUCCESS;
}

double randomdouble(){
	int x=rand();
	double y=(double)x /(double) RAND_MAX;
	return y;
}

void gen_permutation(int p[], int psize) {
	int i,j, itmp;

	for (i=0; i<psize; i++) p[i]=i;

	for (i=0; i<psize; i++) {
		j=rand() % psize;
		itmp = p[i];
		p[i] = p[j];
		p[j] = itmp;
	}
}
