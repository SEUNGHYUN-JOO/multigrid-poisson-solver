#include <stdio.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <iostream>

using namespace std;


/**************변수 정의**************/
#define NUM_NODES 17        // Number of nodes in the i & j directions
#define MG_CYCLES 1         // Maximum number of multigrid cycles
#define MG_LEVELS 4         // Number of multigrid levels
#define NUM_SWEEP 3         // Number of smoothing sweeps at each stage of the multigrid
#define VISUALIZE 1         // Flag controlling whether to write Tecplot mesh/solution files (Yes=1,No=0)
#define FREQUENCY 1         // Iteration frequency with which to print to console and write output files
#define TOLERANCE 12.0      // Convergence criteria in terms of orders reduction in the L2 norm
#define pow_tol -10.0       // tolerance



/*finest grid 상수*/
const int FINE_MESH = 0;


/**************함수 정의**************/
int allocate_arrays(double ****num_sol, double ***exact_sol, double ****f, double ****x, double ****y, double ****mid_sol, int n_nodes);            //array 동적 메모리 할당
void generate_finest_grid(double **x, double **y, int n_nodes, bool visualize);                                                               //fine grid 생성
void generate_coarse_gird(double ***x, double ***y, int n_nodes, int n_level, int level, bool visualize);                                          //coarse grid 생성
void intitialize_solution(double **num_sol, double **exact_sol, double **f, double **x, double **y, int n_nodes);                               //fine grid에서 초기해, exact_sol 설정
void multigrid_cycle(double ***num_sol, double ***f, double ***mid_sol, int n_nodes, int n_smooth, int n_level, int level);                        //MG 4-level V-cycle 수행
void gauss_seidel(double **num_sol, double **f, double **mid_sol, int n_nodes, int n_smooth, double tolerance);                                                //Gauss-Seidel로 smoothing
void restriction_step(double ***num_sol, double ***f, double ***mid_sol, int n_nodes, int level);                                                   //restriction step - fine -> coarse grid
void prolongate_step(double ***num_sol, double ***mid_sol, int n_nodes, int level);                                                                 //prolongate step - coarse -> fine grid
double compute_residual(double **num_sol, double **f, double **residual, int n_nodes);                                                          //residual 계산
void write_output(double **num_sol, double **exact_sol, double **x, double **y, int n_nodes, int iter, double residual, bool visualize);        //result .dat 파일로 저장
void deallocate_arrays(double ***num_sol, double **exact_sol, double ***f, double ***x, double ***y, double ***mid_sol, int n_nodes, int n_level); //동적 메모리 제거


int allocate_arrays(double ****num_sol, double ***exact_sol, double ****f, double ****x, double ****y, double ****mid_sol, int n_nodes){  
    int n_level = MG_LEVELS;
    int i_nodes = n_nodes;    
    int i_max;
    int j_max;
    double ***num_sol_temp;
    double ***mid_sol_temp;
    double ***f_temp;
    double ***x_temp; 
    double ***y_temp;
    double **exact_sol_temp;      
        
    //n_level 만큼 array 생성
    x_temp       = new double**[n_level];
    y_temp       = new double**[n_level];
    num_sol_temp = new double**[n_level];
    mid_sol_temp = new double**[n_level];
    f_temp       = new double**[n_level];
  
    for (int idx_level = 0; idx_level < n_level; idx_level++){
        //각 level에서 node수만큼  x, y, num_sol, mid_sol, f array 생성
        x_temp[idx_level]       = new double*[i_nodes];
        y_temp[idx_level]       = new double*[i_nodes];
        num_sol_temp[idx_level] = new double*[i_nodes];
        mid_sol_temp[idx_level] = new double*[i_nodes];
        f_temp[idx_level]       = new double*[i_nodes];
            
        if (idx_level == 0){
            exact_sol_temp = new double*[i_nodes];
        }
        
        for (int i = 0; i < i_nodes; i++){
            x_temp[idx_level][i]       = new double[i_nodes];
            y_temp[idx_level][i]       = new double[i_nodes];
            num_sol_temp[idx_level][i] = new double[i_nodes];
            mid_sol_temp[idx_level][i] = new double[i_nodes];
            f_temp[idx_level][i]       = new double[i_nodes];
                  
            if (idx_level == 0) exact_sol_temp[i] = new double[i_nodes];
        }
            
        i_nodes = (i_nodes-1)/2 + 1;    
    }
  
    //동적 메모리 함수 외부에서 사용
    *x         = x_temp;
    *y         = y_temp;
    *num_sol       = num_sol_temp;
    *mid_sol       = mid_sol_temp;
    *f         = f_temp;
    *exact_sol = exact_sol_temp;
      
    return 0;
}


void generate_finest_grid(double **x, double **y, int n_nodes, bool visualize){
    // 0<= x,y <=1의 n_nodes로 가장 fine한 grid 생성
    for (int i = 0; i < n_nodes; i++){
        for (int j = 0; j < n_nodes; j++){
            x[i][j] = (double) i/ (double) (n_nodes-1);
            y[i][j] = (double) j/ (double) (n_nodes-1);
        }
    }
}

void generate_coarse_gird(double ***x, double ***y, int n_nodes, int n_level, int level, bool visualize){
    if (level != n_level-1){
        for (int i = 0; i < n_nodes; i++){
            for (int j = 0; j < n_nodes; j++){        
                //coarse grid는 전 level의 i, j 값을 i/2, j/2에 가짐
                if ((i%2 == 0) && (j%2 == 0)){
                    x[level+1][i/2][j/2] = x[level][i][j];
                    y[level+1][i/2][j/2] = y[level][i][j];
                }        
            }
        }
        //가장 coarse한 gird될 때까지 반복(n_level)
        int n_coarse_nodes = (n_nodes-1)/(2) + 1;
        int levels_coarse = level+1;    
        generate_coarse_gird(x, y, n_coarse_nodes, n_level, levels_coarse, visualize);                
    }  
}

void intitialize_solution(double **num_sol, double **exact_sol, double **f, double **x, double **y, int n_nodes){
    //가장 fine한 grid에서 정의하면 coarse grid는 할 필요 없음
    //문제 2번 BC가 Drichlet이므로 boundary만 값 주고 나머지는 0으로 초기화(문제2는 boundary도 0.0)
    float x_temp;
    float y_temp;
    
    for (int i = 0; i < n_nodes; i++){
        for (int j = 0; j < n_nodes; j++){
            x_temp = x[i][j];
            y_temp = y[i][j];
            
            num_sol[i][j]   = 0.0;
            exact_sol[i][j] = -1 * x_temp*x_temp*y_temp*y_temp*(1-x_temp*x_temp)*(1-y_temp*y_temp);
            f[i][j]         = -2 * ((1-6*x_temp*x_temp)*y_temp*y_temp*(1-y_temp*y_temp) + (1-6*y_temp*y_temp)*x_temp*x_temp*(1-x_temp*x_temp));
        }
    }  
}

void multigrid_cycle(double ***num_sol, double ***f, double ***mid_sol, int n_nodes, int n_smooth, int n_level, int level, double tolerance){
    //Gauss-Seidel로 n_smooth번 smoothing
    gauss_seidel(num_sol[level], f[level], mid_sol[level], n_nodes, n_smooth, tolerance);    
    if (level < n_level-1){
        restriction_step(num_sol, f, mid_sol, n_nodes, level);    
        
        int n_coarse     = (n_nodes-1)/(2) + 1;
        int level_coarse = level + 1;      
        multigrid_cycle(num_sol, f, mid_sol, n_coarse, n_smooth, n_level, level_coarse);
     
        prolongate_step(num_sol, mid_sol, n_nodes, level);    
    }  
    gauss_seidel(num_sol[level], f[level], mid_sol[level], n_nodes, n_smooth, tolerance);  
}

void gauss_seidel(double **num_sol, double **f, double **mid_sol, int n_nodes, int n_smooth, double tolerance){
    double dx = 1.0/((double)n_nodes-1.0);
    double dy = 1.0/((double)n_nodes-1.0);
    double beta = dy/dx;
    double coeff = 1.0/(2.0*(1+beta*beta));
    double sigma = (cos(M_PI*dx) + beta*beta*cos(M_PI*dy))/(1+beta*beta); //for Gauss-Seidel over relaxation    
    double weight = 2.0 / (1.0 + sqrt(1-sigma*sigma)); //for Gauss-Seidel over relaxation    
    double num_sol_old, error;
        
    for (int iter = 0; iter < n_smooth; iter++){
        error = 0.0;
        for (int i = 1; i < n_nodes-1; i++){            
            for (int j = 1; j < n_nodes-1; j++){
                num_sol_old = num_sol[i][j];                
                num_sol[i][j] = coeff*(num_sol[i+1][j] + num_sol[i-1][j] + beta*beta*(num_sol[i][j+1] + num_sol[i][j-1]));
                error += (num_sol[i][j] - num_sol_old) * (num_sol[i][j] - num_sol_old);
            }
        }

        if (error < tolerance){
            break;
        }
        else{


        }

    }
}

void restriction_step(double ***num_sol, double ***f, double ***mid_sol, int n_nodes, int level){  
    
    compute_residual(num_sol[level], f[level], mid_sol[level], n_nodes);  
    
    int n_coarse     = (n_nodes-1)/(2) + 1;
    int level_coarse = level+1;
    int i_fine, j_fine;
  
  //! Initialize the solution guess and forcing term for the coarse level
  
  for (int i = 0; i < n_coarse; i++) {
    for (int j = 0; j < n_coarse; j++) {
      num_sol[level_coarse][i][j] = 0.0;
      f[level_coarse][i][j]   = 0.0;
    }
  }
  
  //! Transfer the forcing term to the coarse mesh. We are
  //! restricting by weighting the values on the fine mesh that
  //! surround the given coarse node. The residual on the boundary
  //! nodes should be zero, so we avoid those points in our loop.
  
  for (int i = 1; i < n_coarse-1; i++) {
    for (int j = 1; j < n_coarse-1; j++) {
      
      //! Calculate the indices on the fine mesh for clarity
      
      i_fine = (i*2); j_fine = (j*2);
      
      //! Perform the restriction operation for this node by injection
      
      f[level_coarse][i][j] = (  mid_sol[level][i_fine-1][j_fine+1]*(1.0/16.0)
                               + mid_sol[level][ i_fine ][j_fine+1]*(1.0/8.0)
                               + mid_sol[level][i_fine+1][j_fine+1]*(1.0/16.0)
                               + mid_sol[level][i_fine-1][ j_fine ]*(1.0/8.0)
                               + mid_sol[level][ i_fine ][ j_fine ]*(1.0/4.0)
                               + mid_sol[level][i_fine+1][ j_fine ]*(1.0/8.0)
                               + mid_sol[level][i_fine-1][j_fine-1]*(1.0/16.0)
                               + mid_sol[level][ i_fine ][j_fine-1]*(1.0/8.0)
                               + mid_sol[level][i_fine+1][j_fine-1]*(1.0/16.0));
      
    }
  }
  
}

void prolongate_step(double ***num_sol, double ***mid_sol, int n_nodes,
                         int level) {
  
  //! Prolongate the solution, i.e., transfer it from a coarse to fine level.
  //! We compute a correction to the fine grid solution and add. To create
  //! the correction from the coarse solution, nodes that lie on top of each
  //! other will be transfered directly from coarse->fine, while neighbors
  //! are interpolated using a weighting of neighbors.
  
  //! Compute some information about the coarse level
  
  int n_coarse     = (n_nodes-1)/(2) + 1;
  int level_coarse = level+1;
  int i_fine, j_fine;
  
  //! Initialize correction to zero, just in case
  
  for (int i = 0; i < n_nodes; i++) {
    for (int j = 0; j < n_nodes; j++) {
      mid_sol[level][i][j] = 0.0;
    }
  }
  
  //! Loop over all coarse points and set the solution on the
  //! coincident set of points on the fine grid. At the same time,
  //! contribute weighted values at coarse grid nodes to the set of
  //! neighboring fine grid nodes that are not coincident.
  
  for (int i = 1; i < n_coarse-1; i++) {
    for (int j = 1; j < n_coarse-1; j++) {
      
      //! Calculate the indices on the fine mesh for clarity
      
      i_fine = i*2; j_fine = j*2;
      
      //! Perform the prolongation operation by copying the value for
      //! a coincident node on the fine mesh and also incrementing the
      //! values for the neighbors.
      
      mid_sol[level][i_fine-1][j_fine+1] += num_sol[level_coarse][i][j]*(1.0/4.0);
      mid_sol[level][ i_fine ][j_fine+1] += num_sol[level_coarse][i][j]*(1.0/2.0);
      mid_sol[level][i_fine+1][j_fine+1] += num_sol[level_coarse][i][j]*(1.0/4.0);
      mid_sol[level][i_fine-1][ j_fine ] += num_sol[level_coarse][i][j]*(1.0/2.0);
      mid_sol[level][ i_fine ][ j_fine ]  = num_sol[level_coarse][i][j];
      mid_sol[level][i_fine+1][ j_fine ] += num_sol[level_coarse][i][j]*(1.0/2.0);
      mid_sol[level][i_fine-1][j_fine-1] += num_sol[level_coarse][i][j]*(1.0/4.0);
      mid_sol[level][ i_fine ][j_fine-1] += num_sol[level_coarse][i][j]*(1.0/2.0);
      mid_sol[level][i_fine+1][j_fine-1] += num_sol[level_coarse][i][j]*(1.0/4.0);
      
    }
  }
  
  //! Finally, add the coarse grid correction to the fine grid,
  //! i.e., num_sol_fine = num_sol_fine + prolong(num_sol_coarse)
  
  for (int i = 0; i < n_nodes; i++) {
    for (int j = 0; j < n_nodes; j++) {
      num_sol[level][i][j] += mid_sol[level][i][j];
    }
  }
  
}

double compute_residual(double **num_sol, double **f, double **residual, int n_nodes){
    //! Compute the residual at each node and then take an L2 norm.
    //! Form R(num_sol) = 0 where R(num_sol) is the discrete laplacian + the source.

    double norm;
    double h2 = pow(1.0/((double)n_nodes-1.0),2.0);
    
    norm = 0.0;
    for (int i = 1; i < n_nodes-1; i++){
        for (int j = 1; j < n_nodes-1; j++){
            residual[i][j] = f[i][j] + (num_sol[i][j-1] + num_sol[i-1][j] + num_sol[i+1][j] + num_sol[i][j+1] - 4.0*num_sol[i][j])/h2;                                  
            norm += residual[i][j]*residual[i][j];
        }
    }
    
    norm = sqrt(norm);
  
    return norm;
}

void write_output(double **num_sol, double **exact_sol, double **x, double **y, int n_nodes, int iter, double residual, bool visualize){
  
    double error      = 0.0;
    double error_L2   = 0.0;
    double center_sol = 0.0;    
    int center = n_nodes/2+1;
    char path[100], file_name[100];
  
    strcpy(path, "./sol/solution");
    sprintf(file_name, "_%d.dat", iter+1);
    strcat(path, file_name);
    ofstream sol_file;
    sol_file.precision(10);
    sol_file.open(path, ios::out);    
    sol_file << "x, y, num_sol, exact_sol, l2norm";
    sol_file << endl;

    for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_nodes; j++) {
        sol_file << x[i][j] << " ";
        }
        sol_file << endl;
    }
  
    for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_nodes; j++) {
        sol_file << y[j][j]  << " ";
        }
        sol_file << endl;
    }
  
    for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_nodes; j++) {
        sol_file << num_sol[i][j]  << " ";
        if (i == center && j == center) center_sol = num_sol[i][j];
        }
        sol_file << endl;
    }
  
    for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_nodes; j++) {
        sol_file << exact_sol[i][j]  << " ";
        }
        sol_file << endl;
    }
  
    for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_nodes; j++) {
            error = exact_sol[i][j] - num_sol[i][j];
            sol_file << fabs(error) << " ";
        }
        sol_file << endl;
    }      
  
    for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_nodes; j++) {
            error = exact_sol[i][j] - num_sol[i][j];
            error_L2 += error * error;
        }
    }    
    error_L2 = sqrt(error_L2);

}

void deallocate_arrays(double ***num_sol, double **exact_sol, double ***f,
                       double ***x, double ***y, double ***mid_sol, int n_nodes,
                       int n_level) {
  
  //! Deallocation of all dynamic memory in the program.
  
  int nodes = n_nodes;
  for (int idx_level = 0; idx_level < n_level; idx_level++) {
    
    //! Delete j dimension
    
    for (int i = 0; i < nodes; i++) {
      delete [] x[idx_level][i];
      delete [] y[idx_level][i];
      delete [] num_sol[idx_level][i];
      delete [] mid_sol[idx_level][i];
      delete [] f[idx_level][i];
      if (idx_level == 0) delete [] exact_sol[i];
    }
    
    //! Delete i dimension
    
    delete [] x[idx_level];
    delete [] y[idx_level];
    delete [] num_sol[idx_level];
    delete [] mid_sol[idx_level];
    delete [] f[idx_level];
    if (idx_level == 0) delete [] exact_sol;
    
    //! Compute number of nodes on the next coarse level
    
    nodes = (nodes-1)/2 + 1;
    
  }
  
  //! Delete levels dimension
  
  delete [] x;
  delete [] y;
  delete [] num_sol;
  delete [] mid_sol;
  delete [] f;
  
}


/**************실행**************/
int main(){      
  
    bool visualize = VISUALIZE;
    bool stop_calc = false;
    int freq = FREQUENCY;
    int n_mgcycles = MG_CYCLES;
    int n_level;
    int n_smooth = NUM_SWEEP;
    int i_mgcycles = 0;
    int n_nodes = NUM_NODES;
    int mg_levels = MG_LEVELS;
    double tolerance = TOLERANCE;
    double residual_0;
    double residual;
    double residual_old;
            
    double ***num_sol;
    double **exact_sol;
    double ***f;
    double ***x;
    double ***y;
    double ***mid_sol;

    n_level = mg_levels;    
    allocate_arrays(&num_sol, &exact_sol, &f, &x, &y, &mid_sol, n_nodes); //메모리 할당
    generate_finest_grid(x[FINE_MESH], y[FINE_MESH], n_nodes, visualize); //가장 fine한 initial grid 생성          
    generate_coarse_gird(x, y, n_nodes, n_level, FINE_MESH, visualize);  //MG level에 맞는 coarse grid 생성    
    intitialize_solution(num_sol[FINE_MESH], exact_sol, f[FINE_MESH], x[FINE_MESH], y[FINE_MESH], n_nodes); //초기해, exact_sol 계산    
    multigrid_cycle(num_sol, f, mid_sol, n_nodes, n_smooth, n_level, FINE_MESH, tolerance); //multigrid 계산
    write_output(num_sol[FINE_MESH], exact_sol, x[FINE_MESH], y[FINE_MESH], n_nodes, i_mgcycles, residual, visualize); //result .dat로 저장
    deallocate_arrays(num_sol, exact_sol, f, x, y, mid_sol, n_nodes, n_level); //메모리 제거
    printf("Done!\n");

    return 0;
}