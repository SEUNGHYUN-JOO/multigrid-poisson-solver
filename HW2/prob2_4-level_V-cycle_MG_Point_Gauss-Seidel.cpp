#include <stdio.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <iostream>

using namespace std;


/**************변수 정의**************/
const int Num_node      = 17;    // 각 축 방향 노드 수
const int MG_level      = 4;     // V-cycle level 수
const int MG_max_CYCLE  = 1;     // 최대 MG V-cycle
const int level_fine    = 0;     // 가장 fine한 grid level
const int Smooth_per_GS = 3;     // GS 반복 횟수
const double Tolerance  = 12.0;  // tol 설정


/**************함수 정의**************/
void allocate_arrays(double ***num_sol, double **exact_sol, double ***source_term, double ***A_matrix, double ***x, double ***y, double ***mid_sol, int n_node);  //array 동적 메모리 할당
void create_finest_grid(double **x, double **y, int n_nodes);                                                                                       //fine grid 생성
void create_coarse_grid(double ***x, double ***y, int n_nodes, int n_level, int level);                                                             //coarse grid 생성
void intitialize_solution(double **num_sol, double **exact_sol, double **source_term, double **x, double **y, int n_nodes);                         //fine grid에서 초기해, exact_sol 설정
void compute_multigrid(double ***num_sol, double ***source_term, double ***mid_sol, int n_nodes, int n_smooth, int n_level, int level);             //MG 4-level V-cycle 수행
void gauss_seidel(double **num_sol, double **exact_sol, double **source_term, double **mid_sol, int n_nodes, int n_smooth, double tolerance);       //Gauss-Seidel로 smoothing
void restriction_step(double ***num_sol, double ***source_term, double ***mid_sol, int n_nodes, int level);                                         //restriction step - fine -> coarse grid
void prolongation_step(double ***num_sol, double ***mid_sol, int n_nodes, int level);                                                               //prolongate step - coarse -> fine grid
void compute_residual(double **num_sol, double **source_term, double **residual, int n_nodes);                                                                //residual 계산
void write_error(double **x, double **y, double **num_sol, double **exact_sol, int n_nodes, int iter);                                              //error 파일 작성
void write_result(double **x, double **y, double **num_sol, int n_nodes);                                                                           //result .dat 파일로 저장

void allocate_arrays(double ***num_sol, double **exact_sol, double ***source_term, double ***A_matrix, double ***x, double ***y, double ***mid_sol, int n_node){  
    int n_level = MG_level;    
    int num_elements;

    double **num_sol_temp;
    double *exact_sol_temp;
    double **mid_sol_temp;
    double **source_term_temp;
    double **x_temp; 
    double **y_temp;
    double **A_matrix_temp;

    // 레벨 수만큼 배열을 생성
    x_temp           = new double*[n_level];
    y_temp           = new double*[n_level];
    num_sol_temp     = new double*[n_level];
    mid_sol_temp     = new double*[n_level];
    source_term_temp = new double*[n_level];
    A_matrix_temp    = new double*[n_level];

    for (int idx_level = 0; idx_level < n_level; idx_level++){
        num_elements = n_node * n_node;  // 각 레벨에서 grid pts 수

        // 가장 fine한 레벨에서만 exact_sol 할당
        if (idx_level == 0){
            exact_sol_temp = new double[num_elements];
        }

        // 각 레벨에서 필요한 만큼 1차원 배열을 할당
        x_temp[idx_level]           = new double[num_elements];
        y_temp[idx_level]           = new double[num_elements];
        num_sol_temp[idx_level]     = new double[num_elements];
        mid_sol_temp[idx_level]     = new double[num_elements];
        source_term_temp[idx_level] = new double[num_elements];
        A_matrix_temp[idx_level]    = new double[num_elements * num_elements];

        // 다음 레벨에서는 노드 수를 절반으로 줄입니다
        n_node = (n_node - 1) / 2 + 1;
    }

    // 동적 메모리를 외부 포인터에 할당
    *x           = x_temp;
    *y           = y_temp;
    *num_sol     = num_sol_temp;
    *mid_sol     = mid_sol_temp;
    *source_term = source_term_temp;
    *exact_sol   = exact_sol_temp;
    *A_matrix    = A_matrix_temp;
}


void create_finest_grid(double *x, double *y, int n_node){
    // 0<=x,y<=1 범위에서 n_node를 갖는 가장 fine한 grid 생성
    for (int j = 0; j < n_node; j++){
        for (int i = 0; i < n_node; i++){
            int idx = j * n_node + i;
            x[idx] = (double) i/ (double) (n_node-1);
            y[idx] = (double) j/ (double) (n_node-1);
        }
    }
}


void create_coarse_grid(double **x, double **y, int n_node, int n_level, int level){
    // 각 coarse level에서 grid 생성
    if (level != n_level - 1) {
        for (int j = 0; j < n_node; j++) {
            for (int i = 0; i < n_node; i++) {
                int idx_fine = j * n_node + i;  // 상위 level의 1차원 인덱스

                // coarse grid는 상위 level의 i, j 값을 i/2, j/2에 가짐
                if (i % 2 == 0 && j % 2 == 0) {
                    int i_coarse = i / 2;
                    int j_coarse = j / 2;
                    int idx_coarse = j_coarse * ((n_node - 1) / 2 + 1) + i_coarse;  // coarse level에서 1차원 인덱스

                    x[level + 1][idx_coarse] = x[level][idx_fine];
                    y[level + 1][idx_coarse] = y[level][idx_fine];
                }
            }
        }

        // 가장 coarse한 grid가 될 때까지 반복
        int coarse_level = level + 1;
        int n_coarse_node = (n_node - 1) / 2 + 1;
        create_coarse_grid(x, y, n_coarse_node, n_level, coarse_level);
    }
}


void intitialize_solution(double *num_sol, double *exact_sol, double *source_term, double *x, double *y, int n_node){
    //가장 fine한 grid에서 정의하면 coarse grid는 할 필요 없음
    //문제 2번 BC가 Drichlet이므로 boundary만 값 주고 나머지는 0으로 초기화 (문제2는 boundary이 모두 0.0)
    float x_temp;
    float y_temp;

    for (int j = 0; j < n_node; j++){
        for (int i = 0; i < n_node; i++){
            
            int idx = j * n_node + i;            
            x_temp = x[idx];
            y_temp = y[idx];
            
            num_sol[idx] = 0.0;

            //문제 2번 포아송 방정식 및 exact sol 참고
            exact_sol[idx]   = -1 * x_temp*x_temp*y_temp*y_temp*(1-x_temp*x_temp)*(1-y_temp*y_temp);
            source_term[idx] = -2 * ((1-6*x_temp*x_temp)*y_temp*y_temp*(1-y_temp*y_temp) + (1-6*y_temp*y_temp)*x_temp*x_temp*(1-x_temp*x_temp));
        }
    }  
}

void compute_multigrid(double ***num_sol, double ***source_term, double ***mid_sol, int n_node, int n_smooth, int n_level, int level, double tolerance){
    //Gauss-Seidel로 n_smooth번 smoothing    
    gauss_seidel(num_sol[level], source_term[level], mid_sol[level], n_node, n_smooth);
    if (level < n_level-1){
        int n_coarse_nodes = (n_node-1)/(2) + 1;
        int level_coarse = level + 1;

        restriction_step(num_sol, source_term, mid_sol, n_node, level);            
        compute_multigrid(num_sol, source_term, mid_sol, n_coarse_nodes, n_smooth, n_level, level_coarse, tolerance);
        
        //level >= n_level-1이 되면 compute_multigrid에서 G-S만 수행하고 prolongation_step 진행해서 finest grid로 돌아옴
        prolongation_step(num_sol, mid_sol, n_node, level);
    }

    //prolongation 된걸로 최종 G-S 수행 -> get Sol
    gauss_seidel(num_sol[level], source_term[level], mid_sol[level], n_node, n_smooth);
}

void gauss_seidel(double *num_sol, double *source_term, int n_node, int n_smooth){
    double dx = 1.0/((double)n_node-1.0);
    double dy = 1.0/((double)n_node-1.0);
    double beta = dy/dx;
    double coeff = 1.0/(2.0*(1+beta*beta));
    //double sigma = (cos(M_PI*dx) + beta*beta*cos(M_PI*dy))/(1+beta*beta); //for Gauss-Seidel over relaxation
    //double weight = 2.0 / (1.0 + sqrt(1-sigma*sigma));                    //for Gauss-Seidel over relaxation
        
    for (int iter = 0; iter < n_smooth; iter++){        
        for (int j = 1; j < n_node-1; j++){
            for (int i = 1; i < n_node-1; i++){
                int idx = j * n_node + i;
                
                int idx_left = j * n_node + (i - 1);
                int idx_right = j * n_node + (i + 1);
                int idx_up = (j - 1) * n_node + i;
                int idx_down = (j + 1) * n_node + i;

                num_sol[idx] = coeff * (num_sol[idx_down] + num_sol[idx_up] + beta * beta * (num_sol[idx_right] + num_sol[idx_left]));
            }
        }
    }
}

void restriction_step(double ***num_sol, double ***source_term, double ***mid_sol, int n_node, int level){  
    
    compute_residual(num_sol[level], source_term[level], mid_sol[level], n_node);  
    
    int n_coarse_node = (n_node-1)/2 + 1;
    int level_coarse   = level+1;
    int i_idx_fine;
    int j_idx_fine;
  
    //coarse grid에서 guess sol, source term(grid func) 초기화
    for (int i = 0; i < n_coarse_node; i++){
        for (int j = 0; j < n_coarse_node; j++){
            num_sol[level_coarse][i][j]     = 0.0;
            source_term[level_coarse][i][j] = 0.0;
        }
    }

    //BC 제외한 grid points에서 주변 9 points stencil에서 수행 / ref: Multigrid P42 / source_term=grid function (guess) / d2h = I_2h_h * dh
    for (int i = 1; i < n_coarse_node-1; i++){
        for (int j = 1; j < n_coarse_node-1; j++){                  
            i_idx_fine = 2*i;
            j_idx_fine = 2*j;
            source_term[level_coarse][i][j] = (1.0/16.0)*(1.0*mid_sol[level][i_idx_fine-1][j_idx_fine-1] + 2.0*mid_sol[level][i_idx_fine-1][j_idx_fine] + 1.0*mid_sol[level][i_idx_fine-1][j_idx_fine+1]
                                                         +2.0*mid_sol[level][i_idx_fine  ][j_idx_fine-1] + 4.0*mid_sol[level][i_idx_fine  ][j_idx_fine] + 2.0*mid_sol[level][i_idx_fine  ][j_idx_fine+1]
                                                         +1.0*mid_sol[level][i_idx_fine+1][j_idx_fine-1] + 2.0*mid_sol[level][i_idx_fine+1][j_idx_fine] + 1.0*mid_sol[level][i_idx_fine+1][j_idx_fine+1]);
        }
    }

    //다시 gauss-seidel로 풀기
}

void prolongation_step(double ***num_sol, double ***mid_sol, int n_node, int level){
    int n_coarse_node = (n_node-1)/2 + 1;
    int level_coarse = level+1;
    int i_idx_fine;
    int j_idx_fine;
    
    //mid_sol 초기화
    for (int i = 0; i < n_node; i++){
        for (int j = 0; j < n_node; j++){
            mid_sol[level][i][j] = 0.0;
        }
    }

    //fine grid에서 coarse grid num_sol로 mid_sol update(residual) interpolation
    for (int i = 1; i < n_coarse_node-1; i++){
        for (int j = 1; j < n_coarse_node-1; j++){
            i_idx_fine = i*2;
            j_idx_fine = j*2;
            
            //ref: Multigrid P43~44 / weight 대각은 1/4, 변은 1/2, 자기자신은 그대로
            mid_sol[level][i_idx_fine-1][j_idx_fine-1] += 1.0/4.0 * num_sol[level_coarse][i][j];
            mid_sol[level][i_idx_fine-1][j_idx_fine  ] += 1.0/2.0 * num_sol[level_coarse][i][j];
            mid_sol[level][i_idx_fine-1][j_idx_fine+1] += 1.0/4.0 * num_sol[level_coarse][i][j];
            mid_sol[level][i_idx_fine  ][j_idx_fine-1] += 1.0/2.0 * num_sol[level_coarse][i][j];
            mid_sol[level][i_idx_fine  ][j_idx_fine  ]  = num_sol[level_coarse][i][j];
            mid_sol[level][i_idx_fine  ][j_idx_fine+1] += 1.0/2.0 * num_sol[level_coarse][i][j];
            mid_sol[level][i_idx_fine+1][j_idx_fine-1] += 1.0/4.0 * num_sol[level_coarse][i][j];
            mid_sol[level][i_idx_fine+1][j_idx_fine  ] += 1.0/2.0 * num_sol[level_coarse][i][j];
            mid_sol[level][i_idx_fine+1][j_idx_fine+1] += 1.0/4.0 * num_sol[level_coarse][i][j];
        }
    }

    //coarse residul를 fine grid num_sol에 더해줌
    for (int i = 0; i < n_node; i++){
        for (int j = 0; j < n_node; j++){
            num_sol[level][i][j] += mid_sol[level][i][j];
        }
    }
}


void compute_residual(double **num_sol, double **source_term_temp, double **residual, int n_node){
    //5 points stencil로 residual 계산 
    //residual = source_term - laplacian(num_sol)
    double h_square = ((double)n_node-1.0)*((double)n_node-1.0);
    
    //ref: Multigrid / P42 Lh 참고
    for (int i = 1; i < n_node-1; i++){
        for (int j = 1; j < n_node-1; j++){
            residual[i][j] = source_term_temp[i][j] - (- 1.0*num_sol[i-1][j] - 1.0*num_sol[i+1][j] - 1.0*num_sol[i][j-1] - 1.0*num_sol[i][j+1] + 4.0*num_sol[i][j])/h_square;
        }
    }
}

void write_error(double **x, double **y, double **num_sol, double **exact_sol, int n_node, int iter){
  
    char file_name[100];
    double error  = 0.0;
    double L2norm = 0.0;
    
    strcpy(file_name, "./error_MG_V_cycle");    
    ofstream error_file;
    error_file.precision(10);
    error_file.open(file_name, ios::out);
    
    if (iter == 0){
        error_file << "iter, l2norm";
        error_file << endl;
    }
      
    for (int i = 0; i < n_node; i++) {
        for (int j = 0; j < n_node; j++) {
            error = exact_sol[i][j] - num_sol[i][j];
            L2norm += error * error;
        }
    }

    L2norm = sqrt(L2norm);

    error_file << iter+1 << " | " << L2norm << " ";
    error_file << endl;
}

void write_result(double **x, double **y, double **num_sol, int n_node){

    char file_name[100];

    strcpy(file_name, "./num_sol_MG_V_cycle");        
    ofstream sol_file;
    sol_file.precision(10);
    sol_file.open(file_name, ios::out);    
    sol_file << "x, y, num_sol";
    sol_file << endl;

    for (int i = 0; i < n_node; i++){
        for (int j = 0; j < n_node; j++){
            sol_file << x[i][j] << " ";
        }
        sol_file << endl;
    }
  
    for (int i = 0; i < n_node; i++){
        for (int j = 0; j < n_node; j++){
            sol_file << y[j][j]  << " ";
        }
        sol_file << endl;
    }
  
    for (int i = 0; i < n_node; i++){
        for (int j = 0; j < n_node; j++){
            sol_file << num_sol[i][j]  << " ";
        }
        sol_file << endl;
    }
}

/**************실행**************/
int main(){      
    int n_node = Num_node;
    int n_level = MG_level;
    int idx_vcycle;
    int n_smooth = Smooth_per_GS;
    double tolerance = Tolerance;
    double residual_0;
    double residual;
    double residual_old;

    double ***num_sol;
    double **exact_sol;
    double ***f;
    double ***x;
    double ***y;
    double ***mid_sol;

    allocate_arrays();              //메모리 할당  allocate_arrays(&num_sol, &exact_sol, &f, &x, &y, &mid_sol, n_nodes);     
    create_finest_grid(x[level_fine], y[level_fine], n_node);   //가장 fine한 initial grid 생성              
    create_coarse_grid()  //MG level에 맞는 coarse grid 생성                        
    intitialize_solution(num_sol[level_fine], exact_sol, f[level_fine], x[level_fine], y[level_fine], n_node);       //초기해, exact_sol 계산
  
    for (idx_vcycle = 1; idx_vcycle <= MG_max_CYCLE; idx_vcycle++){
        compute_multigrid(num_sol, f, mid_sol, n_node, n_smooth, n_level, level_fine, tolerance);   
        write_error(x[level_fine], y[level_fine], num_sol[level_fine], exact_sol, n_node, iter);   

        if (//L2norm <= pow((1.0/10.0),tol){
            write_result(x[level_fine]   , y[level_fine]   , num_sol[level_fine]   , n_node);   
            printf("Converged in %d's cycle", idx_vcycle);
        }
    }
    return 0;
}