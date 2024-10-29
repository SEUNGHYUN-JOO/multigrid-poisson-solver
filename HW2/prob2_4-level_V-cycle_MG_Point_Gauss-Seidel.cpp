#include <stdio.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <iostream>

using namespace std;


/**************변수 정의**************/
const int Num_node      = 9;    // 각 축 방향 노드 수
const int MG_level      = 4;     // V-cycle level 수
const int MG_max_CYCLE  = 30;     // 최대 MG V-cycle
const int level_fine    = 0;     // 가장 fine한 grid level
const int Smooth_per_GS = 100;     // GS 반복 횟수
const double Tolerance  = 5.0;  // tol 설정


/**************함수 정의**************/
void allocate_arrays(double ***num_sol_new, double ***num_sol_old, double **exact_sol, double ***source_term, double ***A_matrix, double ***x, double ***y, double ***mid_sol, int n_node);  //array 동적 메모리 할당
void create_finest_grid(double *x, double *y, double *A_matrix, int n_node);                                                                                       //fine grid 생성
void create_coarse_grid(double **x, double **y, double **A_matrix, int n_node, int n_level, int level);                                                             //coarse grid 생성
void intitialize_solution(double *num_sol_new, double *num_sol_old, double *exact_sol, double *source_term, double *x, double *y, int n_node);                         //fine grid에서 초기해, exact_sol 설정
double compute_multigrid(double **num_sol_new, double **num_sol_old, double **source_term, double **A_matrix, double **mid_sol, int n_node, int n_smooth, int n_level, int level);             //MG 4-level V-cycle 수행
double gauss_seidel(double *num_sol_new, double *num_sol_old, double *source_term, double *A_matrix, int n_node, int n_smooth);         //Gauss-Seidel로 smoothing
void restriction_step(double **num_sol_new, double **source_term, double **mid_sol, int n_node, int level);                                         //restriction step - fine -> coarse grid
void prolongation_step(double **num_sol_new, double **num_sol_old, double **mid_sol, int n_node, int level);                                                               //prolongate step - coarse -> fine grid
void compute_residual(double *num_sol, double *source_term, double *residual, int n_node);                                                                //residual 계산
void write_error(int iter, double L2norm);                                                                                                  //error 파일 작성
void write_result(double **x, double **y, double **num_sol, int n_node);                                                                           //result .dat 파일로 저장

void allocate_arrays(double ***num_sol_new, double ***num_sol_old, double **exact_sol, double ***source_term, double ***A_matrix, double ***x, double ***y, double ***mid_sol, int n_node){  
    int n_level = MG_level;    
    int num_elements;

    double **num_sol_new_temp;
    double **num_sol_old_temp;
    double *exact_sol_temp;
    double **mid_sol_temp;
    double **source_term_temp;
    double **x_temp; 
    double **y_temp;
    double **A_matrix_temp;

    // 레벨 수만큼 배열을 생성
    x_temp           = new double*[n_level];
    y_temp           = new double*[n_level];
    num_sol_new_temp = new double*[n_level];
    num_sol_old_temp = new double*[n_level];
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
        num_sol_new_temp[idx_level] = new double[num_elements];
        num_sol_old_temp[idx_level] = new double[num_elements];
        mid_sol_temp[idx_level]     = new double[num_elements];
        source_term_temp[idx_level] = new double[num_elements];
        A_matrix_temp[idx_level]    = new double[num_elements * num_elements];

        // 다음 레벨에서는 노드 수를 절반으로 줄입니다
        n_node = (n_node - 1) / 2 + 1;
    }

    // 동적 메모리를 외부 포인터에 할당
    *x           = x_temp;
    *y           = y_temp;
    *num_sol_new = num_sol_new_temp; //u_new
    *num_sol_old = num_sol_old_temp; //u_old
    *mid_sol     = mid_sol_temp;
    *source_term = source_term_temp; //b
    *exact_sol   = exact_sol_temp;   //u_exact
    *A_matrix    = A_matrix_temp;    //m
}


void create_finest_grid(double *x, double *y, double *A_matrix, int n_node){
    // 0<=x,y<=1 범위에서 n_node를 갖는 가장 fine한 grid 생성
    for (int j = 0; j < n_node; j++){
        for (int i = 0; i < n_node; i++){
            int idx = j * n_node + i;
            x[idx] = (double) i/ (double) (n_node-1);
            y[idx] = (double) j/ (double) (n_node-1);

            
            A_matrix[idx * n_node * n_node + idx] = -4.0;

            if (i>0){
                int idx_left = j * n_node + (i - 1);
                A_matrix[idx * n_node * n_node + idx_left] = 1.0;
                A_matrix[idx_left * n_node * n_node + idx] = 1.0;
            }
            if (j>0){
                int idx_up = (j-1) * n_node + i;
                A_matrix[idx * n_node * n_node + idx_up] = 1.0;
                A_matrix[idx_up * n_node * n_node + idx] = 1.0;
            }
        }
    }
}


void create_coarse_grid(double **x, double **y, double **A_matrix, int n_node, int n_level, int level){
    //double dx = 1.0/((double)n_node-1.0);
    //double dy = 1.0/((double)n_node-1.0);
    //double beta = dy/dx;
    //double coeff = -(2.0*(1+beta*beta));
    //double sigma = (cos(M_PI*dx) + beta*beta*cos(M_PI*dy))/(1+beta*beta); //for Gauss-Seidel over relaxation
    //double weight = 2.0 / (1.0 + sqrt(1-sigma*sigma));                    //for Gauss-Seidel over relaxation
    
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
                    
                    //printf("x[%d] = %f, y[%d] = %f\n", i, x[level + 1][idx_coarse], i, y[level + 1][idx_coarse]);

                }
            }
        }

        // coarse level에서 A_matrix 설정
        int n_coarse_node = (n_node - 1) / 2 + 1;
        for (int j = 0; j < n_coarse_node; j++) {
            for (int i = 0; i < n_coarse_node; i++) {
                int idx_coarse = j * n_coarse_node + i;

                // 대각 요소 설정
                A_matrix[level + 1][idx_coarse * n_coarse_node * n_coarse_node + idx_coarse] = -4.0;

                // 왼쪽 요소 설정
                if (i > 0) {
                    int idx_left = j * n_coarse_node + (i - 1);
                    A_matrix[level + 1][idx_coarse * n_coarse_node * n_coarse_node + idx_left] = 1.0;
                    A_matrix[level + 1][idx_left * n_coarse_node * n_coarse_node + idx_coarse] = 1.0;
                }

                // 위쪽 요소 설정
                if (j > 0) {
                    int idx_up = (j - 1) * n_coarse_node + i;
                    A_matrix[level + 1][idx_coarse * n_coarse_node * n_coarse_node + idx_up] = 1.0;
                    A_matrix[level + 1][idx_up * n_coarse_node * n_coarse_node + idx_coarse] = 1.0;
                }
            }
        }        
        
        for (int j = 0; j < n_coarse_node; j++) {
            for (int i = 0; i < n_coarse_node; i++) {
                int idx_coarse = j * n_coarse_node + i;
                printf("%8.3f ", A_matrix[level + 1][idx_coarse]);  // 각 요소를 8.3 형식으로 출력
            }
            printf("\n");  // 각 행이 끝날 때 줄바꿈
        }
        printf("\n");

        // 다음 coarse 레벨로 재귀 호출
        int coarse_level = level + 1;
        create_coarse_grid(x, y, A_matrix, n_coarse_node, n_level, coarse_level);
    }
}


void intitialize_solution(double *num_sol_new, double *num_sol_old, double *exact_sol, double *source_term, double *x, double *y, int n_node){
    //가장 fine한 grid에서 정의하면 coarse grid는 할 필요 없음
    //문제 2번 BC가 Drichlet이므로 boundary만 값 주고 나머지는 0으로 초기화 (문제2는 boundary이 모두 0.0)
    float x_temp;
    float y_temp;

    for (int j = 0; j < n_node; j++){
        for (int i = 0; i < n_node; i++){
            
            int idx = j * n_node + i;            
            x_temp = x[idx];
            y_temp = y[idx];
            
            num_sol_new[idx] = 0.0;
            num_sol_old[idx] = 0.0;            

            //문제 2번 포아송 방정식 및 exact sol 참고
            exact_sol[idx]   = -1 * x_temp*x_temp*y_temp*y_temp*(1-x_temp*x_temp)*(1-y_temp*y_temp);
            source_term[idx] = -2 * ((1-6*x_temp*x_temp)*y_temp*y_temp*(1-y_temp*y_temp) + (1-6*y_temp*y_temp)*x_temp*x_temp*(1-x_temp*x_temp));
        }
    }
}

double compute_multigrid(double **num_sol_new, double **num_sol_old, double **source_term, double **A_matrix, double **mid_sol, int n_node, int n_smooth, int n_level, int level){
    double L2norm;
    //Gauss-Seidel로 n_smooth번 smoothing            
    L2norm = gauss_seidel(num_sol_new[level], num_sol_old[level], source_term[level], A_matrix[level], n_node, n_smooth);
    if (level < n_level-1){
        int n_coarse_nodes = (n_node-1)/(2) + 1;
        int level_coarse = level + 1;

        restriction_step(num_sol_new, source_term, mid_sol, n_node, level);            
        compute_multigrid(num_sol_new, num_sol_old, source_term, A_matrix, mid_sol, n_coarse_nodes, n_smooth, n_level, level_coarse);
        
        //level >= n_level-1이 되면 compute_multigrid에서 G-S만 수행하고 prolongation_step 진행해서 finest grid로 돌아옴
        prolongation_step(num_sol_new, num_sol_old, mid_sol, n_node, level);
    }

    //prolongation 된걸로 최종 G-S 수행 -> get Sol
    L2norm = gauss_seidel(num_sol_new[level], num_sol_old[level], source_term[level], A_matrix[level], n_node, n_smooth);

    return L2norm;
}


double gauss_seidel(double *num_sol_new, double *num_sol_old, double *source_term, double *A_matrix, int n_node, int n_smooth){
    double error;

    for (int iter = 0; iter < n_smooth; iter++){
        error = 0.0;
        
        // Gauss-Seidel 반복
        for (int j = 1; j < n_node - 1; j++) {
            for (int i = 1; i < n_node - 1; i++){
                int idx = j * n_node + i;

                int idx_left = j * n_node + (i - 1);
                int idx_right = j * n_node + (i + 1);
                int idx_up = (j - 1) * n_node + i;
                int idx_down = (j + 1) * n_node + i;

                num_sol_new[idx] = (1.0 / A_matrix[idx * n_node * n_node + idx]) * (source_term[idx] 
                                                                                - A_matrix[idx * n_node * n_node + idx_up   ] * num_sol_old[idx_up]
                                                                                - A_matrix[idx * n_node * n_node + idx_left ] * num_sol_old[idx_left]
                                                                                - A_matrix[idx * n_node * n_node + idx_right] * num_sol_old[idx_right]
                                                                                - A_matrix[idx * n_node * n_node + idx_down ] * num_sol_old[idx_down]);
            }
        }        

        // 오류 계산 (L2 norm)
        double norm = 0.0;
        for (int k = 0; k < n_node * n_node; k++) {
            double diff = num_sol_new[k] - num_sol_old[k];
            norm += diff * diff;
        }
        error = sqrt(norm);

        // num_sol_old 업데이트
        for (int k = 0; k < n_node * n_node; k++) {
            num_sol_old[k] = num_sol_new[k];
        }
    }
    return error;
}

void restriction_step(double **num_sol_new, double **source_term, double **mid_sol, int n_node, int level){  
    
    compute_residual(num_sol_new[level], source_term[level], mid_sol[level], n_node);  
    
    int n_coarse_node = (n_node-1)/2 + 1;
    int level_coarse   = level+1;
    int i_idx_fine;
    int j_idx_fine;
  
    //coarse grid에서 guess sol, source term(grid func) 초기화
    for (int j = 0; j < n_coarse_node; j++){
        for (int i = 0; i < n_coarse_node; i++){
            int idx_coarse = j * n_coarse_node + i;
            source_term[level_coarse][idx_coarse] = 0.0;
        }
    }

    //BC 제외한 grid points / ref: Multigrid P42 / source_term=grid function (guess) / d2h = I_2h_h * dh
    for (int j = 1; j < n_coarse_node-1; j++){
        for (int i = 1; i < n_coarse_node-1; i++){
            int i_idx_fine = 2*i;
            int j_idx_fine = 2*j;
            int idx_coarse = j * n_coarse_node + i;
            int idx_fine   = j_idx_fine * n_node + i_idx_fine;

            source_term[level_coarse][idx_coarse] = (1.0/16.0)*(1.0*mid_sol[level][idx_fine - n_node - 1] + 2.0*mid_sol[level][idx_fine - n_node] + 1.0*mid_sol[level][idx_fine - n_node + 1]
                                                               +2.0*mid_sol[level][idx_fine - 1         ] + 4.0*mid_sol[level][idx_fine         ] + 2.0*mid_sol[level][idx_fine + 1]
                                                               +1.0*mid_sol[level][idx_fine + n_node - 1] + 2.0*mid_sol[level][idx_fine + n_node] + 1.0*mid_sol[level][idx_fine + n_node + 1]);
        }
    }

    //다시 gauss-seidel로 풀기
}

void prolongation_step(double **num_sol_new, double **num_sol_old, double **mid_sol, int n_node, int level){
    int n_coarse_node = (n_node-1)/2 + 1;
    int level_coarse = level + 1;
    
    //mid_sol 초기화
    for (int i = 0; i < n_node * n_node; i++){        
        mid_sol[level][i] = 0.0;
    }

    //fine grid에서 coarse grid num_sol로 mid_sol update(residual) interpolation
    for (int j = 1; j < n_coarse_node-1; j++){
        for (int i = 1; i < n_coarse_node-1; i++){
            int i_idx_fine = i*2;
            int j_idx_fine = j*2;
            int idx_fine   = j_idx_fine * n_node + i_idx_fine;
            int idx_coarse = j * n_coarse_node + i;
            
            //ref: Multigrid P43~44 / weight 대각은 1/4, 변은 1/2, 자기자신은 그대로
            mid_sol[level][idx_fine - n_node - 1] += 1.0/4.0 * num_sol_new[level_coarse][idx_coarse];
            mid_sol[level][idx_fine - n_node    ] += 1.0/2.0 * num_sol_new[level_coarse][idx_coarse];
            mid_sol[level][idx_fine - n_node + 1] += 1.0/4.0 * num_sol_new[level_coarse][idx_coarse];
            mid_sol[level][idx_fine - 1         ] += 1.0/2.0 * num_sol_new[level_coarse][idx_coarse];
            mid_sol[level][idx_fine             ]  = num_sol_new[level_coarse][idx_coarse];
            mid_sol[level][idx_fine + 1         ] += 1.0/2.0 * num_sol_new[level_coarse][idx_coarse];
            mid_sol[level][idx_fine + n_node - 1] += 1.0/4.0 * num_sol_new[level_coarse][idx_coarse];
            mid_sol[level][idx_fine + n_node    ] += 1.0/2.0 * num_sol_new[level_coarse][idx_coarse];
            mid_sol[level][idx_fine + n_node + 1] += 1.0/4.0 * num_sol_new[level_coarse][idx_coarse];
        }
    }

    //coarse residul를 fine grid num_sol에 더해줌
    for (int i = 0; i < n_node * n_node; i++){        
        num_sol_old[level][i] += mid_sol[level][i];
    }
}


void compute_residual(double *num_sol_new, double *source_term_temp, double *residual, int n_node){
    //5 points stencil로 residual 계산 
    //residual = source_term - laplacian(num_sol)
    double h_square = ((double)n_node-1.0)*((double)n_node-1.0);
    
    //ref: Multigrid / P42 Lh 참고
    for (int j = 1; j < n_node-1; j++){
        for (int i = 1; i < n_node-1; i++){
            int idx = j * n_node + i;
            residual[idx] = source_term_temp[idx] - (- 1.0*num_sol_new[idx - n_node] - 1.0*num_sol_new[idx + n_node] - 1.0*num_sol_new[idx - 1] - 1.0*num_sol_new[idx + 1] + 4.0*num_sol_new[idx])/h_square;
        }
    }
}

void write_error(int iter, double L2norm){
  
    char file_name[100];        
    
    strcpy(file_name, "./error_MG_V_cycle.dat");    
    ofstream error_file;
    error_file.precision(10);
    if (iter == 1){
        error_file.open(file_name, ios::out);
    }
    else{
        error_file.open(file_name, ios::app);
    }
    if (iter == 0){
        error_file << "iter, l2norm";
        error_file << endl;
    }
      
    error_file << iter << " / " << L2norm << " ";
    error_file << endl;
}

void write_result(double **x, double **y, double **num_sol_new, int n_node){

    char file_name[100];

    strcpy(file_name, "./num_sol_MG_V_cycle.dat");        
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
            sol_file << num_sol_new[i][j]  << " ";
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
    double L2norm;

    // 배열 포인터들 선언
    double **num_sol_new;
    double **num_sol_old;
    double *exact_sol;
    double **source_term;
    double **x;
    double **y;
    double **mid_sol;
    double **A_matrix;

    // 메모리 할당
    allocate_arrays(&num_sol_new, &num_sol_old, &exact_sol, &source_term, &A_matrix, &x, &y, &mid_sol, n_node);

    // 가장 fine한 initial grid 생성
    create_finest_grid(x[level_fine], y[level_fine], A_matrix[level_fine], n_node);

    // MG level에 맞는 coarse grid 생성
    create_coarse_grid(x, y, A_matrix, n_node, n_level, level_fine);

    // 초기해 및 exact_sol 계산
    intitialize_solution(num_sol_new[level_fine], num_sol_old[level_fine], exact_sol, source_term[level_fine], x[level_fine], y[level_fine], n_node);

    for (idx_vcycle = 1; idx_vcycle <= MG_max_CYCLE; idx_vcycle++) {
        // Multi-grid 수행
        L2norm = compute_multigrid(num_sol_new, num_sol_old, source_term, A_matrix, mid_sol, n_node, n_smooth, n_level, level_fine);

        // 오차 출력
        write_error(idx_vcycle, L2norm);

        // 수렴 체크
        if (L2norm <= pow(10.0, -tolerance)){
            // 최종 결과 저장
            write_result(x, y, num_sol_new, n_node);
            printf("Converged in %d cycle(s)\n", idx_vcycle);            
            break;
        }
    }
    if (L2norm > pow(10.0, -tolerance)){
        printf("Failed to converge\n");        
    }
    system("pause");
    return 0;
}