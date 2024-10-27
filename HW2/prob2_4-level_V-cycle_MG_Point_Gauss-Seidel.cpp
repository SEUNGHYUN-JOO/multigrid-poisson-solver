#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/*Finest grid level 정의*/
const int FINE_GRID = 0;

/*Parameter 지정*/
struct SimulationConfig {
    int grid_size = 17;
    int max_cycles = 4;
    int smoothing_sweeps = 3;
    bool multigrid_enabled = true;    
    bool output_enabled = true;
    int output_frequency = 1;
    double tolerance = 1e-10;
    double relaxation_factor = 1.7;
} config;


/*변수 및 함수 형태 정의*/
int allocateMemory(double ****solution, double ***exact_solution, double ****source, double ****x_vals, double ****y_vals, double ****temp, int grid_nodes);
double calculateResidual(double **solution, double **source, double **residual, int grid_size);
void displayConfig(const SimulationConfig &config);
void setupFineGrid(double **x_vals, double **y_vals, int grid_size, bool output_flag);
void generateCoarseLevels(double ***x_vals, double ***y_vals, int grid_size, int levels, int current_level, bool output_flag);
void initializeProblem(double **solution, double **exact_solution, double **source, double **x_vals, double **y_vals, int grid_size);
void multigridCycle(double ***solution, double ***source, double ***temp, int grid_size, int sweep_count, int levels, int current_level);
void gaussSeidelSmoothing(double **solution, double **source, double **temp, int grid_size, int sweep_count);
void applyRestriction(double ***solution, double ***source, double ***temp, int grid_size, int level);
void applyProlongation(double ***solution, double ***temp, int grid_size, int level);


/*Main*/
int main(int argc, char* argv[]) {

    if (argc > 1) config.grid_size = atoi(argv[1]);
    if (argc > 2) config.max_cycles = atoi(argv[2]);

    bool stop_calculation = false;
    int levels, cycle_count = 0;
    double initial_residual, residual, previous_residual;

    double ***solution, **exact_solution, ***source, ***x_vals, ***y_vals, ***temp;

    levels = allocateMemory(&solution, &exact_solution, &source, &x_vals, &y_vals, &temp, config.grid_size);
    displayConfig(config);
    setupFineGrid(x_vals[0], y_vals[0], config.grid_size, config.output_enabled);
    generateCoarseLevels(x_vals, y_vals, config.grid_size, levels, 0, config.output_enabled);
    initializeProblem(solution[0], exact_solution, source[0], x_vals[0], y_vals[0], config.grid_size);
    initial_residual = calculateResidual(solution[0], source[0], temp[0], config.grid_size);

    for (cycle_count = 1; cycle_count <= config.max_cycles; cycle_count++) {
        multigridCycle(solution, source, temp, config.grid_size, config.smoothing_sweeps, levels, 0);
        residual = calculateResidual(solution[0], source[0], temp[0], config.grid_size);

        if (fabs(residual - previous_residual) < config.tolerance) stop_calculation = true;
        previous_residual = residual;

        if (stop_calculation) break;
    }
    return 0;
}

/* 초기 설정 출력을 위한 함수 */
void displayConfig(const SimulationConfig &config) {
    double spacing = 1.0 / (config.grid_size - 1);
    printf("## Poisson Equation Solver with Multigrid ##\n");
    printf("Grid Size: %d x %d\n", config.grid_size, config.grid_size);
    printf("Multigrid Levels: %d\n", config.multigrid_enabled ? "Enabled" : "Disabled");
    printf("Max Cycles: %d\n", config.max_cycles);
    printf("Grid Spacing: %e\n\n", spacing);
}


/*초기 메모리 할당*/
int allocateMemory(double ****solution, double ***exact_solution,
                   double ****source, double ****x_vals, double ****y_vals,
                   double ****temp, int grid_nodes) {
    
    int level_count = 1;
    bool can_coarsen = true;
    int current_nodes = grid_nodes;

    while (can_coarsen) {
        if (((current_nodes - 1) % 2 == 0) && ((current_nodes - 1) / 2 + 1 >= 5)) {
            current_nodes = (current_nodes - 1) / 2 + 1;
            level_count++;
        } else {
            can_coarsen = false;
        }
    }

    *x_vals = new double**[level_count];
    *y_vals = new double**[level_count];
    *solution = new double**[level_count];
    *temp = new double**[level_count];
    *source = new double**[level_count];

    current_nodes = grid_nodes;
    for (int level = 0; level < level_count; level++) {
        (*x_vals)[level] = new double*[current_nodes];
        (*y_vals)[level] = new double*[current_nodes];
        (*solution)[level] = new double*[current_nodes];
        (*temp)[level] = new double*[current_nodes];
        (*source)[level] = new double*[current_nodes];

        if (level == 0) *exact_solution = new double*[current_nodes];

        for (int i = 0; i < current_nodes; i++) {
            (*x_vals)[level][i] = new double[current_nodes];
            (*y_vals)[level][i] = new double[current_nodes];
            (*solution)[level][i] = new double[current_nodes];
            (*temp)[level][i] = new double[current_nodes];
            (*source)[level][i] = new double[current_nodes];

            if (level == 0) (*exact_solution)[i] = new double[current_nodes];
        }

        current_nodes = (current_nodes - 1) / 2 + 1;
    }
    return level_count;
}

/*초기 그리드 생성*/
void setupFineGrid(double **x_vals, double **y_vals, int grid_size, bool output_flag) {
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            x_vals[i][j] = (double)i / (double)(grid_size - 1);
            y_vals[i][j] = (double)j / (double)(grid_size - 1);
        }
    }
    printf("Fine grid generated (%d x %d)...\n", grid_size, grid_size);
}

/*초기화*/
void initializeProblem(double **solution, double **exact_solution, double **source,
                       double **x_vals, double **y_vals, int grid_size) {
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            if (i == 0 || j == 0 || i == grid_size - 1 || j == grid_size - 1) {
                solution[i][j] = log(1 + sin(M_PI * x_vals[i][j] * x_vals[i][j])) *
                                 (cos(sin(x_vals[i][j])) - 1) * sin(M_PI * y_vals[i][j]);
            } else {
                solution[i][j] = 0.0;
            }
            exact_solution[i][j] = log(1 + sin(M_PI * x_vals[i][j] * x_vals[i][j])) *
                                   (cos(sin(x_vals[i][j])) - 1) * sin(M_PI * y_vals[i][j]);
        }
    }
    printf("Problem initialized with boundary conditions...\n");
}

/*Multi-grid 사이클*/
void multigridCycle(double ***solution, double ***source, double ***temp,
                    int grid_size, int sweep_count, int levels, int current_level) {
    
    // 현재 레벨에서 잔차 계산 및 출력
    double residual = calculateResidual(solution[current_level], source[current_level], temp[current_level], grid_size);
    printf("Level %d Residual before smoothing: %e\n", current_level, residual);

    // 스무딩 연산 수행
    gaussSeidelSmoothing(solution[current_level], source[current_level], temp[current_level], grid_size, sweep_count);
    
    // 마지막 레벨이 아니라면, restriction을 적용하여 다음 레벨로 이동
    if (current_level < levels - 1) {
        applyRestriction(solution, source, temp, grid_size, current_level);
        int coarser_grid_size = (grid_size - 1) / 2 + 1;
        
        // 다음 레벨에서 multigridCycle 재귀 호출
        multigridCycle(solution, source, temp, coarser_grid_size, sweep_count, levels, current_level + 1);
        
        // 보간(prolongation)을 통해 다시 상위 레벨로 이동
        applyProlongation(solution, temp, grid_size, current_level);
    }

    // 스무딩 후 잔차 계산 및 출력
    residual = calculateResidual(solution[current_level], source[current_level], temp[current_level], grid_size);
    printf("Level %d Residual after smoothing: %e\n", current_level, residual);
    
    // 후처리로 다시 스무딩 연산 수행
    gaussSeidelSmoothing(solution[current_level], source[current_level], temp[current_level], grid_size, sweep_count);
    
    // 스무딩 후 최종 잔차 계산 및 출력
    residual = calculateResidual(solution[current_level], source[current_level], temp[current_level], grid_size);
    printf("Level %d Final Residual: %e\n", current_level, residual);
}

/*Gauss-Seidel smoothing 함수*/
void gaussSeidelSmoothing(double **solution, double **source, double **temp,
                          int grid_size, int sweep_count) {
    double h2 = pow(1.0 / (grid_size - 1), 2.0);
    for (int iter = 0; iter < sweep_count; iter++) {
        for (int i = 1; i < grid_size - 1; i++) {
            for (int j = 1; j < grid_size - 1; j++) {
                solution[i][j] = (solution[i][j - 1] + solution[i - 1][j] +
                                  solution[i + 1][j] + solution[i][j + 1] + h2 * source[i][j]) / 4.0;
            }
        }
    }
}

/*Coars level grid 만들기*/
void generateCoarseLevels(double ***x_vals, double ***y_vals, int grid_size, 
                          int levels, int current_level, bool output_flag) {
    if (current_level == levels - 1) {
        printf("Reached coarsest grid level.\n");
        return;
    } else {
        int coarse_grid_size = (grid_size - 1) / 2 + 1;
        for (int i = 0; i < grid_size; i++) {
            for (int j = 0; j < grid_size; j++) {
                if (i % 2 == 0 && j % 2 == 0) {
                    x_vals[current_level + 1][i / 2][j / 2] = x_vals[current_level][i][j];
                    y_vals[current_level + 1][i / 2][j / 2] = y_vals[current_level][i][j];
                }
            }
        }
        printf("Generated coarse grid (Level %d)...\n", current_level + 1);
        generateCoarseLevels(x_vals, y_vals, coarse_grid_size, levels, current_level + 1, output_flag);
    }
}

/*Restriction 적용 - up->down*/
void applyRestriction(double ***solution, double ***source, double ***temp,
                      int grid_size, int level) {
    int coarse_grid_size = (grid_size - 1) / 2 + 1;
    for (int i = 1; i < coarse_grid_size - 1; i++) {
        for (int j = 1; j < coarse_grid_size - 1; j++) {
            temp[level + 1][i][j] = (solution[level][2 * i - 1][2 * j - 1] +
                                     solution[level][2 * i - 1][2 * j + 1] +
                                     solution[level][2 * i + 1][2 * j - 1] +
                                     solution[level][2 * i + 1][2 * j + 1]) * 0.25;
        }
    }
    printf("Applied restriction to coarser grid level %d.\n", level + 1);
}

/*Prolongation 적용 - down->up*/
void applyProlongation(double ***solution, double ***temp, int grid_size, int level) {
    int fine_grid_size = (grid_size - 1) * 2 + 1;
    for (int i = 1; i < grid_size - 1; i++) {
        for (int j = 1; j < grid_size - 1; j++) {
            solution[level][2 * i][2 * j] += temp[level + 1][i][j];
            solution[level][2 * i + 1][2 * j] += 0.5 * (temp[level + 1][i][j] + temp[level + 1][i + 1][j]);
            solution[level][2 * i][2 * j + 1] += 0.5 * (temp[level + 1][i][j] + temp[level + 1][i][j + 1]);
            solution[level][2 * i + 1][2 * j + 1] += 0.25 * (temp[level + 1][i][j] + temp[level + 1][i + 1][j] +
                                                              temp[level + 1][i][j + 1] + temp[level + 1][i + 1][j + 1]);
        }
    }
    printf("Applied prolongation to finer grid level %d.\n", level);
}

/*Residual 계산*/
double calculateResidual(double **solution, double **source, double **residual, int grid_size) {
    double h2 = pow(1.0 / (grid_size - 1), 2.0);
    double residual_sum = 0.0;
    for (int i = 1; i < grid_size - 1; i++) {
        for (int j = 1; j < grid_size - 1; j++) {
            residual[i][j] = source[i][j] - (solution[i - 1][j] + solution[i + 1][j] +
                                             solution[i][j - 1] + solution[i][j + 1] -
                                             4.0 * solution[i][j]) / h2;
            residual_sum += residual[i][j] * residual[i][j];
        }
    }
    return sqrt(residual_sum);
}