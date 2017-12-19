/*
 * main.c:
 *  This is the main procedures of a general EMO algorithm (generational evolution model).
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Renzhi Chen, Ke Li
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

# include "header/rand.h"
# include "header/metaheuristics.h"

/* common paramters */
int run_index;
int run_index_begin;
int run_index_end;
int max_evaluation;              // maximum number of evaluations (stopping criterion)
int evaluation_count;            // evaluation counter
int popsize;                     // population size
int number_variable;             // number of variables
int number_objective;            // number of objectives
double* ideal_point;             // ideal point
double* nadir_point;             // nadir point
double* variable_lowerbound;     // variable lower bound
double* variable_upperbound;     // variable upper bound
char dummy[BUFSIZE_S];
char problem_name[BUFSIZE_S];
char algorithm_name[BUFSIZE_S];
char analyse_stream[BUFSIZE_L];
char problem_param_stream[BUFSIZE_L];
/* crossover and mutation */
double eta_c;                    // eta_c in SBX
double eta_m;                    // eta_m in polynomial mutation
double pcross_real;              // crossover rate for real encoded
double pmut_real;                // mutation rate for real encoded
double CR;                       // CR in DE
double F;                        // F in DE
double K;

/* performance metrics */
int PF_size;                 // size of the true Pareto-optimal Front
double **PF_data;            // true Pareto-optimal front data
double *ref_point;           // reference point for Hypervolume calculation

/* MOEA/D variants */
int neighbor_size;                           // neighborhood length
int number_weight;                           // number of weight vectors
char *weight_file;
int function_type;                           // type of the aggregation function
int maximumNumberOfReplacedSolutions;        // the maximum replacement number of a superior offspring
double neighborhood_selection_probability;   // probability to replace in the neighborhood
double **lambda;                             // weight vectors
int **neighborhood;                          // neighborhood structure
int *permutation;                            // subproblem index permutation
int *frequency;                              // subproblem usages counter arrary
double *utility;                             // subproblem utility array
struct int_vector *selected;
struct int_vector *candidate;

/* analysis platform */
int runtime_output;
int output_interval;
int analyse_list[BUFSIZE_S];
FILE *pythonplot;

char* global_base_dir;
double* true_nadir;
double* global_reference_point;
double region_of_interest = 2;

double* test_problem_intersection(char* test_problem, double *vector, int number_objectives);

int main(int argc, char *argv[])
{
    if(argc < 4) exit(-1);

    int test_idx =         (int) strtod(argv[1],NULL);
    int algorithm_idx =    (int) strtod(argv[2],NULL);
    int dimension_idx =    (int) strtod(argv[3],NULL);
    run_index =            (int) strtod(argv[4],NULL);

    population_real *parent_pop;
    population_real *offspring_pop;
    population_real *mixed_pop;

    randomize ();

    double non_dominance_threshold = 0.0001;
    double specificity = 0.02;
    double radius = 2;
    double epsilon = 0.001;

    int dimensions[] = {2, 3, 5, 8, 10, 15};
//    int dimensions_length = 3;

    int test_number_variables[] = {
            30, 30, 30, 10,
            7, 12, 12, 12};

    char *test_problems[] = {
            "ZDT1", "ZDT2", "ZDT4", "ZDT6",
            "DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4",
//            "WFG1", "WFG2", "WFG4", "WFG5", "WFG7", "WFG8", "WFG9"
    };

//    int test_problems_length = 10;

    char* algorithm_names[] = {
            "NSGA2", "MOEAD", "IBEA", "rNSGA2", "RNSGA2", "gNSGA2", "RMEAD", "PBEA"
    };

//    int algorithm_length = 8;

    int pop_sizes[] =       {100, 120, 120, 120, 100, 120};
    int moead_pop_sizes[] = {105, 120, 126, 120, 220, 120};

    number_objective = dimensions[dimension_idx];
    sprintf(problem_name, test_problems[test_idx]);
    sprintf(algorithm_name, algorithm_names[algorithm_idx]);

    if(algorithm_idx == 1)
        popsize = moead_pop_sizes[dimension_idx];
    else
        popsize = pop_sizes[dimension_idx];

    // Was 600 for 2, 3, 5
    max_evaluation = 1000 * popsize;

    double* reference_points[(number_objective+1) * 5];
    for(int i = 0; i < (number_objective+1) * 5; i++)
        reference_points[i] = malloc(sizeof(double) * number_objective);

    double weights[number_objective];
    for(int i = 0; i < number_objective; i++)
        weights[i] = 1/(double) number_objective;

    for(int r = 0; r < number_objective + 1; r++) {

        int ref_point_idx = r * 5;

        if (test_idx < 4 && r > 2) break;

        // All others
        if(r < number_objective)
            for(int i = 0; i < number_objective; i++) {
                reference_points[ref_point_idx][i] = 0;

                if(i == r)
                    reference_points[ref_point_idx][i] = 1;

                // ZDT3
                if(test_idx == 3 && r == 1 && i != 1)
                    reference_points[ref_point_idx][i] = -1;

                // DTLZ1
                if(test_idx == 5 && i == r)
                    reference_points[ref_point_idx][i] = 0.5;

                // WFG
    //            if(m == k && i > 10)
    //                reference_points[i][j][k][m] = pow(2, m);
            }

        if(r == number_objective) {
            for(int i = 0; i < number_objective; i++) {

                double sum = 0;
                for(int j = 0; j < number_objective; j++)
                    sum += reference_points[i][j];

                reference_points[ref_point_idx][i] = sum / number_objective;
            }

            reference_points[ref_point_idx]
                    = test_problem_intersection(test_problems[test_idx], reference_points[ref_point_idx], number_objective);
        }

        for(int i = 0; i < number_objective; i++) {
            reference_points[ref_point_idx + 1][i] = 0.5 * reference_points[ref_point_idx][i];  // Infeasible middle
            reference_points[ref_point_idx + 2][i] = 0.75 * reference_points[ref_point_idx][i]; // Infeasible near
            reference_points[ref_point_idx + 3][i] = 1.25 * reference_points[ref_point_idx][i]; // Feasible near
            reference_points[ref_point_idx + 4][i] = 1.5 * reference_points[ref_point_idx][i];  // Feasible middle
        }
    }

    global_base_dir = malloc(sizeof(char*) * BUFSIZE_L);

    // run experiments
    number_variable = test_number_variables[test_idx];

    initialization_real (argc,argv);

    number_objective = dimensions[dimension_idx];

    parent_pop    = (population_real *) malloc (sizeof(population_real));
    offspring_pop = (population_real *) malloc (sizeof(population_real));
    mixed_pop     = (population_real *) malloc (sizeof(population_real));

    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (offspring_pop, popsize);
    allocate_memory_pop (mixed_pop, 2 * popsize);

    true_nadir = malloc(sizeof(double) * number_objective);

    for(int i = 0; i < number_objective; i++)
        true_nadir[i] = 1;

    for(int i = 0; i < number_objective + 1; i++) {
        for(int j = 0; j < 5; j++) {

            if(algorithm_idx < 3 && (i + j) > 0 )
                continue;

            global_reference_point = reference_points[(i * 5) + j];

            char* distance[] = {"Perfect", "InMid", "InNear", "FeaNear", "FeaMid"};
            char* placement = i == number_objective ? "Cent" : "Obj";

            if(i == number_objective)
                sprintf(global_base_dir, "./out/%s/%i/%s/%s %s/%i", test_problems[test_idx], dimensions[dimension_idx], algorithm_names[algorithm_idx], distance[j], placement, run_index + 1);
            else
                sprintf(global_base_dir, "./out/%s/%i/%s/%s %s %i/%i", test_problems[test_idx], dimensions[dimension_idx], algorithm_names[algorithm_idx], distance[j], placement, i + 1, run_index + 1);

            _mkdir(global_base_dir);

            switch (algorithm_idx) {
                case 0 : NSGA2(parent_pop, offspring_pop, mixed_pop); break;
                case 1 : MOEAD(parent_pop, offspring_pop, mixed_pop); break;
                case 2 : IBEA(parent_pop, offspring_pop, mixed_pop); break;
                case 3 : rNSGA2(parent_pop, offspring_pop, mixed_pop, reference_points[i], weights, non_dominance_threshold); break;
                case 4 : RNSGA2(parent_pop, offspring_pop, mixed_pop, reference_points[i], weights, epsilon); break;
                case 5 : gNSGA2(parent_pop, offspring_pop, mixed_pop, reference_points[i]); break;
                case 6 : RMEAD2(parent_pop, offspring_pop, mixed_pop, reference_points[i], radius); break;
                case 7 : PBEA(parent_pop, offspring_pop, mixed_pop, reference_points[i], weights, specificity); break;
                default: exit(-1);
            }

            char population_file[BUFSIZE_L];
            sprintf(population_file, "%s/pop_final.out", global_base_dir);

            print_objective(population_file, parent_pop);
        }
    }

    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (offspring_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2 * popsize);

    free(parent_pop);
    free(offspring_pop);
    free(mixed_pop);

    free(true_nadir);

    printf("\n");
    printf ("TEST COMPLETE \n");

    // free memory
    if (number_variable != 0)
    {
        free (variable_lowerbound);
        free (variable_upperbound);
    }

    for (int i = 0; i < PF_size; i++)
        free (PF_data[i]);
    free (PF_data);

    return 0;
}

double* test_problem_intersection(char* test_problem, double *vector, int number_objectives) {

    double* intersection = malloc(sizeof(double) * number_objectives);

    if(strcmp(test_problem, "ZDT1") == 0 || strcmp(test_problem, "ZDT4") == 0) {
        double m = vector[1] / vector[0];
        double c = (-1 + sqrt(1 + (4 * m))) / (2 * m);

        intersection[0] = pow(c, 2);
        intersection[1] = 1 - c;

        return intersection;
    }

    if(strcmp(test_problem, "ZDT2") == 0 || strcmp(test_problem, "ZDT6") == 0) {
        double m = vector[1] / vector[0];

        intersection[0] = (-m + sqrt(m*m + 4)) * 0.5;
        intersection[1] = 1 - intersection[0]*intersection[0];

        return intersection;
    }

    if(strcmp(test_problem, "DTLZ1") == 0) {

        double sum = 0;
        for(int i = 0; i < number_objectives; i++)
            sum += vector[i];

        for(int i = 0; i < number_objectives; i++)
            intersection[i] = (0.5 * vector[i]) / sum;

        return intersection;
    }

    if(strcmp(test_problem, "DTLZ2") == 0 || strcmp(test_problem, "DTLZ3") == 0 || strcmp(test_problem, "DTLZ4") == 0) {

        double squared_sum = 0;
        for(int i = 0; i < number_objectives; i++)
            squared_sum += vector[i] * vector[i];

        for(int i = 0; i < number_objectives; i++)
            intersection[i] = vector[i] / sqrt(squared_sum);
    }

    if(strcmp(test_problem, "ZDT3") == 0 || strcmp(test_problem, "DTLZ7") == 0) {
        // Point in dominated region
        exit(1);
    }


    return vector;
}

