/*
 * initialization.c:
 *  This file contains the functions to perform initialization operations, mostly for reading parameters.
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

# include "../header/initialization.h"
# include "../header/utility.h"

char line[BUFSIZE_L];

int initialization_real (int argc, char** argv)
{
    int i, j;
    int random;
    char configFileName[BUFSIZE_S];
    char PF_name[BUFSIZE_S];

    int c;
    int pos;

    FILE *PF     = NULL;
    FILE *config = NULL;

    weight_file = NULL;

    // SBX parameter settings
    pcross_real = 0.9;
    eta_c       = 20.0;

    // polynomial mutation parameter settings
    pmut_real = 1.0 / number_variable;
    eta_m     = 20.0;

    // differential evolution parameter settings
    CR = 0.5;
    F  = 0.5;
    K  = 0.5;

    // intrisic parameters used in MOEA/D variants
    neighbor_size = 20;
    function_type = ITCH;
    neighborhood_selection_probability = 0.9;
    maximumNumberOfReplacedSolutions = 2;

    // set the reference point for Hypervolume calculation
    ref_point = (double *) malloc (number_objective * sizeof(double));
    for(i = 0; i< number_objective; i++)
        ref_point[i] = 4.0;

    // calculate the number of points in the PF data
    sprintf (PF_name, "PF/%s.%dD.pf", problem_name, number_objective);
    PF = fopen (PF_name, "r");
    //print_error (PF == NULL, 2, "Fail to open PF: ", PF_name);
    if(PF != NULL){
        PF_size = 0;
        while (fgets (line, BUFSIZE_L, PF) != NULL)
            PF_size++;

        // read the PF data
        rewind (PF);
        PF_data = (double **) malloc (PF_size * sizeof(double *));
        for (i = 0; i < PF_size; i++)
            PF_data[i] = (double *) malloc (number_objective * sizeof(double));
        for (i = 0; i < PF_size; i++)
            for (j = 0; j < number_objective; j++)
                fscanf (PF, "%lf", &PF_data[i][j]);
    }
    // boundary settings
    variable_lowerbound = (double *) malloc (number_variable * sizeof(double));
    variable_upperbound = (double *) malloc (number_variable * sizeof(double));
    if (!strcmp(problem_name, "ZDT4"))
    {
        variable_lowerbound[0] = 0.0;
        variable_upperbound[0] = 1.0;
        for (i = 1; i < number_variable; i++)
        {
            variable_lowerbound[i] = -5.0;
            variable_upperbound[i] = 5.0;
        }
    }
    else
    {
        for (i = 0; i < number_variable; i++)
        {
            variable_lowerbound[i] = 0.0;
            variable_upperbound[i] = 1.0;
        }
    }

    // initialize a random seed
    srand ((unsigned) time (NULL));
    random = rand () % 1000;
    seed   = (float) random / 1000.0;
    print_error (seed <= 0.0 || seed >= 1.0, 1, "Entered seed value is wrong, seed value must be in (0,1)!");

    return 0;
}

/* Initialize the ideal point */
void initialize_idealpoint (void * pop)
{
    int i;

    population_real * ptr = (population_real*) pop;
    ideal_point = (double *) malloc (sizeof(double) * number_objective);    // free in moead_free

    for (i = 0; i < number_objective; i++)
        ideal_point[i] = INF;
        //ideal_point[i] = 0;

    for (i = 0 ;i < popsize ; i ++)
        update_ideal_point (&(ptr->ind[i]));

    return;
}

/* Initialize the nadir point */
void initialize_nadirpoint (void * pop)
{
    int i;

    population_real * ptr = (population_real*) pop;
    nadir_point = (double *) malloc (sizeof(double) * number_objective);    // free in moead_free

    for (i = 0; i < number_objective; i++)
        nadir_point[i] = -INF;

    for (i = 0 ;i < popsize ; i ++)
        update_nadir_point (&(ptr->ind[i]));

    return;
}
