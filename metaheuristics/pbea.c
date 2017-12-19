/*
 * pbea.c:
 *  This file contains the main procedures of PBEA. It is based on the following reference
 *
 *  L. Thiele, K. Miettinen, P. J. Korhonen, and J. Molina. A preference-based evolutionary algorithm for
 *  multi-objective optimization. Evolutionary computation, 17(3):411–436, 2009.
 *
 * Authors:
 *  Joe Billingsley <jb931@exeter.ac.uk>
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

# include "../header/metaheuristics.h"

void PBEA (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop, double* reference_point, double* weights, double specificity)
{
    int i;
    int generation;
    int *flag;
    double **fitcomp;

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop);
    evaluate_population (parent_pop);

    flag    = (int *) malloc (2 * popsize * sizeof(int));
    fitcomp = (double **) malloc (2 * popsize * sizeof(double*));
    for (i = 0; i < 2 * popsize; i++)
        fitcomp[i] = (double *) malloc (2 * popsize * sizeof(double*));

    // track the current evolutionary progress, including population and metrics
    pref_track_evolution (parent_pop, generation);
    while (evaluation_count < max_evaluation)
    {
        generation++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_real (parent_pop, offspring_pop);
        mutation_real (offspring_pop);
        evaluate_population (offspring_pop);

        // environmental selection
        merge (parent_pop, offspring_pop, mixed_pop);

        pbea_selection (mixed_pop, parent_pop, flag, fitcomp, weights, reference_point, specificity);

        // track the current evolutionary progress, including population and metrics
        pref_track_evolution (parent_pop, generation);
    }

    // release memory
    for (i = 0; i < 2 * popsize; i++)
        free (fitcomp[i]);
    free (fitcomp);
    free (flag);

    return;
}
