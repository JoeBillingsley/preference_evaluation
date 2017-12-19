/*
 * g-NSGA2.c:
 *  This file implements the main procedures of g-NSGAII. It is based on the following reference:
 *
 *  J. Molina, L. V. Santana, A. G. Hernández-Dı́az, C. A. C. Coello, and R. Caballero. g-dominance: Reference point
 *  based dominance for multiobjective metaheuristics. European Journal of Operational Research, 197(2):685–692, 2009.
 *
 * Authors:
 *  Joe Billingsley <jb931@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Renzhi Chen, Ke Li, Joe Billingsley
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
#include "../header/global.h"

void gNSGA2(population_real *parent_pop, population_real *offspring_pop, population_real *mixed_pop,
            double *reference_point)
{
    int generation;

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop);
    evaluate_population (parent_pop);

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

        initialize_idealpoint(mixed_pop);
        initialize_nadirpoint(mixed_pop);

        fill_g_nondominated_sort (parent_pop, mixed_pop, reference_point);

        // track the current evolutionary progress, including population and metrics
        pref_track_evolution (parent_pop, generation);
    }
}

