/*
 * r_metric.c:
 *  This file contains the functions to calculate the r-metric modified versions of GD, IGD and HV.
 *
 * Authors:
 *  Joe Billingsley <jb931@exeter.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
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

# include "../header/indicators.h"
#include "../header/dominance.h"
#include "../header/memory.h"

population_real* restrict_to_range(population_real *pop, int popsize, double *centroid, double range, int *new_size);

double calculate_r_mod(population_real *pop, double (*metric)(void *), const double *reference_point,
                       const double *weights, double roi) {

    double worst_point[number_objective];

    for(int i = 0; i < number_objective; i++)
        worst_point[i] = reference_point[i] + (2.0 * weights[i]);

    // Process population to remove dominated solutions
    int new_size = 0;

    population_real *nondominated_pop = get_nondominated(pop, popsize, &new_size);

    // Find the centroid of the processed pop
    double centroid[number_objective];

    for(int i = 0; i < number_objective; i++) {
        centroid[i] = 0;

        for(int j = 0; j < new_size; j++)
            centroid[i] += pop->ind[j].obj[i];

        centroid[i] = centroid[i] / new_size;
    }

    // Find the closest point to the centre
    int min_idx = 0;
    double min_dist = INF;
    for(int i = 0; i < popsize; i++) {
        double dist = euclidian_distance(pop->ind[i].obj, centroid, number_objective);

        if(dist < min_dist) {
            min_dist = dist;
            min_idx = i;
        }
    }

    individual_real centroid_ind = pop->ind[min_idx];

    // Only consider the points inside the region of interest
    population_real *reduced_pop = restrict_to_range(nondominated_pop, new_size, centroid_ind.obj, roi / 2, &new_size);
    free(nondominated_pop);

    int k = 0;
    double max_k = -INF;
    for(int i = 0; i < number_objective; i++) {
        double k_i = (centroid_ind.obj[i] - reference_point[i]) / (worst_point[i] - reference_point[i]);

        if(k_i > max_k) {
            max_k = k_i;
            k = i;
        }
    }

    double dir_vector[number_objective], shift[number_objective];
    for(int i = 0; i < number_objective; i++) {

        double grad = (centroid_ind.obj[k] - reference_point[k]) / (worst_point[k] - reference_point[k]);
        dir_vector[i] = reference_point[i] + grad * (worst_point[i] - reference_point[i]);

        shift[i] = dir_vector[i] - centroid_ind.obj[i];
    }

    for(int i = 0; i < new_size; i++)
        for(int j = 0; j < number_objective; j++)
            reduced_pop->ind[i].obj[j] += shift[j];

    int temp_popsize = popsize;
    popsize = new_size;

    double output = metric(reduced_pop);

    popsize = temp_popsize;

    return output;
}

population_real* restrict_to_range(population_real *pop, int popsize, double *centroid, double range, int *new_size) {

    int outside_range[popsize];

    for(int i = 0; i < popsize; i++)
        outside_range[i] = 0;

    int in_range_count = popsize;

    for(int i = 0; i < popsize; i++) {

        int in_range = euclidian_distance(pop->ind[i].obj, centroid, number_objective) < range;

        if(in_range) {
            outside_range[i] = 1;
            in_range_count--;

            break;
        }
    }

    *new_size = in_range_count;

    population_real *restricted_pop = malloc(sizeof(population_real) * in_range_count);
    allocate_memory_pop(restricted_pop, in_range_count);

    int cnt = 0;
    for(int i = 0; i < popsize; i++) {
        if(!outside_range[i])
            copy_ind(&(pop->ind[i]), &(restricted_pop->ind[cnt++]));
    }

    return restricted_pop;
}