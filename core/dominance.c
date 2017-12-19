/*
 * dominance.c:
 *  This file contains the functions to perform the dominance check.
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

#include "../header/global.h"
#include "../header/population.h"
#include "../header/memory.h"

int check_dominance (individual_real *a, individual_real *b)
{
    int i;
    int flag1;
    int flag2;

    flag1 = flag2 = 0;
    for (i = 0; i < number_objective; i++)
    {
        if (a->obj[i] < b->obj[i])
            flag1 = 1;
        else
        {
            if (a->obj[i] > b->obj[i])
                flag2 = 1;
        }
    }
    if (flag1 == 1 && flag2 == 0)
        return (1);
    else
    {
        if (flag1 == 0 && flag2 == 1)
            return (-1);
        else
            return (0);
    }
}

int check_g_dominance_flag(individual_real ind_a, const double* reference_point);
int check_g_dominance(individual_real ind_a, individual_real ind_b, double* reference_point) {

    int flag_a = check_g_dominance_flag(ind_a, reference_point);
    int flag_b = check_g_dominance_flag(ind_b, reference_point);

    if(flag_a > flag_b)
        return 1;
    if(flag_a < flag_b)
        return -1;

    return check_dominance(&ind_a, &ind_b);
}

int check_g_dominance_flag(individual_real ind_a, const double* reference_point) {
    int a = 0, r = 0;

    for(int k = 0; k < number_objective; k++)
        if(ind_a.obj[k] <= reference_point[k])
            a++;
        else
            r++;

    return ((a || r) && !(a && r));
}

list** nondominated_sort_idxs(population_real* pop, int popsize) {

    int domination_counter[popsize];

    list* dominates[popsize];
    for(int i = 0; i < popsize; i++) {
        dominates[i] = NULL;
        domination_counter[i] = 0;
    }

    int* fronts_size = malloc(sizeof(int) * popsize);
    list** fronts = malloc(sizeof(list*) * popsize);

    for(int i = 0; i < popsize; i++) {
        fronts[i] = NULL;
        fronts_size[i] = 0;
    }

    for(int i = 0; i < popsize; i++) {

        for(int j = 0; j < popsize; j++) {
            int dominance = check_dominance(&pop->ind[i], &pop->ind[j]);

            if(dominance == 1) { // i dominates j
                list* dominated_ind = malloc(sizeof(list));
                dominated_ind->index = j;
                dominated_ind->child = NULL;

                if(dominates[i] == NULL)
                    dominates[i] = dominated_ind;
                else
                    append(dominated_ind, dominates[i]);

            } else if (dominance == -1) // j dominates i
                domination_counter[i] += 1;
        }

        if(domination_counter[i] == 0) {
            list* ind = malloc(sizeof(list));
            ind->index = i;
            ind->child = NULL;

            fronts_size[0] += 1;

            if(fronts[0] == NULL)
                fronts[0] = ind;
            else
                append(ind, fronts[0]);
        }
    }

    int curr_front = 0;
    int next_front = 1;

    while(fronts_size[curr_front] > 0) {

        for(int i = 0; i < fronts_size[curr_front]; i++) {

            list* dominated_inds = dominates[get_item(fronts[curr_front], i)->index];

            for(int j = 0; j < length(dominated_inds); j++) {

                list* dominated_ind = get_item(dominated_inds, j);

                domination_counter[dominated_ind->index] -= 1;

                if(domination_counter[dominated_ind->index] == 0) {

                    if(fronts[next_front] == NULL)
                        fronts[next_front] = dominated_ind;
                    else
                        append(dominated_ind, fronts[next_front]);

                    fronts_size[next_front] += 1;
                }

            }

        }

        curr_front++;
        next_front++;
    }

    free(fronts_size);

    return fronts;
}

population_real* get_nondominated(population_real *pop, int popsize, int *new_size) {

    int is_dominated[popsize];

    for(int i = 0; i < popsize; i++)
        is_dominated[i] = 0;

    int non_dominated_count = popsize;

    for(int i = 0; i < popsize; i++) {
        if(is_dominated[i]) continue;

        for(int j = 0; j < popsize; j++){
            if(is_dominated[j]) break;

            int dominates = check_dominance(&pop->ind[i], &pop->ind[j]);

            if(dominates) {
                is_dominated[j] = 1;

                non_dominated_count--;

                break;
            }
        }
    }

    *new_size = non_dominated_count;

    population_real *processed_pop = malloc(sizeof(population_real) * non_dominated_count);
    allocate_memory_pop(processed_pop, non_dominated_count);

    int cnt = 0;
    for(int i = 0; i < popsize; i++)
        if(!is_dominated[i])
            copy_ind(&(pop->ind[i]), &(processed_pop->ind[cnt++]));

    return processed_pop;
}