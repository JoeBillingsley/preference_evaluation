/*
 * dominance.h:
 *  This is the header file for dominance check.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
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

# ifndef Samaritan_DOMINANCE_H
# define Samaritan_DOMINANCE_H

int check_dominance(individual_real *a, individual_real *b);
int check_g_dominance(individual_real ind_a, individual_real ind_b, double* reference_point);

list** nondominated_sort_idxs(population_real* pop, int popsize);
population_real* get_nondominated(population_real *pop, int popsize, int *new_size);

# endif // Samaritan_DOMINANCE_H
