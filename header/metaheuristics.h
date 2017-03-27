/*
 * metaheuristics.h:
 *  This is the header file for metaheuristics.
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

#ifndef SAMARITAN_METAHEURISTICS_H
#define SAMARITAN_METAHEURISTICS_H

<<<<<<< HEAD
void NSGA2(population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
void MOEAD(population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
void MOEAD_DRA(population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
void MOEAD_STM(population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
void MOEAD_STM_DRA(population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
void SMSEMOA(population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
=======
/* external variables used by MOEA/D-STM variants */
extern int *idx ;
extern int *next;
extern int *statusWoman;
extern double *nicheCount;
extern double **distMatrix;
extern double **fitnessMatrix;
extern struct double_with_index **solMatrix;
extern struct double_with_index **subpMatrix;

void NSGA2 (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
void MOEAD (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
void MOEAD_DRA (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
void MOEAD_STM (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
void MOEAD_STM_DRA (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop);
>>>>>>> 785547837a435a72717cea48ad8d584f4d0ac66f

#endif //SAMARITAN_METAHEURISTICS_H
