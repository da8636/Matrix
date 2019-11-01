#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <string.h>
#include <math.h>


#include "matrix.h"

//static int g_seed = 0;

//static ssize_t g_width = 0;
//static ssize_t g_height = 0;
static ssize_t g_side = 0;
static ssize_t g_elements = 0;

static ssize_t g_nthreads = 1;


typedef struct {
  const size_t id;
  const float * matrix;
  float * value;
  float * thread_freq_results;
}arg_freq;

typedef struct {
  const size_t id;
  const float * matrix;
  float * thread_sum_results;
}arg_sum;

typedef struct {
  const size_t id;
  //const float * matrix;
  float * value;
  float * result_matrix;
}arg_uniform;

typedef struct {
  const size_t id;
  const float * matrix;
  float resultvalue;
  //float * result_matrix;
}arg_minmax;

typedef struct {
  const size_t id;
  //const float * matrix;
  float * value;
  float * result_matrix;
}arg_scalarmult;

typedef struct {
	const size_t id;
	float* result;
	const float* a;
	const float* b;
} arg_mMult;

typedef struct {
  const size_t id;
  //const float * matrix;
  float value;
  float step;
  float * result_matrix;
}arg_sequence;

typedef struct {
	float min;
	float max;
	float sum;
	//if i get time float * sorted_matrix
}arg_cache;
////////////////////////////////
///     UTILITY FUNCTIONS    ///
////////////////////////////////

/**
 * Returns pseudorandom number determined by the seed.
 */
// int fast_rand(void) {
//
// 	g_seed = (214013 * g_seed + 2531011);
// 	return (g_seed >> 16) & 0x7FFF;
//
// }

/**
 * Sets the seed used when generating pseudorandom numbers.
 */
// void set_seed(int seed) {
//
// 	g_seed = seed;
// }

/**
 * Sets the number of threads available.
 */
void set_nthreads(ssize_t count) {

	g_nthreads = count;
}

/**
 * Sets the dimensions of the matrix.
 */
void set_dimensions(ssize_t order) {


	g_side = order;

	g_elements = g_side * g_side;
	if (g_nthreads > g_elements){
		g_nthreads = g_elements;
	}
}

/**
 * Displays given matrix.
 */
void display(const float* matrix) {

	for (ssize_t y = 0; y < g_side; y++) {
		for (ssize_t x = 0; x < g_side; x++) {
			if (x > 0) printf(" ");
			printf("%.2f", matrix[y * g_side + x]);
		}

		printf("\n");
	}
}

/**
 * Displays given matrix row.
 */
void display_row(const float* matrix, ssize_t row) {
	for (ssize_t x = 0; x < g_side; x++) {
		if (x > 0) printf(" ");
		printf("%.2f", matrix[row * g_side + x]);
	}

	printf("\n");
}

/**
 * Displays given matrix column.
 */
void display_column(const float* matrix, ssize_t column) {

	for (ssize_t i = 0; i < g_side; i++) {
		printf("%.2f\n", matrix[i * g_side + column]);
	}
}

/**
 * Displays the value stored at the given element index.
 */
void display_element(const float* matrix, ssize_t row, ssize_t column) {
	printf("%.2f\n", matrix[row * g_side + column]);
}

////////////////////////////////
///   MATRIX INITALISATIONS  ///
////////////////////////////////

/**
 * Returns new matrix with all elements set to zero.
 */
float* new_matrix(void) {

	return calloc(1, sizeof(float)*g_elements + sizeof(arg_cache));
}

/**
 * Returns new identity matrix.
 */
float* identity_matrix(void) {

	float* result = new_matrix();

	int step = g_side+1;
	for (int i = 0; i < g_elements; i += step){
		result[i] = 1.00;
	}

	arg_cache* m = (arg_cache*)(result+g_elements);
	m->min = 0;
	m->max = 1;
	m->sum = g_side;

	/*
		TODO

		1 0
		0 1
	*/

	return result;
}

/**
 * Returns new matrix with elements generated at random using given seed.
 */
float* random_matrix(int seed) {

	float* result = new_matrix();

	for (ssize_t i = 0; i < g_elements; i++) {
		seed = (214013 * seed + 2531011);
		result[i] = (seed >> 16) & 0x7FFF;

	}

	return result;
}

void * worker_uniform(void * arg) {

	arg_uniform* argument = (arg_uniform *) arg;

	const size_t start = argument->id * (g_elements/g_nthreads);
	const size_t end = argument->id == g_nthreads - 1 ? g_elements : (argument->id + 1) * (g_elements/g_nthreads);

	// target value
	float value = * (argument->value);

	float* m = argument->result_matrix;
	for (size_t i = start; i < end; i++) {

		m[i] = value;

	}
	m[argument->id] = *m;


	//argument->thread_uniform_results[argument->id] = count;
	return NULL;
}



/**
 * Returns new matrix with all elements set to given value.
 */
float* uniform_matrix(float value) {

	float* result = new_matrix();

	if (g_elements > 500){

		pthread_t thread_ids[g_nthreads];

		float * thread_uniform_results = calloc(g_nthreads, sizeof(float));

		arg_uniform * args = malloc(sizeof(arg_uniform) * g_nthreads);

		for (size_t i = 0; i < g_nthreads; i++) {
	    	args[i] = (arg_uniform) {
	        	.result_matrix = result,
				.value = &value,
	        	.id = i
	    	};
			pthread_create(thread_ids + i, NULL, worker_uniform, args + i);

		}

		for (size_t i = 0; i < g_nthreads; ++i) {
	        // wait for threads to finish
	        pthread_join(thread_ids[i], NULL);
	    }

	    for (int i = 0; i < g_nthreads; ++i) {
			args->result_matrix[i] = value;
	    	//count += thread_uniform_results[i];
	    }

		//args->thread_freq_results[args->id] = count;

		free(args);
	    free(thread_uniform_results);

	}
	else{

		for (int i = 0; i < g_elements; i++){
			result[i] = value;
		}
	}
	//memset(result, value, sizeof(float)*g_elements);
	arg_cache* m = (arg_cache*)(result+g_elements);
	m->min = value;
	m->max = value;
	m->sum = value*g_elements;
	/*
		TODO

		     1 1
		1 => 1 1
	*/

	return result;
}



void * worker_sequence(void * arg) {

	arg_sequence* argument = (arg_sequence *) arg;

	const size_t start = argument->id * (g_elements/g_nthreads);
	const size_t end = argument->id == g_nthreads - 1 ? g_elements : (argument->id + 1) * (g_elements/g_nthreads);

	// target value
	float step = (argument->step);
	float value = (argument->value + step*start);

	float* m = argument->result_matrix;
	for (size_t i = start; i < end; i++) {
		m[i] = value + step*i;

	}
	m[argument->id] = *m;


	//argument->thread_uniform_results[argument->id] = count;
	return NULL;
}


/**
 * Returns new matrix with elements in sequence from given start and step
 */
float* sequence_matrix(float start, float step) {

	float* result = new_matrix();

	if (g_side > 22){
		pthread_t thread_ids[g_nthreads];

		//float * thread_uniform_results = calloc(g_nthreads, sizeof(float));

		arg_sequence * args = malloc(sizeof(arg_sequence) * g_nthreads);

		for (size_t i = 0; i < g_nthreads; i++) {
			args[i] = (arg_sequence) {
				.result_matrix = result,
				.value = start,
				.step = step,
				.id = i
			};
			pthread_create(thread_ids + i, NULL, worker_sequence, args + i);

		}

		for (size_t i = 0; i < g_nthreads; ++i) {
			// wait for threads to finish
			pthread_join(thread_ids[i], NULL);
		}

		free(args);
		//free(thread_uniform_results);
	}
	else{
		for (int i = 0; i < g_elements; i++){
			result[i] = start;
			start += step;
		}
	}
	arg_cache* m = (arg_cache*)(result+g_elements);


	if (step >= 0){
		m->min= result[0];
		m->max= result[g_elements-1];
	}
	else{
		m->min = result[g_elements-1];
		m->max = result[0];
	}
	if (g_elements%2 == 0){
		m->sum = (result[0]+result[g_elements-1])*(g_elements/2);
	}
	else {
		m->sum = (result[0]+result[g_elements-1])*(g_elements/2) + result[(g_elements/2)];
	}

	/*
		TODO

		       1 2
		1 1 => 3 4
	*/

	return result;
}

////////////////////////////////
///     MATRIX OPERATIONS    ///
////////////////////////////////

/**
 * Returns new matrix with elements cloned from given matrix.
 */
//TODO memcopy
float* cloned(const float* matrix) {

	float* result = new_matrix();
	memcpy(result, matrix, sizeof(float)*g_elements + sizeof(arg_cache));

	return result;
}

int comparsionfunction (const void * x1, const void * x2)
{
   return ( *(float*)(x1) - *(float*)(x2) );
}

/**
 * Returns new matrix with elements in ascending order.
 */
float* sorted(const float* matrix) {

	float* result = cloned(matrix);

	qsort(result, g_elements, sizeof(float), comparsionfunction);

	/*
		TODO

		3 4    1 2
		2 1 => 3 4

	*/

	return result;
}

/**
 * Returns new matrix with elements rotated 90 degrees clockwise.
 */
float* rotated(const float* matrix) {

	float* result = new_matrix();

	for (int i = 0; i < g_side; i++){
		for (int j = 0; j < g_side; j++){
			result[g_side*i+j] = matrix[g_side * (g_side - 1 - j)+i];

		}
	}

	/*
		TODO

		1 2    3 1
		3 4 => 4 2
	*/

	return result;
}

/**
 * Returns new matrix with elements ordered in reverse.
 */
float* reversed(const float* matrix) {

	float* result = new_matrix();
	for (int i = 0; i < (g_elements); i++){
		result[i] = matrix[g_elements -1 -i];
	}

	/*
		TODO

		1 2    4 3
		3 4 => 2 1
	*/

	return result;
}

/**
 * Returns new transposed matrix.
 */
float* transposed(const float* matrix) {

	float* result = new_matrix();

	for (int i = 0; i < g_side; i++){
		for (int j = 0; j < g_side; j++){
			result[g_side*i+j] = matrix[g_side*j+i];
			//result[g_width*i+j] = matrix[g_width*j+i];
		}
	}

	/*
		TODO

		1 2    1 3
		3 4 => 2 4
	*/

	return result;
}

void * worker_scalaradd(void * arg) {
	//float count = 0;
	arg_scalarmult* argument = (arg_scalarmult *) arg;

	const size_t start = argument->id * (g_elements/g_nthreads);
	const size_t end = argument->id == g_nthreads - 1 ? g_elements : (argument->id + 1) * (g_elements/g_nthreads);

	// target value
	float value = * (argument->value);

	float* m = argument->result_matrix;
	for (size_t i = start; i < end; i++) {

		//if (argument->matrix[i] == value){
			//argument->thread_freq_results += 1;
		m[i] += value;
		//}
	}
	//argument->thread_uniform_results[argument->id] = count;
	return NULL;
}


/**
 * Returns new matrix with scalar added to each element.
 */
float* scalar_add(const float* matrix, float value) {

	float* result = cloned(matrix);

	if (true){

		pthread_t thread_ids[g_nthreads];

		//float * thread_scalaradd_results = calloc(g_nthreads, sizeof(float));

		arg_scalarmult * args = malloc(sizeof(arg_scalarmult) * g_nthreads);

		for (size_t i = 0; i < g_nthreads; i++) {
	    	args[i] = (arg_scalarmult) {
	        	.result_matrix = result,
				.value = &value,
	        	.id = i
	    	};
			pthread_create(thread_ids + i, NULL, worker_scalaradd, args + i);

		}

		for (size_t i = 0; i < g_nthreads; ++i) {
	        // wait for threads to finish
	        pthread_join(thread_ids[i], NULL);
	    }


		free(args);


	} else{
		for (int i = 0; i < g_elements; i++){
			result[i] += value;

		}
	}
	//parallel_args* args = parallelFunction(matrix, parallel_scalar_add, scalar, result);
	return result;




	/*
		TODO

		1 0        2 1
		0 1 + 1 => 1 2

		1 2        5 6
		3 4 + 4 => 7 8
	*/
	return result;
}


void * worker_scalarmult(void * arg) {
	//float count = 0;
	arg_scalarmult* argument = (arg_scalarmult *) arg;

	const size_t start = argument->id * (g_elements/g_nthreads);
	const size_t end = argument->id == g_nthreads - 1 ? g_elements : (argument->id + 1) * (g_elements/g_nthreads);

	// target value
	float value = * (argument->value);

	float* m = argument->result_matrix;
	for (size_t i = start; i < end; i++) {

		//if (argument->matrix[i] == value){
			//argument->thread_freq_results += 1;
		m[i] *= value;
		//}
	}
	//argument->result_matrix[argument->id] = *argument->result_matrix;

	//argument->thread_uniform_results[argument->id] = count;
	return NULL;
}


//TO DO
/**
 * Returns new matrix with scalar multiplied to each element.
 */
float* scalar_mul(const float* matrix, float value) {


	float* result = cloned(matrix);

	if (true){

		pthread_t thread_ids[g_nthreads];

		//float * thread_scalaradd_results = calloc(g_nthreads, sizeof(float));

		arg_scalarmult * args = malloc(sizeof(arg_scalarmult) * g_nthreads);

		for (size_t i = 0; i < g_nthreads; i++) {
	    	args[i] = (arg_scalarmult) {
	        	.result_matrix = result,
				.value = &value,
	        	.id = i
	    	};
			pthread_create(thread_ids + i, NULL, worker_scalarmult, args + i);

		}

		for (size_t i = 0; i < g_nthreads; ++i) {
	        // wait for threads to finish
	        pthread_join(thread_ids[i], NULL);
	    }

	    // for (int i = 0; i < g_nthreads; ++i) {
		// 	args->result_matrix[i] = value;
	    // 	//count += thread_uniform_results[i];
	    // }

		//args->thread_freq_results[args->id] = count;

		free(args);


	} else{

		for (int i = 0; i < g_elements; i++){
			result[i] *= value;

		}
	}

	arg_cache* m = (arg_cache*)(result+g_elements);
	if (value >= 0){
		m->min = m->min*value;
		m->max = m->max*value;
		m->sum = m->sum*value;
	}
	else{
		m->min = m->max*value;
		m->max = m->min*value;
		m->sum = m->sum*value;
	}


	/*
		TODO

		1 0        2 0
		0 1 x 2 => 0 2

		1 2        2 4
		3 4 x 2 => 6 8
	*/

	return result;
}

/**
 * Returns new matrix that is the result of
 * adding the two given matrices together.
 */
float* matrix_add(const float* matrix_a, const float* matrix_b) {

	float* result = new_matrix();

	for (int i = 0; i < g_elements; i++){
		result[i] = matrix_a[i]+matrix_b[i];

	}


	/*
		TODO

		1 0   0 1    1 1
		0 1 + 1 0 => 1 1

		1 2   4 4    5 6
		3 4 + 4 4 => 7 8
	*/

	return result;
}


void* worker_mMult(void* args) {

	arg_mMult* wargs = (arg_mMult*) args;


	const size_t start = wargs->id * (g_side/g_nthreads);
	const size_t end = wargs->id == g_nthreads - 1 ? g_side : (wargs->id + 1) * (g_side/g_nthreads);

	float* m = wargs->result;
	const float* m1 = wargs->a;
	const float* m2 = wargs->b;

	for (size_t i = start; i < end; i++) {
		for (size_t k = 0; k < g_side; k++) {
			for (size_t j = 0; j < g_side; j++) {

				m[i * g_side + j] += m1[(i) * g_side + (k)] * m2[(k) * g_side + (j)];
			}
		}
	}

	return NULL;
}



/**
 * Returns new matrix that is the result of
 * multiplying the two matrices together.
 */
float* matrix_mul(const float* matrix_a, const float* matrix_b) {

	float* result = new_matrix();
	if (true){
		arg_mMult args[g_nthreads];
		pthread_t thread_ids[g_nthreads];

		for (size_t i = 0; i < g_nthreads; i++) {
			args[i] = (arg_mMult) {
				.a = matrix_a,
				.b = matrix_b,
				.id = i,
				.result = result,
			};
		}

		// Launch threads
		for (size_t i = 0; i < g_nthreads; i++) {
			pthread_create(thread_ids + i, NULL, worker_mMult, args + i);
		}

		// Wait for threads to finish
		for (size_t i = 0; i < g_nthreads; i++) {
			pthread_join(thread_ids[i], NULL);
		}

		// return result;
	}
	else{
		for (int i = 0; i < g_side; i++){
			for (int j = 0; j < g_side; j++){

				float totalmult = 0;
				for (int l = 0; l < g_side; l++){
					totalmult += matrix_a[g_side*i+l]*matrix_b[g_side*l+j];
				}
				result[g_side*i+j] = totalmult;

			}
		}
	}






	/*
		TODO

		(1*1)+(2*0) (1*0)+(2*1)
		(3*1)+(4*0) (3*0)+(4*1)



		1 2   1 0    1 2
		3 4 x 0 1 => 3 4

		1 2   5 6    19 22
		3 4 x 7 8 => 43 50
	*/

	return result;
}

/**
 * Returns new matrix that is the result of
 * powering the given matrix to the exponent.
 */
float* matrix_pow(const float* matrix, int exponent) {

	float* result;

	if (g_elements == 1){
		result = identity_matrix();
		result[0] = (float) pow(matrix[0], exponent);
	}

	if (exponent == 0){
		result = identity_matrix();
	}

	// while (exponent != 0){
	//
	// }
	// return 0;


	else{
		result = identity_matrix();
		float* squareResult = cloned(matrix);

		int total = exponent;
		while (total > 0){

			if (total%2 == 1){
				float* temp = result;
				result = matrix_mul(squareResult, result);
				free(temp);
			}
			total = total/2;
			float* temp = squareResult;
			squareResult = matrix_mul(squareResult, squareResult);
			free(temp);
		}
		free(squareResult);
	}
	/*
		TODO

		1 2        1 0
		3 4 ^ 0 => 0 1

		1 2        1 2
		3 4 ^ 1 => 3 4

		1 2        199 290
		3 4 ^ 4 => 435 634
	*/

	return result;
}

/**
 * Returns new matrix that is the result of
 * convolving given matrix with a 3x3 kernel matrix.
 */
float* matrix_conv(const float* matrix, const float* kernel) {

	float* result = new_matrix();

	if (g_elements == 1){
		result[0]
			= matrix[0] * (kernel[0]
			+ kernel[1]
			+ kernel[2]
			+ kernel[3]
			+ kernel[4]
			+ kernel[5]
			+ kernel[6]
			+ kernel[7]
			+ kernel[8]);
		return result;
	}

	//TOP LEFT CORNER
	result[0]
		= matrix[0] * (kernel[0] + kernel[1] + kernel[3] + kernel[4])
		+ matrix[1] * (kernel[2] + kernel[5])
		+ matrix[g_side] * (kernel[6] + kernel[7])
		+ matrix[g_side + 1] * kernel[8];

	//TOP RIGHT CORNER

	result[g_side-1]
		= matrix[g_side - 1] * (kernel[1] + kernel[2] + kernel[4] + kernel[5])
		+ matrix[g_side - 2] * (kernel[0] + kernel[3])
		+ matrix[(g_side * 2) - 1] * (kernel[7] + kernel[8])
		+ matrix[(g_side * 2) - 2] * kernel[6];


	//BOTTOM LEFT CORNER
	result[g_side * (g_side-1)]
		= matrix[g_side * (g_side-1)] * (kernel[6] + kernel[7] + kernel[3] + kernel[4])
		+ matrix[(g_side * (g_side-1)) + 1] * (kernel[8] + kernel[5])
		+ matrix[(g_side * (g_side-2))] * (kernel[0] + kernel[1])
		+ matrix[(g_side * (g_side-2)) + 1] * kernel[2];

	//BOTTOM RIGHT CORNER
	result[(g_side * g_side) - 1]
		= matrix[(g_side * g_side) - 1] * (kernel[8] + kernel[7] + kernel[5] + kernel[4])
		+ matrix[(g_side * g_side) - 2] * (kernel[6] + kernel[3])
		+ matrix[(g_side * (g_side-1)) - 1] * (kernel[2] + kernel[1])
		+ matrix[(g_side * (g_side-1)) - 2] * kernel[0];


	//top edge minus corners
	for (int i = 1; i < g_side-1; i++){
		result[i]
			= matrix[i - 1] * (kernel[0] + kernel[3])
			+ matrix[i] * (kernel[1] + kernel[4])
			+ matrix[i + 1] * (kernel[2] + kernel[5])
			+ matrix[i + g_side - 1] * (kernel[6])
			+ matrix[i + g_side] * (kernel[7])
			+ matrix[i + g_side + 1] * (kernel[8]);


	}

	//bottom edge minus corners
	for (int i = g_elements - g_side + 1; i < g_elements - 1; i++){

		result[i]
			= matrix[i - 1] * (kernel[6] + kernel[3])
			+ matrix[i] * (kernel[7] + kernel[4])
			+ matrix[i + 1] * (kernel[8] + kernel[5])
			+ matrix[i - g_side - 1] * (kernel[0])
			+ matrix[i - g_side] * (kernel[1])
			+ matrix[i - g_side + 1] * (kernel[2]);
	}

	//left edge minus corners
	for (int i = g_side; i < g_elements - g_side; i += g_side){
		result[i]
			= matrix[i - g_side] * (kernel[0] + kernel[1])
			+ matrix[i] * (kernel[3] + kernel[4])
			+ matrix[i + g_side] * (kernel[6] + kernel[7])
			+ matrix[i - g_side + 1] * kernel[2]
			+ matrix[i + 1] * kernel[5]
			+ matrix[i + g_side + 1] * kernel[8];
	}

	//right edge minus corners
	for (int i = g_side + g_side - 1; i < g_elements - 1; i += g_side){
		result[i]
			= matrix[i - g_side] * (kernel[2] + kernel[1])
			+ matrix[i] * (kernel[5] + kernel[4])
			+ matrix[i + g_side] * (kernel[8] + kernel[7])
			+ matrix[i - g_side - 1] * kernel[0]
			+ matrix[i - 1] * kernel[3]
			+ matrix[i + g_side - 1] * kernel[6];
	}

	//inside minus the frame
	for (int i = 1; i < g_side-1; i++){
		for (int j = 1; j < g_side-1; j++){
			int setPos = i*g_side+j;
			result[setPos]
				= matrix[setPos - g_side - 1] * kernel[0]
				+ matrix[setPos - g_side] * kernel[1]
				+ matrix[setPos - g_side + 1] * kernel[2]
				+ matrix[setPos - 1] * kernel[3]
				+ matrix[setPos] * kernel[4]
				+ matrix[setPos + 1] * kernel[5]
				+ matrix[setPos + g_side - 1] * kernel[6]
				+ matrix[setPos + g_side] * kernel[7]
				+ matrix[setPos + g_side + 1] * kernel[8];
		}
	}


	/*
		TODO

		Convolution is the process in which the values of a matrix are
		computed according to the weighted sum of each value and it's
		neighbours, where the weights are given by the kernel matrix.
	*/

	return result;
}

////////////////////////////////
///       COMPUTATIONS       ///
////////////////////////////////

void * worker_sum(void * arg) {
	float count = 0;
	//reuse
	arg_sum* argument = (arg_sum *) arg;

	const size_t start = argument->id * (g_elements/g_nthreads);
	const size_t end = argument->id == g_nthreads - 1 ? g_elements : (argument->id + 1) * (g_elements/g_nthreads);

	const float* m = argument->matrix;
	for (size_t i = start; i < end; i++) {
		count += m[i];
	}
	//argument->thread_sum_results[argument->id] = count;



	// for (size_t i = start; i < end; i++) {
	// 	if (argument->matrix[i] == value){
    //                    // This line should be moved after this for loop is done
	// 		argument->thread_freq_results += 1;
	// 		count++;
	// 	}
	// }

	// arg->thread_freq_results[arg->id] = count;
	argument->thread_sum_results[argument->id] = count;
	return NULL;
}

/**
 * Returns the sum of all elements.
 */
float get_sum(const float* matrix) {
	arg_cache* m = (arg_cache*)(matrix+g_elements);
	if (m->sum != 0){
		return m->sum;
	}
	float count = 0;
	if (g_elements > 500){

		pthread_t thread_ids[g_nthreads];

		float * thread_sum_results = calloc(g_nthreads, sizeof(float));

		arg_sum * args = malloc(sizeof(arg_sum) * g_nthreads);

		for (size_t i = 0; i < g_nthreads; i++) {
	    	args[i] = (arg_sum) {
	        	.matrix = matrix,
	        	.thread_sum_results = thread_sum_results,
	        	.id = i
	    	};
			pthread_create(thread_ids + i, NULL, worker_sum, args + i);

		}

		for (size_t i = 0; i < g_nthreads; ++i) {
	        // wait for threads to finish
	        pthread_join(thread_ids[i], NULL);
	    }

	    for (int i = 0; i < g_nthreads; ++i) {
	    	count += thread_sum_results[i];
	    }

		//args->thread_freq_results[args->id] = count;

		free(args);
	    free(thread_sum_results);

	}

	else {
		for (int i = 0; i < g_elements; i++){
			count += matrix[i];
		}
	}
	return count;

	/*
		TODO

		2 1
		1 2 => 6

		1 1
		1 1 => 4
	*/
	//
	// return sum;
}


/**
 * Returns the trace of the matrix.
 */
float get_trace(const float* matrix) {

	//int step = g_width + 1;
	int step = g_side + 1;
	float total = 0;
	for (int i = 0; i < g_elements; i += step){
		total += matrix[i];
	}
	/*


		TODO
		1 0 0 1
		1 0
		0 1 => 2

		2 1
		1 2 => 4
	*/

	return total;
}

void * worker_min(void * arg) {

	//reuse
	arg_minmax* argument = (arg_minmax *) arg;

	const size_t start = argument->id * (g_elements/g_nthreads);
	const size_t end = argument->id == g_nthreads - 1 ? g_elements : (argument->id + 1) * (g_elements/g_nthreads);

	// target value
	const float* m = argument->matrix;
	float count = m[start];
	for (size_t i = start+1; i < end; i++) {

		if (m[i] < count){
			count = m[i];
		}

	}


	// float min = matrix[0];
	//
	// for (int i = 1; i < g_elements; i++){
	// 	if (matrix[i] < min){
	// 		min = matrix[i];
	// 	}
	// }
	argument->resultvalue = count;
	return NULL;
}






/**
 * Returns the smallest value in the matrix.
 */
float get_minimum(const float* matrix) {

	arg_cache* m = (arg_cache*)(matrix+g_elements);
	if (m->min != 0){
		return m->min;
	}

	float count = 0;
	if (true){

		pthread_t thread_ids[g_nthreads];

		//float * thread_min_results = calloc(g_nthreads, sizeof(float));

		arg_minmax * args = malloc(sizeof(arg_minmax) * g_nthreads);

		for (size_t i = 0; i < g_nthreads; i++) {
	    	args[i] = (arg_minmax) {
	        	.matrix = matrix,
	        	//.thread_minmax_results = thread_minmax_results,
				.resultvalue = 0,
	        	.id = i
	    	};
			pthread_create(thread_ids + i, NULL, worker_min, args + i);

		}
		for (size_t i = 0; i < g_nthreads; ++i) {
	        // wait for threads to finish
	        pthread_join(thread_ids[i], NULL);
	    }
		count = args[0].resultvalue;

	    for (int i = 1; i < g_nthreads; ++i) {
			if (args[i].resultvalue < count){
				count = args[i].resultvalue;
			}
	    	//count += thread_minmax_results[i];
	    }
		// args->resultvalue = count;
		//args->thread_freq_results[args->id] = count;

		free(args);


	}
	else{
		count = matrix[0];
		for (int i = 1; i < g_elements; i++){
			if (matrix[i] < count){
				count = matrix[i];
			}
		}
	}



	/*
		TODO

		1 2
		3 4 => 1

		4 3
		2 1 => 1
	*/

	return count;
}

void * worker_max(void * arg) {

	//reuse
	arg_minmax* argument = (arg_minmax *) arg;

	const size_t start = argument->id * (g_elements/g_nthreads);
	const size_t end = argument->id == g_nthreads - 1 ? g_elements : (argument->id + 1) * (g_elements/g_nthreads);

	// target value
	const float* m = argument->matrix;
	float count = argument->matrix[start];
	for (size_t i = start+1; i < end; i++) {

		if (m[i] > count){
			count = m[i];
		}

	}


	// float min = matrix[0];
	//
	// for (int i = 1; i < g_elements; i++){
	// 	if (matrix[i] < min){
	// 		min = matrix[i];
	// 	}
	// }
	argument->resultvalue = count;
	return NULL;
}



/**
 * Returns the largest value in the matrix.
 */
float get_maximum(const float* matrix) {
	arg_cache* m = (arg_cache*)(matrix+g_elements);
	if (m->max != 0){
		return m->max;
	}
	float count = 0;
	if (true){

		pthread_t thread_ids[g_nthreads];

		//float * thread_min_results = calloc(g_nthreads, sizeof(float));

		arg_minmax * args = malloc(sizeof(arg_minmax) * g_nthreads);

		for (size_t i = 0; i < g_nthreads; i++) {
	    	args[i] = (arg_minmax) {
	        	.matrix = matrix,
	        	//.thread_minmax_results = thread_minmax_results,
				.resultvalue = 0,
	        	.id = i
	    	};
			pthread_create(thread_ids + i, NULL, worker_max, args + i);

		}
		for (size_t i = 0; i < g_nthreads; ++i) {
	        // wait for threads to finish
	        pthread_join(thread_ids[i], NULL);
	    }
		count = args[0].resultvalue;

	    for (int i = 1; i < g_nthreads; ++i) {
			if (args[i].resultvalue > count){
				count = args[i].resultvalue;
				//printf("%.2f is greater than %.2f\n", args[i].resultvalue, count);
			}
	    	//count += thread_minmax_results[i];

	    }
		// args->resultvalue = count;
		//args->thread_freq_results[args->id] = count;

		free(args);


	}
	else{
		count = matrix[0];

		for (int i = 1; i < g_elements; i++){
			if (matrix[i] > count){
				count = matrix[i];
			}
		}
	}




	/*
		TODO

		1 2
		3 4 => 4

		4 3
		2 1 => 4
	*/

	return count;
}


float det(const float* matrix, int width){
	if (width == 2){
		return (matrix[3]*matrix[0] - matrix[2]*matrix[1]);
	}
	float sum = 0;
	float* submatrix = new_matrix();
	for (int i = 0; i < width; i++){
		for (int j = 1; j < width; j++){
			for (int k = 0; k < i; k++){
				submatrix[(width-1)*(j-1)+k] = matrix[width*j+k];
			}
			for (int k = i+1; k < width; k++){
				submatrix[(width-1)*(j-1)+k-1] = matrix[width*j+k];
			}
		}
		if (i % 2 == 1){
			sum -= matrix[i] * det(submatrix, width-1);
		}
		else{
			sum += matrix[i] * det(submatrix, width-1);
		}

	}
	free(submatrix);
	return sum;
}


/**
 * Returns the determinant of the matrix.
 */
float get_determinant(const float* matrix) {
	if (g_elements == 1){
		return matrix[0];
	}
	else{
		//return det(matrix, g_width);
		return det(matrix, g_side);
	}
}


void * worker_freq(void * arg) {
	float count = 0;
	arg_freq* argument = (arg_freq *) arg;

	const size_t start = argument->id * (g_elements/g_nthreads);
	const size_t end = argument->id == g_nthreads - 1 ? g_elements : (argument->id + 1) * (g_elements/g_nthreads);

	// target value
	float value = * (argument->value);

	const float* m = argument->matrix;
	for (size_t i = start; i < end; i++) {
		if (m[i] == value){
			//argument->thread_freq_results += 1;
			count++;
		}
	}
	//argument->thread_freq_results[argument->id] = count;



	// for (size_t i = start; i < end; i++) {
	// 	if (argument->matrix[i] == value){
    //                    // This line should be moved after this for loop is done
	// 		argument->thread_freq_results += 1;
	// 		count++;
	// 	}
	// }


	// arg->thread_freq_results[arg->id] = count;
	argument->thread_freq_results[argument->id] = count;
	return NULL;
}



/**
 * Returns the frequency of the given value in the matrix.
 */
ssize_t get_frequency(const float* matrix, float value){

	ssize_t count = 0;
	if (g_elements > 500){

		pthread_t thread_ids[g_nthreads];

		float * thread_freq_results = calloc(g_nthreads, sizeof(float));

		arg_freq * args = malloc(sizeof(arg_freq) * g_nthreads);

		for (size_t i = 0; i < g_nthreads; i++) {
	    	args[i] = (arg_freq) {
	        	.matrix = matrix,
	        	.value = &value,
	        	.thread_freq_results = thread_freq_results,
	        	.id = i
	    	};
			pthread_create(thread_ids + i, NULL, worker_freq, args + i);

		}

		for (size_t i = 0; i < g_nthreads; ++i) {
	        // wait for threads to finish
	        pthread_join(thread_ids[i], NULL);
	    }

	    for (int i = 0; i < g_nthreads; ++i) {
	    	count += thread_freq_results[i];
	    }

		//args->thread_freq_results[args->id] = count;

		free(args);
	    free(thread_freq_results);

	} else{
		for (int i = g_elements; --i >= 0;){
			if (matrix[i] == value){
				count++;
			}
		}
	}

	return count;


	/*
		TODO

		1 1
		1 1 :: 1 => 4

		1 0
		0 1 :: 2 => 0
	*/

	// return freq;
}
