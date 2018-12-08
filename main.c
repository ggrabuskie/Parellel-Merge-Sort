/*
Georden Grabuskie
Anette Ulrichsen
CPSC 474
Project 2
Based off mpi hello world code from Eclipse PTP.
 */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <stdlib.h>
#include <time.h>

void merge(int arr[], int l, int m, int r);
void mergeSort(int arr[], int l, int r);
void final_merge(int arr[], int p, int size);
void sort_extra(int arr[], int size, int p);

int main(int argc, char* argv[]){
	int  my_rank; /* rank of process */
	int  p;       /* number of processes */
	/* start up MPI */
	
	MPI_Init(&argc, &argv); //takes up to one argument which would be the number of random elements to generate and sort
	
	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	
	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &p); 
	
	int size = 64;
	if(argc > 1) //for running from CLI
	{
		size = atoi(argv[1]);
	}

	int subsize = (size / p); //number of values to send to each process


	int* main_array = NULL;
	int* sub_array = malloc(sizeof(int) * subsize); //every process gets a subarray

	int sort_start_time;
	if (my_rank == 0) //only the root process needs to generate values
	{
		printf("Processes: %d\nSize: %d\nSubsize: %d\nRemainder: %d\n", p, size, subsize, size%p);//display once

		main_array = malloc(sizeof(int) * size); //dynamic array for random ints
		srand(time(NULL));
		for (int i = 0; i < size; ++i)
		{
			main_array[i] = (rand() % 10); //keeping it to single digits
		}
		printf("Initial array: ");
		for (int i = 0; i < size; ++i)
		{
			printf("%d ", main_array[i]);
		}
		printf("\n");
		sort_start_time = MPI_Wtime();


	}


	MPI_Scatter(main_array, subsize, MPI_INT, sub_array, subsize, MPI_INT, 0, MPI_COMM_WORLD); //now each process has its own piece of the main array

	mergeSort(sub_array, 0, subsize-1); //every process sorts its own sub array

	if(my_rank == 0)
	{
		sort_extra(main_array, size, p); //The remaining values that didn't get scattered are sorted in place in the main array on root 0
	}

	MPI_Gather(sub_array, subsize, MPI_INT, main_array, subsize, MPI_INT, 0, MPI_COMM_WORLD); //the sorted sub arrays are gathered in pieces bak into the main array on process 0

	if(my_rank == 0)
	{

		final_merge(main_array, p, size);
		printf("Sorted array:  ");
		for (int i = 0; i < size; ++i)
		{
			printf("%d ", main_array[i]);
		}
		printf("\nDuration: %f seconds", (MPI_Wtime() - sort_start_time));
		free(main_array);
	}				//Don't forget to free memory.
	free(sub_array);


	/* shut down MPI */
	MPI_Finalize(); 
	
	
	return 0;

}

//sorts the extra values at the end of the array that are left over after scatter
void sort_extra(int arr[], int size, int p)
{
	int remainder = size % p;
	if(remainder)
	{
		mergeSort(arr, (size - remainder), (size - 1));
	}
}

//merges the p sorted subarrays in place
//Root process calls merge iteratively on the sorted blocks of the main array
void final_merge(int arr[], int p, int size)
{
	int subsize = size / p;
	int shift_size = (2 * subsize);
	int left, middle, right;
	while (shift_size <= (size))
	{
		left = 0;
		middle = (subsize -1);
		right = ((2 * subsize) - 1);

		while (right <= (size-1))
		{
			merge(arr, left, middle, right);
			left += shift_size;
			middle += shift_size;
			 right += shift_size;
		}
		subsize *= 2;
		shift_size *= 2;

	}
	int remainder = (size % p);
	merge(arr, 0, (size - remainder - 1), (size - 1) );
	return;
}

//Basic merge sort for integers taken from https://www.geeksforgeeks.org/merge-sort/
// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(int arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 =  r - m;

    /* create temp arrays */
    int L[n1], R[n2];

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1+ j];

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort(int arr[], int l, int r)
{
    if (l < r)
    {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l+(r-l)/2;

        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m+1, r);

        merge(arr, l, m, r);
    }
}
