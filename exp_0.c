/**
 *  @author nirellabs
 *  @details The target of this experiment is to demonstrate that in reality it is possible to improve the sorting time 
 *  while having the same theory complexity time by adding some heuristics function.
 *  Here we analyze the time complexity of the following comparison-based sorting algorithms:
 *  Insertion Sort, Merge Sort, Hybrid Sort, Quick Sort and some heuristics function like Median Of Three.
 *  This script prints as output the average time that choosen algorithm takes to sort the array for n dimension,
 *  and a boolean as result of the antagonist function to verify if every array is sorted correctly.
 * 
**/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

time_t seed = 13;               
const int min = 10;             
const int max = 10000;           
const int nExp = 10;            
const int maxRand = 1000000;    
const int granularity = 10;     
const int thresholdHS = 100;      
const int thresholdQS = 30;

typedef struct result{
    clock_t time;
    bool isSorted;
} result;

bool isSorted(const int* A, const int n)
{
    for (int i = 0; i < n - 1; i++) 
        if (A[i] > A[i + 1]) 
            return false;                   
    return true;
}

void generateRanomArray(int* A, const int n)
{
    for (int i = 0; i < n; i++)
        A[i] = rand() % maxRand + 1;
}

void printArray(const int* A, const int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", A[i]);
    printf("\n");
}

int copyfrom(int* A, int i)
{
    return A[i];
}

void swapValue(int* A, const int i, const int j)
{
    int temp = A[i];
    A[i] = A[j];
    A[j] = temp;
}
/*
    SORTING ALGORITHMS
*/
//O(n^2)
void insertionSort(int* A, const int n)
{
    int i, j, key;

    for (i = 1; i < n; i++)
    {
        key = A[i]; 
        j = i - 1;  
        while ((j >= 0) && (A[j] > key))
        {
            A[j + 1] = A[j];
            j--;
        }
        A[j + 1] = key;
    }
}

void adaptedInsertionSort(int* A, const int left, const int right)
{
    int i, j, key;

    for (i = left + 1; i <= right; i++)
    {
        key = A[i];
        j = i - 1;
        while ((j >= left) && (A[j] > key))
        {
            A[j + 1] = A[j];
            j--;
        }
        A[j + 1] = key;
    }
}

//O(n*log(n))
void merge(int* A, const int left, const int mid, const int right)
{   
    int i, j, key,
        n1 = mid - left + 1, 
        n2 = right - mid,    
        L[n1], R[n2];        

    for (i = 0; i < n1; i++)
        L[i] = A[left + i];

    for (j = 0; j < n2; j++)
        R[j] = A[mid + j + 1];
   
    i = 0;
    j = 0;

    for(key = left; key <= right; key++)
    {
        if(i < n1)
        {
            if(j < n2)
            {
                if(L[i] <= R[j])
                {
                    A[key] = copyfrom(L, i);
                    i++;
                }
                else
                {
                    A[key] = copyfrom(R, j); 
                    j++;
                }
            }
            else
            {
                A[key] = copyfrom(L, i);
                i++;
            }
        }
        else
        {
            A[key] = copyfrom(R, j);
            j++;
        }
    }
}

void mergeSort(int* A, const int p, const int r)
{
    if(p < r)
    {
        int q = p + (r - p) / 2;
        mergeSort(A, p, q);
        mergeSort(A, q + 1, r);
        merge(A, p, q, r);
    }
}
//Insertion Sort and Merge Sort combined
void hybridSort(int* A, const int p, const int r)
{
    if (r - p < thresholdHS)
        adaptedInsertionSort(A, p, r);
    else
    {
        int q = (r + p) / 2;
        hybridSort(A, p, q);
        hybridSort(A, q + 1, r);
        merge(A, p, q, r);
    }
}
//Worst case O(n^2), best and average case O(n*log(n))
int partition(int* A, const int p, const int r)
{
    int j, 
        temp = A[r],  
        i = p - 1;
    
    for(j = p; j <= r - 1; j++)
    {
        if(A[j] <= temp)
        {
            i = i + 1;
            swapValue(A, i, j);
        }
    }
    swapValue(A, i + 1, r);
    return i + 1;
}

int medianOfThree(int* A, const int i, const int j, const int k)
{
    if(A[i] > A[j])
    {
        if(A[j] >= A[k])
            return j;
        else if(A[i] < A[k])
            return i;
        else
            return k;
    }
    else
    {
        if(A[i] >= A[k])
            return i;
        else if(A[j] < A[k])
            return j;
        else
            return k;
    }
}

int medianOfThreePartition(int* A, const int p, const int r)
{
    int i, q, temp;

    q = medianOfThree(A, p, r, (p + r) / 2);
    swapValue(A, q, r);
    temp = A[r];
    i = p - 1;
    
    for(int j = p; j <= r - 1; j++)
    {
        if(A[j] <= temp)
        {
            i = i + 1;
            swapValue(A, i, j);
        }
    }
    swapValue(A, i + 1, r);
    return i + 1;
}

void medianOfThreeQuickSort(int* A, const int p, const int r)
{
    if(p < r)
    {
        int q = medianOfThreePartition(A, p, r);
        medianOfThreeQuickSort(A, p, q - 1);
        medianOfThreeQuickSort(A, q + 1, r);
    }
}

void tailQuickSort(int* A, int p, const int r)
{
    if(r - p < thresholdQS)
         adaptedInsertionSort(A, p, r);
    else
        while (p < r)
        {
            int q = medianOfThreePartition(A, p, r);
            tailQuickSort(A, p, q - 1);
            p = q + 1;
        }
}

result sortingAlgorithm(const int* A, const int dim, const char* algorithm)
{
    result result = {0, true};

    clock_t startTime = 0;
    clock_t endTime = 0;
    clock_t totTime = 0;

    int* randomArray = malloc(dim * sizeof(int));

    for(int i = 0; i < dim; i++)
        randomArray[i] = A[i];
    
    startTime = clock();

    if(strcmp(algorithm, "insertionSort") == 0)
        insertionSort(randomArray, dim);
    else if(strcmp(algorithm, "mergeSort") == 0)
        mergeSort(randomArray, 0, dim - 1);
    else if(strcmp(algorithm, "hybridSort") == 0)
        hybridSort(randomArray, 0, dim - 1);
    else if(strcmp(algorithm, "quickSort") == 0)
        medianOfThreeQuickSort(randomArray, 0, dim - 1);
    else if(strcmp(algorithm, "tailQuickSort") == 0)
        tailQuickSort(randomArray, 0, dim - 1);
    else {
        printf("Error: %s algorithm not found", algorithm);
        exit(1);
    }

    endTime = clock();

    totTime = endTime - startTime;

    result.time = totTime;

    if (!isSorted(randomArray, dim))
        result.isSorted = false;

    //printArray(randomArray, dim);

    free(randomArray);

    return result;
}

int main()
{
    srand(seed);

    FILE* res = fopen("result2.csv","w");

    bool isSortedIS = true;
    bool isSortedMS = true;
    bool isSortedHS = true;
    bool isSortedQS = true;
    bool isSortedTQ = true;

    clock_t timeIS = 0;
    clock_t timeMS = 0;
    clock_t timeHS = 0;
    clock_t timeQS = 0;
    clock_t timeTQ = 0;

    result resultIS;
    result resultMS;
    result resultHS;
    result resultQS;
    result resultTQ;

    int* array = malloc(max * sizeof(int));

    printf("+-----------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+\n");
    printf("| ######### | InsertionSort                 | MergeSort                     | HybridSort                    | QuickSort with MedianOfThree  | TailQuickSort                 |\n");
    printf("+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+\n");
    printf("| Dimension | Time              | isSorted? | Time              | isSorted? | Time              | isSorted? | Time              | isSorted? | Time              | isSorted? |\n");
    printf("+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+-------------------------------+-------------------+-----------+\n");
    
    fprintf(res, "Dimension\tInsertionSort\tMegeSort\tHybridSort\tQuickSortMOT\tTailQuickSort\n");
       
    for(int dim = min; dim <= max; dim += granularity)
    {
        timeIS = 0;
        timeMS = 0;
        timeHS = 0;
        timeQS = 0;
        timeTQ = 0;

        isSortedIS = true;  
        isSortedMS = true;  
        isSortedHS = true;  
        isSortedQS = true;
        isSortedTQ = true;

        for(int exp = 0; exp < nExp; exp++)
        {
            generateRanomArray(array, dim);

            resultIS = sortingAlgorithm(array, dim, "insertionSort");
            timeIS += resultIS.time;
            isSortedIS = resultIS.isSorted;

            resultMS = sortingAlgorithm(array, dim, "mergeSort");
            timeMS += resultMS.time;
            isSortedMS = resultMS.isSorted;

            resultHS = sortingAlgorithm(array, dim, "hybridSort");
            timeHS += resultHS.time;
            isSortedHS = resultHS.isSorted;

            resultQS = sortingAlgorithm(array, dim, "quickSort");
            timeQS += resultQS.time;
            isSortedQS = resultQS.isSorted;

            resultTQ = sortingAlgorithm(array, dim, "tailQuickSort");
            timeTQ += resultTQ.time;
            isSortedTQ = resultTQ.isSorted;            
        }

        printf("| %9d | %17f | %9s | %17f | %9s | %17f | %9s | %17f | %9s | %17f | %9s |\n", dim, 
        (double) timeIS/nExp, isSortedIS ? "true" : "false",
        (double) timeMS/nExp, isSortedMS ? "true" : "false",
        (double) timeHS/nExp, isSortedHS ? "true" : "false",
        (double) timeQS/nExp, isSortedQS ? "true" : "false",
        (double) timeTQ/nExp, isSortedTQ ? "true" : "false");
       
        fprintf(res, "%d\t%f\t%f\t%f\t%f\t%f\n", dim, 
        (double) timeIS/nExp,
        (double) timeMS/nExp,
        (double) timeHS/nExp,
        (double) timeQS/nExp,
        (double) timeTQ/nExp);
    }
    free(array);
    fclose(res);
    printf("+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+\n");
    
    return 0;
}