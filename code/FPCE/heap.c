#include<stdio.h>
#include<stdlib.h>

#define swap(u,v) do{\
    u ^= v; \
    v ^= u; \
    u ^= v; \
}while(0)

typedef struct{
    int *arr;
    int size;
} Heap;

/*
int left(int k) 
{ return 2*k+1; }

int right(int k)
{ return 2*k + 2; }
*/

void heapify(Heap *hp, int i)
{
    int l, r, minIndex=i;
    l = 2*i+1;//left(i);
    r = 2*i+2;//right(i);
    if(l < hp->size && hp->arr[l] < hp->arr[minIndex])
        minIndex = l;
    if(r < hp->size && hp->arr[r] < hp->arr[minIndex])
        minIndex = r;

    if(minIndex != i){
        swap(hp->arr[i], hp->arr[minIndex]);
        heapify(hp, minIndex);
    }
}

void createHeap(Heap *hp, int *arr, int n)
{
    int i;
    hp->size = n;
    hp->arr = arr;//(int *)malloc(sizeof(int) * n);
    for(i = n / 2; i >= 0; i--){
        heapify(hp, i);
    }
}

int extractMin(Heap *hp)
{
    if(hp->size < 1){
        printf("Heap underflow");
        exit(1);
    }
    int min = hp->arr[0];
    hp->arr[0] = hp->arr[hp->size-1];
    hp->size = hp->size - 1;
    heapify(hp, 0);
    return min;
}
/*
void main()
{
    Heap hp;
    int a[] = {6, 5, 3, 9, 2, 11, 7, 1};
    createHeap(&hp, a,8);
    printf("min: %d\n", extractMin(&hp));
    printf("min: %d\n", extractMin(&hp));
    printf("min: %d\n", extractMin(&hp));
    printf("min: %d\n", extractMin(&hp));
}
*/