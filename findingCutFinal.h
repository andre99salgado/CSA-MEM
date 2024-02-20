#ifndef INC_FINDINGCUTFINAL_H // Inclusion guard.
#define INC_FINDINGCUTFINAL_H // Inclusion guard.

// Define a structure to represent a MEM block with start and end positions.
struct Block
{
    size_t start;
    size_t end;
};

struct CombinedMEM
{
    int posRef;
    int *posQuery;
    int len;
};

struct CombinedREAD
{
    char name[256];
    int len;
    struct CombinedMEM **mems;
};

// Function to compute the intersection of two blocks.
struct Block get_intersection(struct Block b1, struct Block b2)
{
    struct Block intersection;

    // Determine the start and end of the intersection
    intersection.start = (b1.start < b2.start) ? b2.start : b1.start;
    intersection.end = (b1.end < b2.end) ? b1.end : b2.end;

    // If there's no overlap, adjust the intersection to be empty
    if (intersection.start > intersection.end)
        intersection.start = intersection.end = 0;

    return intersection;
}

void removeDuplicates(struct READ *organism, const char *removeType)
{
    int n = organism->len;                     // Number of elements in the mems array
    struct MEM **sortedArray = organism->mems; // Assuming you store an array of pointers to struct MEM
    if (n <= 0)
    {
        return; // No elements or invalid input
    }

    struct MEM **uniqueArray = malloc(n * sizeof(struct MEM *));
    if (uniqueArray == NULL)
    {
        printf("Memory allocation failed\n");
        return;
    }

    int uniqueCount = 0;

    uniqueArray[uniqueCount++] = sortedArray[0];

    if (strcmp(removeType, "posQuery") == 0)
    {
        for (int i = 1; i < n; i++)
        {
            struct Block b1 = {sortedArray[i]->posQuery, (sortedArray[i]->posQuery + sortedArray[i]->len)},
                         b2 = {uniqueArray[uniqueCount - 1]->posQuery, (uniqueArray[uniqueCount - 1]->posQuery + uniqueArray[uniqueCount - 1]->len)};
            struct Block intersection = get_intersection(b1, b2);
            if (intersection.end == 0 && intersection.start == 0)
            {
                uniqueArray[uniqueCount++] = sortedArray[i];
            }
            else if (sortedArray[i]->len > uniqueArray[uniqueCount - 1]->len)
            {
                uniqueArray[uniqueCount - 1] = sortedArray[i];
            }
        }
    }
    else if (strcmp(removeType, "posRef") == 0)
    {
        for (int i = 1; i < n; i++)
        {
            struct Block b1 = {sortedArray[i]->posRef, (sortedArray[i]->posRef + sortedArray[i]->len)},
                         b2 = {uniqueArray[uniqueCount - 1]->posRef, (uniqueArray[uniqueCount - 1]->posRef + uniqueArray[uniqueCount - 1]->len)};
            struct Block intersection = get_intersection(b1, b2);
            if (intersection.end == 0 && intersection.start == 0)
            {
                uniqueArray[uniqueCount++] = sortedArray[i];
            }
            else if (sortedArray[i]->len > uniqueArray[uniqueCount - 1]->len)
            {
                uniqueArray[uniqueCount - 1] = sortedArray[i];
            }
        }
    }
    else
    {
        printf("Invalid remove type argument\n");
        return;
    }

    // Update the unique elements back to organism->mems
    struct MEM **newMems = malloc(uniqueCount * sizeof(struct MEM *));
    if (newMems == NULL)
    {
        printf("Memory allocation failed\n");
        return;
    }

    for (int i = 0; i < uniqueCount; i++)
    {
        newMems[i] = malloc(sizeof(struct MEM));
        *(newMems[i]) = *(uniqueArray[i]);
    }

    // Free the old mems array
    for (int i = 0; i < n; i++)
    {
        free(sortedArray[i]);
    }
    free(sortedArray);

    // Update the organism's mems array
    organism->mems = newMems;
    organism->len = uniqueCount;

    // // Print the unique elements
    // for (int i = 0; i < uniqueCount; i++)
    // {
    //     printf("{%d, %d, %d}\n", newMems[i]->posRef, newMems[i]->posQuery, newMems[i]->len);
    // }

    free(uniqueArray); // Free the allocated memory
}

// Comparison function for qsort
int compareMems(const void *a, const void *b)
{
    struct MEM *memA = *(struct MEM **)a;
    struct MEM *memB = *(struct MEM **)b;
    
    return memA->posRef - memB->posRef;
}

// Comparison function for qsort
int compareCombinedMems(const void *a, const void *b)
{
    struct CombinedMEM *memA = *(struct CombinedMEM **)a;
    struct CombinedMEM *memB = *(struct CombinedMEM **)b;
    
    return memA->posRef - memB->posRef;
}


void sortMemsByPosRef(struct READ *organism)
{
    qsort(organism->mems, organism->len, sizeof(struct MEM *), compareMems);
}

void sortCombinedMemsByPosRef(struct CombinedREAD *organism)
{
    qsort(organism->mems, organism->len, sizeof(struct CombinedMEM *), compareCombinedMems);
}

struct CombinedMEM * init_CombinedMEM(size_t posRef, size_t len, int numReads){
    struct CombinedMEM *temp = (struct CombinedMEM *)malloc(sizeof(struct CombinedMEM));
    temp->posRef = posRef;
    temp->len = len;
    temp->posQuery = (int *)malloc(numReads*sizeof(int));
    return temp;
}

struct CombinedREAD * initCombinedREAD(struct READ *organism, int mem_number, int numReads ){
    struct CombinedREAD *combinedRead = (struct CombinedREAD *)malloc(sizeof(struct CombinedREAD));
    combinedRead->mems = (struct CombinedMEM **)malloc(mem_number * sizeof(struct CombinedMEM*));
    combinedRead->len = mem_number;
    for(size_t i = 0; i < mem_number; i++)
        combinedRead->mems[i] = init_CombinedMEM(0, 0,numReads);
    
    size_t combinedRead_counter = 0;
    for(size_t i = 0; i < 2; i++){
        for(size_t j = 0; j < organism[i].len; j++){
            combinedRead->mems[combinedRead_counter]->posRef = organism[i].mems[j]->posRef;
            combinedRead->mems[combinedRead_counter]->len = organism[i].mems[j]->len;
            combinedRead->mems[combinedRead_counter]->posQuery[i] = organism[i].mems[j]->posQuery;
            combinedRead_counter++;
        }
    }
    return combinedRead;
}

struct CombinedREAD * addCombinedREAD(struct READ *organism, struct CombinedREAD * combinedRead, int numReads, int index ){
    size_t i = combinedRead->len;
    combinedRead->len += organism->len;
    combinedRead->mems= (struct CombinedMEM **)realloc(combinedRead->mems,combinedRead->len*sizeof(struct CombinedMEM*));
    for(size_t j = 0; j < organism->len; j++){
        combinedRead->mems[i] = init_CombinedMEM( organism->mems[j]->posRef, organism->mems[j]->len, numReads);
        combinedRead->mems[i]->posQuery[index] = organism->mems[j]->posQuery;
        i++;
    }
        
    return combinedRead;
}

// struct CombinedREAD * combineOverlapping(struct CombinedREAD * combinedRead, size_t numReads, size_t index){
//     struct CombinedREAD * temp_combinedRead = (struct CombinedREAD *)malloc(sizeof(struct CombinedREAD));
//     temp_combinedRead->mems = (struct CombinedMEM **)malloc(sizeof(struct CombinedMEM*));
//     temp_combinedRead->len = 0;
//     temp_combinedRead->mems[0] = init_CombinedMEM(0, 0, 0);
//     size_t counter = 0;
//     //for (size_t i = 0; i < numReads; i++)
//     //    temp_combinedRead->mems[0]->posQuery[i] = combinedRead->mems[0]->posQuery[i];
    
//     for (size_t i = 1; i < combinedRead->len; i++)
//     {
//         /* code */
//         struct Block b1 = {combinedRead->mems[i-1]->posRef, (combinedRead->mems[i-1]->posRef + combinedRead->mems[i-1]->len)},
//                  b2 = {combinedRead->mems[i]->posRef, (combinedRead->mems[i]->posRef + combinedRead->mems[i]->len)},
//                  b3 = {temp_combinedRead->mems[counter]->posRef, temp_combinedRead->mems[counter]->posRef + temp_combinedRead->mems[counter]->len};
//         struct Block intersection = get_intersection(b2, b3);
//         if(intersection.end != 0 && intersection.start != 0){
//             temp_combinedRead->mems[counter] = init_CombinedMEM(intersection.start,intersection.end-intersection.start,numReads);
//             for(size_t j=0;j<numReads;j++){
//                 if(temp_combinedRead->mems[counter]->posQuery[j] != 0)
//                     temp_combinedRead->mems[counter]->posQuery[j] += b2.start - b3.start;
//                 if(combinedRead->mems[i]->posQuery[j] != 0)
//                     temp_combinedRead->mems[counter]->posQuery[j] = combinedRead->mems[i]->posQuery[j];
//             }
//             //adiciona à lista combinedtemp_combinedReadRead sem aumentar o counter
//         }
//         intersection = get_intersection(b2, b1);
//         if(intersection.end != 0 && intersection.start != 0){
//             temp_combinedRead->len += 1;
//             temp_combinedRead->mems = (struct CombinedMEM **)realloc(temp_combinedRead->mems, temp_combinedRead->len * sizeof(struct CombinedMEM*));
//             counter = temp_combinedRead->len - 1;
//             temp_combinedRead->mems[temp_combinedRead->len-1] = init_CombinedMEM(intersection.start,intersection.end-intersection.start,numReads);
//             for(size_t j=0;j<numReads;j++){
//                 if(combinedRead->mems[i-1]->posQuery[j] != 0)
//                     temp_combinedRead->mems[temp_combinedRead->len-1]->posQuery[j] = combinedRead->mems[i-1]->posQuery[j] + b2.start - b1.start;
                    
//                 if(combinedRead->mems[i]->posQuery[j] != 0)
//                     temp_combinedRead->mems[temp_combinedRead->len-1]->posQuery[j] =  combinedRead->mems[i]->posQuery[j];
//             }
            
//             //adiciona à lista combinedtemp_combinedReadRead aumentando antes o counter
//         } 
//     }
//     return temp_combinedRead;
// }


struct CombinedREAD * combineOverlapping(struct CombinedREAD * combinedRead, size_t numReads, size_t index){
    struct CombinedREAD * temp_combinedRead = (struct CombinedREAD *)malloc(sizeof(struct CombinedREAD));
    temp_combinedRead->mems = (struct CombinedMEM **)malloc(sizeof(struct CombinedMEM*));
    temp_combinedRead->len = 0;
    temp_combinedRead->mems[0] = init_CombinedMEM(0, 0, 0);
    size_t counter = 0;
    //for (size_t i = 0; i < numReads; i++)
    //    temp_combinedRead->mems[0]->posQuery[i] = combinedRead->mems[0]->posQuery[i];
    
    for (size_t i = 1; i < combinedRead->len; i++)
    {
        /* code */
        struct Block b1 = {combinedRead->mems[i-1]->posRef, (combinedRead->mems[i-1]->posRef + combinedRead->mems[i-1]->len)},
                 b2 = {combinedRead->mems[i]->posRef, (combinedRead->mems[i]->posRef + combinedRead->mems[i]->len)};
        struct Block intersection = get_intersection(b2, b1);
        if(intersection.end != 0 && intersection.start != 0){
            temp_combinedRead->len += 1;
            temp_combinedRead->mems = (struct CombinedMEM **)realloc(temp_combinedRead->mems, temp_combinedRead->len * sizeof(struct CombinedMEM*));
            temp_combinedRead->mems[temp_combinedRead->len-1] = init_CombinedMEM(intersection.start,intersection.end-intersection.start,numReads);
            for(size_t j=0;j<numReads;j++){
                if(combinedRead->mems[i-1]->posQuery[j] != 0)
                    temp_combinedRead->mems[temp_combinedRead->len-1]->posQuery[j] = combinedRead->mems[i-1]->posQuery[j] + b2.start - b1.start;
                if(combinedRead->mems[i]->posQuery[j] != 0)
                    temp_combinedRead->mems[temp_combinedRead->len-1]->posQuery[j] =  combinedRead->mems[i]->posQuery[j];                    
            }
        } 
    }
    return temp_combinedRead;
}

struct CombinedMEM * getHeaviestPath(struct CombinedREAD * combinedRead, size_t numReads) {
    struct CombinedMEM* starting_mem = (struct CombinedMEM*)malloc(sizeof(struct CombinedMEM));
    starting_mem = init_CombinedMEM(0,0,numReads); 
    size_t heaviest_sum = 0, current_sum = 0;

    for (size_t i = 0; i < combinedRead->len; i++)
    {  
        current_sum = combinedRead->mems[i]->len;
        size_t last_i= i;
        for (size_t j = i+1; j < combinedRead->len; j++)
        {  
            int flag = 0;
            for(size_t k = 0; k < numReads; k++){
                if(combinedRead->mems[j]->posQuery[k] < combinedRead->mems[j-1]->posQuery[k]){
                    flag = 1;
                    break;
                }
            }

            if(flag == 1){
                i = j;
                break;
            }
            current_sum += combinedRead->mems[j]->len;
        }
        if(current_sum > heaviest_sum){
            starting_mem = combinedRead->mems[last_i];
            heaviest_sum = current_sum;
        }
    } 
    printf("\n[%d", starting_mem->posRef);
    for(size_t i = 0; i < numReads; i++){
        printf(",%d",starting_mem->posQuery[i]);
    }
    printf("]\n");
    return starting_mem;
}

size_t getCutPosition(struct READ *organism, size_t numReads)
{   
    size_t number_of_unique_mems = 0;
    struct CombinedREAD * c_read = NULL;
    for(size_t i = 0; i < numReads; i++){
        //removeDuplicates(&organism[i], "posQuery");
        sortMemsByPosRef(&organism[i]);
        removeDuplicates(&organism[i], "posRef");
        number_of_unique_mems += organism[i].len;
        if(i == 1){
            c_read = initCombinedREAD(organism, number_of_unique_mems, numReads);
            sortCombinedMemsByPosRef(c_read);
        }
        if(i > 1){
            addCombinedREAD(&organism[i], c_read, numReads, i);
            sortCombinedMemsByPosRef(c_read);
        }
        if(i > 0){
            c_read = combineOverlapping(c_read, numReads,i);
            int oioioi = 0;
            getHeaviestPath(c_read, numReads);
        }
        
    }
    getHeaviestPath(c_read, numReads);
    return 0;
}

#endif // Inclusion guard.