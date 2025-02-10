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



///////////////////////////////////////////////

struct Node {
    int * value;
    int weight;
    struct Node *pointer;
    long index;
};

// Function to insert a new object into the list
void inserttt(struct Node **L, int * o, int wi, long index, size_t numbReads){
    struct Node *newNode = (struct Node *)malloc(sizeof(struct Node));
    newNode->value = (int *)malloc(numbReads*sizeof(int));
    for(size_t i=0; i<numbReads; i++)
        newNode->value[i] = o[i];
        
    newNode->weight = wi;
    newNode->pointer = NULL;
    newNode->index = index;

    if (*L == NULL) {
        *L = newNode;
    } else {
        struct Node *current = *L;
        
        int bol = 0;
        for(size_t i = 0; i < numbReads; i++)
            if(current->value[i] <= newNode->value[i])
                bol = 1;
        
        if(bol == 0) {
                newNode->pointer = current;
                *L = newNode;
                return ;
            }

        while (current->pointer != NULL) {
            bol = 0;
            for(size_t i = 0; i < numbReads; i++)
                if(current->pointer->value[i] <= o[i])
                    bol = 1;

            if(bol == 0){
                newNode->pointer = current->pointer;
                current->pointer = newNode;
                return ;
            }
            current = current->pointer;
            
        }
        current->pointer = newNode;
    }
}

// Function to delete an object from the list
void delete(struct Node **L, int *o, int wi, size_t numbReads) {
    struct Node *current = *L;
    struct Node *prev = NULL;

    int bol = 0;
        for(size_t i = 0; i < numbReads; i++)
            if(current->value[i] != o[i])
                bol = 1;

    while (current != NULL && (bol == 1 || current->weight != wi)) {
        prev = current;
        current = current->pointer;

        bol = 0;
        for(size_t i = 0; i < numbReads; i++)
            if(current->value[i] != o[i])
                bol = 1;
    }

    if (current == NULL) {
        return; // Object not found
    }

    if (prev == NULL) {
        *L = current->pointer;
    } else {
        prev->pointer = current->pointer;
    }

    free(current);
}

// Function to find the least element strictly larger than o in the list
struct Node *next(struct Node *L, int *o, int wi, size_t numbReads) {
    struct Node *current = L;
    struct Node *result = NULL;

    while (current != NULL) {
        int bol = 0;
        for(size_t i = 0; i < numbReads; i++)
            if(current->value[i] <= o[i])
                bol = 1;

        if (bol == 0) {
            return current;
        }
        current = current->pointer;
    }

    return result;
}

// Function to find the largest element strictly smaller than o in the list
struct Node *prev(struct Node *L, int *o, int wi, size_t numbReads) {
    struct Node *current = L;
    struct Node *result = NULL;
    
    int bol = 0;
    if(current != NULL){
        for(size_t i = 0; i < numbReads; i++)
            if(current->value[i] >= o[i])
                return NULL;
    }
        
    while (current != NULL) {
        
        if (current->pointer == NULL)
            return current;

        bol = 0;
        for(size_t i = 0; i < numbReads; i++)
            if(current->pointer->value[i] > o[i])
                return current;

        current = current->pointer;
    }

    return result;
}

// Function to find the largest element strictly smaller than o in the list
struct Node *new_prev(struct Node *L, int *o, int wi, size_t numbReads) {
    struct Node *current = L;
    struct Node *result = NULL;

    int bol = 0;
    while (current != NULL) {
    
        bol = 0;
        size_t median = 0;
        for(size_t i = 0; i < numbReads; i++){
            if(current->value[i] > o[i]){
                bol = 1;
                break;
            }  
            if(i > 1 && ( (o[i] - current->value[i]) < median/2 ||  (o[i] - current->value[i]) > ((median/2)+median) )){
               bol = 1;
               break;
            }

            if(bol == 0)median = o[i] - current->value[i];
        }
            

        if(bol == 0){
            result = current;
        }

            

        current = current->pointer;
    }

    return result;
}


// Function to find the maximal element in the list
struct Node *max(struct Node *L) {
    struct Node *current = L;
    struct Node *result = NULL;

    size_t max_weigth = 0;
    
    while (current != NULL) {
        if(current->weight > max_weigth){
            max_weigth = current->weight;
            result = current;
        }
        current = current->pointer;
    }
    return result;
}

// Function to find the minimal element in the list
struct Node *min(struct Node *L) {
    return L;
}

// Function to free the memory of the linked list
void freeList(struct Node *L) {
    while (L != NULL) {
        struct Node *temp = L;
        L = L->pointer;
        free(temp);
    }
}

// Function to construct a new node
struct Node *newNode(int *value, int weight, struct Node *pointer, size_t numbReads) {
    struct Node *newNode = (struct Node *)malloc(sizeof(struct Node));
    newNode->weight = weight;
    newNode->value = (int *)malloc(numbReads*sizeof(int));
    for (size_t i = 0; i < numbReads; i++)
    {
        newNode->value[i] = value[i];
    }
    newNode->pointer = pointer;
    if(pointer != NULL)
        newNode->index = pointer->index;
    else newNode->index = -1;
    return newNode;
}

size_t new_hcis(struct Node *nodes, size_t n, size_t numbReads) {
    struct Node *L = NULL;
    struct Node **node = (struct Node **)malloc(n * sizeof(struct Node *));
    
    for (int i = 0; i < n; i++) {
        node[i] = NULL; // Initialize the node array
    }

    for (int i = 0; i < n; i++) {

        int *oi = (int*)malloc(numbReads * sizeof(int));
        int wi = 0;
        for (size_t ji = 0; ji < numbReads; ji++)
        {
            oi[ji] = nodes[i].value[ji];
        }
        wi = nodes[i].weight;
        
        struct Node *s = new_prev(L, oi, wi, numbReads);
        struct Node *t = NULL;

        if (s == NULL) {
            t = min(L);
        }else{
            t = next(L, s->value, s->weight, numbReads);
        }        

        int temp_weight = 0;
            if(s != NULL) 
                temp_weight = s->weight;

        printf("added: %d\n", (wi + temp_weight));
        inserttt(&L, oi,wi + temp_weight, i, numbReads);
        node[i] = newNode(oi,wi + temp_weight, s, numbReads);
    
    }

    // Find the maximal element in L
    struct Node *maxNode = max(L);
    
    // Recover the LIS using the node array
    long index = maxNode->index;
    size_t sum = maxNode->weight;
    size_t return_pos = 0;
    size_t coun_sum = 0;
    while (index != -1) {
        printf("%d\n", index);

        printf("\n[_,");
        for(size_t itee=0;itee<numbReads;itee++){
            printf("%zu,", nodes[index].value[itee]);
        }
        coun_sum += nodes[index].weight;
        printf("]\n");
        
        if(node[index]->index == -1){
            printf("Cut position with Heaviest increasing subsequence: %zu\n", index);
            return_pos = index;
        }
            
        index = node[index]->index;
    }
    printf("A soma: %zu\n", sum);

    // Free the memory of the list and node array
    freeList(L);
    free(node);

    return return_pos;
}

size_t his(struct Node *nodes, size_t n, size_t numbReads) {
    struct Node *L = NULL;
    struct Node **node = (struct Node **)malloc(n * sizeof(struct Node *));
    
    for (int i = 0; i < n; i++) {
        node[i] = NULL; // Initialize the node array
    }

    for (int i = 0; i < n; i++) {

        int *oi = (int*)malloc(numbReads * sizeof(int));
        int wi = 0;
        for (size_t ji = 0; ji < numbReads; ji++)
        {
            oi[ji] = nodes[i].value[ji];
        }
        wi = nodes[i].weight;
        
        struct Node *s = prev(L, oi, wi, numbReads);
        struct Node *t = NULL;

        if (s == NULL) {
            t = min(L);
        }else{
            t = next(L, s->value, s->weight, numbReads);
        }

        while(t != NULL){
            int temp_weight = 0;
            if(s != NULL) temp_weight = s->weight;
            if(temp_weight + wi < t->weight)
                break;
            
            //ver se tenho de fazer o for para copiar
            int *last_value = t->value;
            int last_weight = t->weight;

            delete(&L, t->value, t->weight, numbReads);
            t = next(L, last_value, last_weight, numbReads);
        }
       

        int bol = 0;
        long long median = 0;

        if(t != NULL){
            for (size_t jiterador = 0; jiterador < numbReads; jiterador++){
                if(t->value[jiterador] < oi[jiterador]){
                    median = -1;
                    break;
                }
                median += t->value[jiterador] - oi[jiterador];
                printf("%zu,",t->value[jiterador] - oi[jiterador]);
            }
            printf("]\n");
            if(median != -1)
                median = median / numbReads;

            for(size_t ite2 = 0; ite2 < numbReads; ite2++)
                if( (oi[ite2] >= t->value[ite2] || median == -1 || t->value[ite2] - oi[ite2] < median/2 || t->value[ite2] - oi[ite2] > median + (median/2)))
                    bol = 1;
        }
        
        int flag_deep = 0;
        if(s != NULL)
        if(s->weight > wi){
            flag_deep = 1;
            for (size_t jiterador = 0; jiterador < numbReads; jiterador++){
                if(s->value[jiterador] > oi[jiterador]){
                    median = -1;
                    break;
                }
                median += oi[jiterador] -  s->value[jiterador];
                printf("%zu,", oi[jiterador] - s->value[jiterador]);
            }
            printf("]\n");
            if(median != -1)
                median = median / numbReads;

            for(size_t ite2 = 0; ite2 < numbReads; ite2++)
                if( median == -1 ||  oi[ite2] - s->value[ite2] < median/2 ||  oi[ite2] - s->value[ite2] > median + (median/2))
                    bol = 1;
        }
        
        

        if(t == NULL && bol == 0){
            int temp_weight = 0;
            if(s != NULL) 
                temp_weight = s->weight;

            if(flag_deep==0){
                inserttt(&L, oi,wi, i, numbReads);
                node[i] = newNode(oi,wi, NULL, numbReads);
            }else{
                inserttt(&L, oi,wi + temp_weight, i, numbReads);
                node[i] = newNode(oi,wi + temp_weight, s, numbReads);
            }

            
        }

        int teste = 0;
    }

    // Find the maximal element in L
    struct Node *maxNode = max(L);
    
    // Recover the LIS using the node array
    long index = maxNode->index;
    size_t sum = maxNode->weight;
    size_t return_pos = 0;
    while (index != -1) {
        printf("%d\n", index);
        if(node[index]->index == -1){
            printf("Cut position with Heaviest increasing subsequence: %zu\n", index);
            return_pos = index;
        }
            
        index = node[index]->index;
    }
    printf("A soma: %zu\n", sum);

    // Free the memory of the list and node array
    freeList(L);
    free(node);

    return return_pos;
}


///////////////////////////////////////////////

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
        printf(",%zu",starting_mem->posQuery[i]);
    }
    printf("]\n");
    return starting_mem;
}


size_t get_sequenced(struct Node *nodes, size_t n, size_t numbReads){
    size_t heaviest_sum = 0, current_sum = 0;
    size_t chosen_pos = 0;
    size_t current_start = 0;
    for (size_t i = 0; i < n; i++)
    {   
        current_sum += nodes[i].weight;
        if(i == n - 1){
            if(current_sum > heaviest_sum){
                    heaviest_sum = current_sum;
                    chosen_pos = current_start;
                }
            
            printf("Heaviest sum: %zu\n", heaviest_sum);
            printf("Chosen pos with heaviest increasing consecutive sequence: %zu\n", chosen_pos);
            return chosen_pos;
        }

        long long median = 0;
        //printf("\n[");

        for (size_t j = 0; j < numbReads; j++){
            if(nodes[i+1].value[j] < nodes[i].value[j]){
                median = -1;
                break;
            }
            median += nodes[i+1].value[j] - nodes[i].value[j];
            //printf("%zu,",nodes[i+1].value[j] - nodes[i].value[j]);
        }
        //printf("]\n");
        if(median != -1)
            median = median / numbReads;

        for (size_t j = 0; j < numbReads; j++)
        {   
            if(nodes[i].value[j] >= nodes[i+1].value[j] || median == -1 || nodes[i+1].value[j] - nodes[i].value[j] < median/2 || nodes[i+1].value[j] - nodes[i].value[j] > median + (median/2)){
                if(current_sum > heaviest_sum){
                    heaviest_sum = current_sum;
                    chosen_pos = current_start;
                }
                current_sum = 0;
                current_start = i+1;
                break;
            }
            /* code */
        }
        /* code */
    }
    
}

size_t getCutPosition(struct READ *organism, size_t numReads)
{   
    size_t number_of_unique_mems = 0;
    struct CombinedREAD * c_read = NULL;
    for(size_t i = 0; i < numReads; i++){
        removeDuplicates(&organism[i], "posQuery");
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
        }
    }
    struct Node nodeList[c_read->len];
    size_t sumtemp = 0;
    for(size_t iterador=0;iterador<c_read->len;iterador++){
        nodeList[iterador].index = iterador;
    	nodeList[iterador].pointer = NULL;
        nodeList[iterador].value = (int*)malloc(numReads*sizeof(int));

        printf("\n[ ");
        for(size_t itee=0;itee<numReads;itee++){
            nodeList[iterador].value[itee] = c_read->mems[iterador]->posQuery[itee];
            printf("%zu,", c_read->mems[iterador]->posQuery[itee]);
        }
        sumtemp +=c_read->mems[iterador]->len;
        printf("]\n");
        printf("%zu\n", c_read->mems[iterador]->len);
        printf("temp_sum:%zu\n", sumtemp);
    	    
    	nodeList[iterador].weight = c_read->mems[iterador]->len;

    }

    //size_t his_pos = his(nodeList, c_read->len, numReads);

    size_t his_pos = new_hcis(nodeList, c_read->len, numReads);
    //size_t numNodes = sizeof(nodeList) / sizeof(nodeList[0]);
    printf("\n[%zu,", c_read->mems[his_pos]->posRef);
        for(size_t itee=0;itee<numReads;itee++){
            nodeList[his_pos].value[itee] = c_read->mems[his_pos]->posQuery[itee];
            printf("%zu,", c_read->mems[his_pos]->posQuery[itee]);
        }
        sumtemp +=c_read->mems[his_pos]->len;
        printf("]\n");
    
    //getHeaviestPath(c_read, numReads);
    size_t not_his_pos = get_sequenced(nodeList, c_read->len, numReads);
    printf("\n[%zu,", c_read->mems[not_his_pos]->posRef);
        for(size_t itee=0;itee<numReads;itee++){
            nodeList[not_his_pos].value[itee] = c_read->mems[not_his_pos]->posQuery[itee];
            printf("%zu,", c_read->mems[not_his_pos]->posQuery[itee]);
        }
        sumtemp +=c_read->mems[not_his_pos]->len;
        printf("]\n");
    
    return 0;
}

#endif // Inclusion guard.
