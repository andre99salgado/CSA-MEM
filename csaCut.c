#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "Reads.h"
#include "Trie.h"
//#include "findingCut.h"
#include "findingCutFinal.h"

/*
Coloco cada READ de MEMs numa lista
Construo uma Trie (pouco memory efficient) para a READ mais pequena
Para procurar as common sequences repito para todas as outras READS:
    Procuro todas as sequencias possiveis na READ
    Procuro-as na TRIE
    Se tiverem presentes marco como presentes
    Removo todos os nodes da trie que não tenham dado match
Por fim, calculo o heaviest path da TRIE
*/

int main(void)
{
    printf("Opening the mem file!\n");

    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    // fp = fopen("e_coli_1303_250-mem.txt", "r");
    //fp = fopen("csa_cut_example.txt", "r");
    //fp = fopen("slaMem/slaMEM/fasta-mems.txt", "r");
    fp = fopen("csa_cut_mem.txt", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    printf("Success!\n");

    int countReads = 0;
    // Allocate space for a single READ object
    struct READ *organism = malloc(sizeof(struct READ));
    int currentRead = 0;

    int number_mems = 0;
    struct AuxStruct *auxStruct = malloc(sizeof(struct AuxStruct));
    auxStruct->maxNumberOfMems = 0;
    auxStruct->idOfRead = 0;

    int tamMin = 2147483647;
    while ((read = getline(&line, &len, fp)) != -1)
    {
        if (line[0] == '>')
        {
            // Found a new read, create a new READ object and add it to the organism
            struct READ newRead;
            init_read(&newRead, line + 1, 0);

            organism = realloc(organism, (currentRead + 1) * sizeof(struct READ));
            organism[currentRead] = newRead;
            currentRead++;

            if (number_mems < tamMin && 0 != number_mems)
            {
                tamMin = number_mems;
            }

            number_mems = 0;
        }
        else
        {
            // Found a MEM, parse it and add it to the current read
            int posRef, posQuery, memLen;
            char *te;
            sscanf(line, "%d %d %d", &posRef, &posQuery, &memLen);

            struct MEM *newMem = (struct MEM *)malloc(sizeof(struct MEM));
            init_mem(newMem, posRef, posQuery, memLen);

            number_mems++;
            if (number_mems > auxStruct->maxNumberOfMems)
            {
                auxStruct->maxNumberOfMems = number_mems;
                auxStruct->idOfRead = currentRead - 1;
            }

            add_mem_to_read(&organism[currentRead - 1], newMem, auxStruct->maxNumberOfMems);

            if (organism[currentRead - 1].mems[organism[currentRead - 1].len - 1]->posRef < 0)
                printf("LINE if NULL value: %s\n", line);
        }
    }

    if (organism[currentRead - 1].len < tamMin)
    {
        tamMin = number_mems;
    }
    printf("Tam max %d \n", organism[auxStruct->idOfRead].len);
    printf("Tam min %d \n", tamMin);

    fclose(fp);
    if (line)
        free(line);

    if (currentRead < 1)
        return 0;

    // struct TrieNode *root = getNode();
    // buildTrie(organism, root);

    // for (int i = 1; i < currentRead; i++)
    // {
    //     compareWithTrie(organism, i, root);
    // }

    // struct MaxPath maxPath2 = findMaxPath(root);


    for (int i = 0; i < currentRead; i++)
    {
        struct READ *currentReadPtr = &organism[i];
        size_t left = 0;
        size_t right = currentReadPtr->len - 1;
        while (left < right)
        {
            // Swap the elements at left and right indices in the mems array
            struct MEM *temp = currentReadPtr->mems[left];
            currentReadPtr->mems[left] = currentReadPtr->mems[right];
            currentReadPtr->mems[right] = temp;

            // Move towards the middle of the array
            left++;
            right--;
        }
    }
     //commonSequences(organism, currentRead);
    getCutPosition(organism, currentRead);

    print_organism_to_file(organism, currentRead, "output.txt");

    // free_organism(organism, currentRead);
    return 1;
}
