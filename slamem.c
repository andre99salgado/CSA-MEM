/* ================================================================= *
 *  slamem.c : Main application                                      *
 *                                                                   *
 *  slaMEM: MUMmer-like tool to retrieve Maximum Exact Matches using *
 *          an FM-Index and a Sampled Longest Common Prefix Array    *
 *                                                                   *
 *  Copyright (c) 2013, Francisco Fernandes <fjdf@kdbio.inesc-id.pt> *
 *  Knowledge Discovery in Bioinformatics group (KDBIO/INESC-ID)     *
 *  All rights reserved                                              *
 *                                                                   *
 *  This file is subject to the terms and conditions defined in the  *
 *  file 'LICENSE', which is part of this source code package.       *
 * ================================================================= */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "sequence.h"
#include "bwtindex.h"
#include "lcparray.h"
#include "graphics.h"

#define VERSION "0.8.2"

//#define BENCHMARK 1
#if defined(unix) && defined(BENCHMARK)
#include "unistd.h"
#include "string.h"
#endif

#ifdef _MSC_VER
#define PAUSE_AT_EXIT 1
#endif

#define MATCH_TYPE_CHAR "EAU"

struct returnBuildBwt
{
    unsigned int first_sa_position;
    char *texto;
};

struct mem
{
    size_t ref_pos;
	size_t query_pos;
	size_t size;
};

struct Node {
    int value;
    int weight;
    struct Node *pointer;
    long index;
};

// Function to insert a new object into the list
void insert(struct Node **L, int o, int wi, long index){
    struct Node *newNode = (struct Node *)malloc(sizeof(struct Node));
    newNode->value = o;
    newNode->weight = wi;
    newNode->pointer = NULL;
    newNode->index = index;

    if (*L == NULL) {
        *L = newNode;
    } else {
        struct Node *current = *L;
        if(current->value > newNode->value) {
                newNode->pointer = current;
                *L = newNode;
                return ;
            }
        while (current->pointer != NULL) {
            if(current->pointer->value >  o){
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
void delete(struct Node **L, int o, int wi) {
    struct Node *current = *L;
    struct Node *prev = NULL;

    while (current != NULL && current->value != o && current->weight != wi) {
        prev = current;
        current = current->pointer;
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
struct Node *next(struct Node *L, int o, int wi) {
    struct Node *current = L;
    struct Node *result = NULL;

    while (current != NULL) {
        if (current->value > o) {
            return current;
        }
        current = current->pointer;
    }

    return result;
}

// Function to find the largest element strictly smaller than o in the list
struct Node *prev(struct Node *L, int o, int wi) {
    struct Node *current = L;
    struct Node *result = NULL;
    
    if(current != NULL)
        if(current->value >= o)
            return NULL;
    

    while (current != NULL) {
        if (current->pointer == NULL)
            return current;

        if (current->pointer->value >= o) {
            return current;
        }
        current = current->pointer;
    }

    return result;
}

// Function to find the maximal element in the list
struct Node *max(struct Node *L) {
    struct Node *current = L;
    struct Node *result = NULL;

    while (current != NULL) {
        if(current->pointer == NULL)return current;
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
struct Node *newNode(int value, int weight, struct Node *pointer) {
    struct Node *newNode = (struct Node *)malloc(sizeof(struct Node));
    newNode->weight = weight;
    newNode->value = value;
    newNode->pointer = pointer;
    if(pointer != NULL)
        newNode->index = pointer->index;
    else newNode->index = -1;
    return newNode;
}

void his(struct Node *nodes, size_t n) {
    struct Node *L = NULL;
    struct Node **node = (struct Node **)malloc(n * sizeof(struct Node *));
    
    for (int i = 0; i < n; i++) {
        node[i] = NULL; // Initialize the node array
    }

    for (int i = 0; i < n; i++) {
        int oi = nodes[i].value;
        int wi = nodes[i].weight;

        struct Node *s = prev(L, oi, wi);
        struct Node *t = NULL;

        if (s == NULL) {
            t = min(L);
        }else{
            t = next(L, s->value, s->weight);
        }

        while(t != NULL){
            int temp_weight = 0;
            if(s != NULL) temp_weight = s->weight;
            if(temp_weight + wi < t->weight)
                break;
            
            int last_value = t->value;
            int last_weight = t->weight;

            delete(&L, t->value, t->weight);
            t = next(L, last_value, last_weight);
        }

        if(t == NULL || oi < t->value){
            int temp_weight = 0;
            if(s != NULL) temp_weight = s->weight;
            insert(&L, oi,wi + temp_weight, i);
            node[i] = newNode(oi,wi + temp_weight, s);
        }

        int teste = 0;
    }

    // Find the maximal element in L
    struct Node *maxNode = max(L);
    
    // Recover the LIS using the node array
    long index = maxNode->index;
    while (index != -1) {
        printf("%d\n", node[index]->value);
        index = node[index]->index;
    }

    // Free the memory of the list and node array
    freeList(L);
    free(node);
}

struct UniqueReferencesResult {
    struct mem **final_common_mems;
    size_t *common_mems_size;
};

struct UniqueReferencesResult unique_queries(struct mem **common_mems,size_t *common_mems_size,size_t numSeqs, size_t numRefs){
	struct mem **final_common_mems = (struct mem **)malloc((numSeqs-numRefs)*sizeof(struct mem*)); 

	for(size_t iter=0; iter<(numSeqs-numRefs); iter++){
		final_common_mems[iter] = (struct mem*)malloc(sizeof(struct mem));
	}

	size_t count = 0;
	for(size_t iterador= numSeqs - numRefs; iterador<(numSeqs - numRefs)*2;iterador++){
		count = 0;
		for(size_t iterador2 = 0; iterador2 < common_mems_size[iterador]; iterador2++){
			if(iterador2 == 0){
				if(common_mems[iterador][iterador2 + 1].query_pos != common_mems[iterador][iterador2].query_pos){
					final_common_mems[iterador-(numSeqs - numRefs)] = (struct mem*)realloc(final_common_mems[iterador-(numSeqs - numRefs)], (count+1) * sizeof(struct mem)); 
					final_common_mems[iterador-(numSeqs - numRefs)][count].query_pos = common_mems[iterador][iterador2].query_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].ref_pos = common_mems[iterador][iterador2].ref_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].size = common_mems[iterador][iterador2].size;
					count++;	
				}
				continue;
			}

			if(iterador2 == common_mems_size[iterador]-1){
				if(common_mems[iterador][iterador2 - 1].query_pos != common_mems[iterador][iterador2].query_pos){
					final_common_mems[iterador-(numSeqs - numRefs)] = (struct mem*)realloc(final_common_mems[iterador-(numSeqs - numRefs)], (count+1) * sizeof(struct mem)); 
					final_common_mems[iterador-(numSeqs - numRefs)][count].query_pos = common_mems[iterador][iterador2].query_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].ref_pos = common_mems[iterador][iterador2].ref_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].size = common_mems[iterador][iterador2].size;
					count++;
				}
				continue;
			}

			if(common_mems[iterador][iterador2 + 1].query_pos != common_mems[iterador][iterador2].query_pos){
				if(common_mems[iterador][iterador2 - 1].query_pos != common_mems[iterador][iterador2].query_pos){
					final_common_mems[iterador-(numSeqs - numRefs)] = (struct mem*)realloc(final_common_mems[iterador-(numSeqs - numRefs)], (count+1) * sizeof(struct mem)); 
					final_common_mems[iterador-(numSeqs - numRefs)][count].query_pos = common_mems[iterador][iterador2].query_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].ref_pos = common_mems[iterador][iterador2].ref_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].size = common_mems[iterador][iterador2].size;
					count++;
				}
			}			
		}
		common_mems_size[iterador-(numSeqs - numRefs)]=count;
	}
	struct UniqueReferencesResult result;
    result.final_common_mems = final_common_mems;
    result.common_mems_size = common_mems_size;

    return result;
}


struct UniqueReferencesResult unique_references(struct mem **common_mems,size_t *common_mems_size,size_t numSeqs, size_t numRefs){
	struct mem **final_common_mems = (struct mem **)malloc(((numSeqs-numRefs))*sizeof(struct mem*)); 

	for(size_t iter=0; iter<(numSeqs-numRefs); iter++){
		final_common_mems[iter] = (struct mem*)malloc(sizeof(struct mem));
	}

	size_t count = 0;
	for(size_t iterador= numSeqs - numRefs; iterador<(numSeqs - numRefs)*2;iterador++){
	//for(size_t iterador= 0; iterador<(numSeqs - numRefs);iterador++){
		count = 0;
		for(size_t iterador2 = 0; iterador2 < common_mems_size[iterador]; iterador2++){
			if(iterador2 == 0){
				if(common_mems[iterador][iterador2 + 1].ref_pos != common_mems[iterador][iterador2].ref_pos){
					final_common_mems[iterador-(numSeqs - numRefs)] = (struct mem*)realloc(final_common_mems[iterador-(numSeqs - numRefs)], (count+1) * sizeof(struct mem)); 
					final_common_mems[iterador-(numSeqs - numRefs)][count].query_pos = common_mems[iterador][iterador2].query_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].ref_pos = common_mems[iterador][iterador2].ref_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].size = common_mems[iterador][iterador2].size;
					count++;	
				}
				continue;
			}

			if(iterador2 == common_mems_size[iterador]-1){
				if(common_mems[iterador][iterador2 - 1].ref_pos != common_mems[iterador][iterador2].ref_pos){
					final_common_mems[iterador-(numSeqs - numRefs)] = (struct mem*)realloc(final_common_mems[iterador-(numSeqs - numRefs)], (count+1) * sizeof(struct mem)); 
					final_common_mems[iterador-(numSeqs - numRefs)][count].query_pos = common_mems[iterador][iterador2].query_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].ref_pos = common_mems[iterador][iterador2].ref_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].size = common_mems[iterador][iterador2].size;
					count++;
				}
				continue;
			}

			if(common_mems[iterador][iterador2 + 1].ref_pos != common_mems[iterador][iterador2].ref_pos){
				if(common_mems[iterador][iterador2 - 1].ref_pos != common_mems[iterador][iterador2].ref_pos){
					final_common_mems[iterador-(numSeqs - numRefs)] = (struct mem*)realloc(final_common_mems[iterador-(numSeqs - numRefs)], (count+1) * sizeof(struct mem)); 
					final_common_mems[iterador-(numSeqs - numRefs)][count].query_pos = common_mems[iterador][iterador2].query_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].ref_pos = common_mems[iterador][iterador2].ref_pos;
					final_common_mems[iterador-(numSeqs - numRefs)][count].size = common_mems[iterador][iterador2].size;
					count++;
				}
			}			
		}
		common_mems_size[iterador-(numSeqs - numRefs)]=count;
	}
	struct UniqueReferencesResult result;
    result.final_common_mems = final_common_mems;
    result.common_mems_size = common_mems_size;

    return result;
}


// Comparison function for qsort
int compare(const void *a, const void *b)
{
    const struct mem *mem1 = (const struct mem *)a;
    const struct mem *mem2 = (const struct mem *)b;
    
    if (mem1->ref_pos < mem2->ref_pos)
        return -1;
    else if (mem1->ref_pos > mem2->ref_pos)
        return 1;
    else
        return 0;
}

struct UniqueReferencesResult getCommonResult(struct UniqueReferencesResult result1, struct UniqueReferencesResult result2, size_t tamanho){
    struct UniqueReferencesResult result3;
	result3.final_common_mems = (struct mem **)malloc(((tamanho))*sizeof(struct mem*));

	for(size_t iter=0; iter<(tamanho); iter++){
		result3.final_common_mems[iter] = (struct mem*)malloc(sizeof(struct mem));
	}

	for(size_t i = 0; i < tamanho; i++){
		size_t j1 = 0, j2 = 0, count = 0;
		while(j1 < result1.common_mems_size[i] && j2 < result2.common_mems_size[i]){
			if(result1.final_common_mems[i][j1].ref_pos == result2.final_common_mems[i][j2].ref_pos){
				result3.final_common_mems[i] = (struct mem*)realloc(result3.final_common_mems[i], (count+1) * sizeof(struct mem)); 
				result3.final_common_mems[i][count].query_pos = result1.final_common_mems[i][j1].query_pos;
				result3.final_common_mems[i][count].ref_pos = result1.final_common_mems[i][j1].ref_pos;
				result3.final_common_mems[i][count].size = result1.final_common_mems[i][j1].size;
				count++;
				j1++;
				j2++;
			}
			if(result1.final_common_mems[i][j1].ref_pos > result2.final_common_mems[i][j2].ref_pos){
				j2++;
			}
			if(result1.final_common_mems[i][j1].ref_pos < result2.final_common_mems[i][j2].ref_pos){
				j1++;
			}
		}
		result1.common_mems_size[i] =  count;
	}

	result3.common_mems_size = result1.common_mems_size;
	return result3;
}

void GetMatches(int numRefs, int numSeqs, int matchType, int minMatchSize, int bothStrands, char *outFilename, int argCircular){
	FILE *matchesOutputFile;
	int i, s, depth, matchSize, numMatches, refId;
	unsigned int j, textsize, refPos;
	long long int sumMatchesSize, totalNumMatches, totalAvgMatchesSize;
	unsigned int topPtr, bottomPtr, prevTopPtr, prevBottomPtr, savedTopPtr, savedBottomPtr, n;
	char c, *text;
	char *refsTexts[1];
	unsigned int refsTextSizes[1];
	unsigned char *lcpArray;
	int progressCounter, progressStep;
	#ifdef DEBUGMEMS
	char *refText;
	int refSize;
	#endif
	#if defined(unix) && defined(BENCHMARK)
	char command[32];
	int commretval;
	#endif
	printf("> Using options: minimum M%cM length = %d ; strand = %s\n", MATCH_TYPE_CHAR[matchType], minMatchSize,(bothStrands==0)?"forward only":"forward + reverse");
	matchesOutputFile=fopen(outFilename,"w");
	FILE *filenon = fopen(strcat(outFilename, "-zero-mems.txt"), "w");
	// Check if the file was opened successfully
    if (filenon == NULL) {
        printf("Error opening the file.\n");
        return ; // Return an error code
    }

	if(matchesOutputFile==NULL){
		printf("\n> ERROR: Cannot create output file <%s>\n",outFilename);
		exit(-1);
	}
	printf("> Building index for reference sequence");
	if(numRefs==1) printf(" \"%s\"", (allSequences[0]->name));
	else printf("s");
	printf(" (%u Mbp) ...\n",(allSequences[0]->size)/1000000U);
	fflush(stdout);
	text=(allSequences[0]->chars);
	textsize=(allSequences[0]->size);
	refId=0;
	refsTexts[0]=text;
	refsTextSizes[0]=textsize;
	lcpArray=NULL;
	struct returnBuildBwt rbwt = FMI_BuildIndex(refsTexts,refsTextSizes,1,&lcpArray,1, text, argCircular);
	text = rbwt.texto;
	//unsigned int first_sa_position =  ;
	i=BuildSampledLCPArray(text,textsize,lcpArray,minMatchSize,1, rbwt.first_sa_position, argCircular);
	if(lcpArray!=NULL) free(lcpArray);
	#ifndef DEBUGMEMS
	FreeSequenceChars(allSequences[0]);
	#else
	refText=(allSequences[0]->chars);
	refSize=(allSequences[0]->size);
	#endif
	#if defined(unix) && defined(BENCHMARK)
	sprintf(command,"memusgpid %d &",(int)getpid());
	commretval=system(command);
	#endif
	printf("> Matching query sequences against index ...\n");
	fflush(stdout);
	totalNumMatches=0;
	totalAvgMatchesSize=0;

	int *marked = (int*)calloc(FMI_GetBWTSize()-1, sizeof(int));
	struct mem *common_mems[(numSeqs-numRefs)*2]; 
	for(size_t iter=0; iter<(numSeqs-numRefs)*2; iter++){
		common_mems[iter] = (struct mem*)malloc(sizeof(struct mem));
	}
	size_t common_mems_size[(numSeqs-numRefs)*2];

	int flag_rep = 0;
	size_t second_lap = 0;

	for(i=numRefs;i<numSeqs;i++){ // process all queries
		LoadSequenceChars(allSequences[i]);
		//edited tinha esta linha
		text=(allSequences[i]->chars);
		textsize=(allSequences[i]->size);
		progressStep=(textsize/10);
		common_mems_size[i-numRefs+second_lap] = 0;
		for(s=0;s<=bothStrands;s++){ // process one or both strands
			if(s==0){ // forward strand
				printf(":: \"%s\" ",(allSequences[i]->name));
				if(flag_rep == 1)
					fprintf(matchesOutputFile,">%s\n",(allSequences[i]->name));
			} else { // reverse strand
				ReverseComplementSequence(text,textsize); // convert to reverse strand
				printf(":: \"%s Reverse\" ",(allSequences[i]->name));
				fprintf(matchesOutputFile,">%s Reverse\n",(allSequences[i]->name));
			}
			fflush(stdout);
			progressCounter=0;
			matchSize=0;
			numMatches=0;
			sumMatchesSize=0;
			depth=0;
			topPtr=0;
			//edited
			if(argCircular == 1)
				bottomPtr=FMI_GetBWTSize()-1;
			else bottomPtr=FMI_GetBWTSize();
			prevTopPtr=topPtr;
			prevBottomPtr=bottomPtr;
			for(j=textsize;j!=0;){
				j--;
				if(progressCounter==progressStep){ // print progress dots
					putchar('.');
					fflush(stdout);
					progressCounter=0;
				} else progressCounter++;
				while( (n=FMI_FollowLetter(text[j],&topPtr,&bottomPtr))==0 ){ // when no match exits, follow prefix links to broaden the interval
					topPtr = prevTopPtr; // restore pointer values, because they got lost when no hits exist
					bottomPtr = prevBottomPtr;
					depth = GetEnclosingLCPInterval(&topPtr,&bottomPtr); // get enclosing interval and corresponding destination depth
					if( depth == -1 ) break; // can happen for example when current seq contains 'N's but the indexed reference does not
					prevTopPtr = topPtr; // save pointer values in case the match fails again
					prevBottomPtr = bottomPtr;
				}
				depth++;
				if( depth >= minMatchSize ){
					if(matchType==1 && n!=1) continue; // not a MAM if we are looking for one
					savedTopPtr = topPtr; // save the original interval to restore after finished processing MEMs
					savedBottomPtr = bottomPtr;
					prevTopPtr = (bottomPtr+1); // to process the first interval entirely
					prevBottomPtr = bottomPtr;
					matchSize = depth;
					if( j != 0 ) c = text[j-1]; // next char to be processed (to the left)
					else c = '\0';
					while( matchSize >= minMatchSize ){ // process all parent intervals down to this size limit
						for( n = topPtr ; n != prevTopPtr ; n++ ){ // from topPtr down to prevTopPtr
							if( FMI_GetCharAtBWTPos(n) != c ){
								refPos = FMI_PositionInText(n, argCircular);
								#ifndef DEBUGMEMS
								if(numRefs!=1){ // multiple refs
									refId = GetSeqIdFromMergedSeqsPos(&refPos); // get ref id and pos inside that ref
									if( strcmp(allSequences[i]->name, allSequences[refId]->name) != 0 )
										fprintf(matchesOutputFile," %s\t",(allSequences[refId]->name));
									else continue;
								}
								
	
								struct mem my_mem;
								my_mem.ref_pos = 0;
								my_mem.query_pos = 0;
								my_mem.size = 0;

								int flag_mem = 0;
								for(size_t ite=0; ite<matchSize; ite++){
									if((i - numRefs + 1 == 1 && flag_rep==0) || marked[ite+(refPos)] == i - numRefs+second_lap || marked[ite+(refPos)] == i - numRefs + 1+second_lap){
										marked[ite+(refPos)]= i - numRefs + 1+second_lap;
										if(flag_mem == 1)
											my_mem.size++;
										if(flag_mem == 0){
											my_mem.ref_pos = (refPos+1+ite);
											my_mem.query_pos = (j+1+ite);
											my_mem.size = 1;
											flag_mem = 1;
										}
										if(ite == matchSize-1){
											if(flag_rep == 1)
												fprintf(matchesOutputFile,"%zu\t%zu\t%zu\n",(my_mem.ref_pos),(my_mem.query_pos),my_mem.size);
											//if(i == numSeqs -1 ){
												common_mems_size[i-numRefs+second_lap]++;
												common_mems[i-numRefs+second_lap] = (struct mem*)realloc(common_mems[i-numRefs+second_lap], (common_mems_size[i-numRefs+second_lap]) * sizeof(struct mem));
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].size = my_mem.size;
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].ref_pos = my_mem.ref_pos;
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].query_pos = my_mem.query_pos;
												flag_mem = 0;
												my_mem.ref_pos = 0;
												my_mem.query_pos = 0;
												my_mem.size = 0;
										}
										
									}else{
										//marked[ite+(refPos)]= 0;
										if(flag_mem == 1 && my_mem.size >= minMatchSize){
											if(flag_rep == 1)
												fprintf(matchesOutputFile,"%zu\t%zu\t%zu\n",(my_mem.ref_pos),(my_mem.query_pos),my_mem.size);
											//if(i == numSeqs -1 ){
												common_mems_size[i-numRefs+second_lap]++;
												common_mems[i-numRefs+second_lap] = (struct mem*)realloc(common_mems[i-numRefs+second_lap], (common_mems_size[i-numRefs+second_lap]) * sizeof(struct mem));
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].size = my_mem.size;
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].ref_pos = my_mem.ref_pos;
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].query_pos = my_mem.query_pos;
												
											//}
										}
											
										flag_mem = 0;
										my_mem.ref_pos = 0;
										my_mem.query_pos = 0;
										my_mem.size = 0;	
									}	
								}
								#else
								fprintf(matchesOutputFile,"%u\t%d\t%d",(refPos+1),(j+1),matchSize);
								fputc('\t',matchesOutputFile);
								fputc((refPos==0)?('$'):(refText[refPos-1]+32),matchesOutputFile);
								fprintf(matchesOutputFile,"%.*s...%.*s",4,(char *)(refText+refPos),4,(char *)(refText+refPos+matchSize-4));
								fputc(((refPos+matchSize)==refSize)?('$'):(refText[refPos+matchSize]+32),matchesOutputFile);
								fputc('\t',matchesOutputFile);
								fputc((j==0)?('$'):(text[j-1]+32),matchesOutputFile);
								fprintf(matchesOutputFile,"%.*s...%.*s",4,(char *)(text+j),4,(char *)(text+j+matchSize-4));
								fputc(((j+matchSize)==textsize)?('$'):(text[j+matchSize]+32),matchesOutputFile);
								fputc('\n',matchesOutputFile);
								#endif
								numMatches++;
								sumMatchesSize += matchSize;
							}
						}
						for( n = bottomPtr ; n != prevBottomPtr ; n-- ){ // from bottomPtr up to prevBottomPtr
							if( FMI_GetCharAtBWTPos(n) != c ){
								refPos = FMI_PositionInText(n, argCircular);
								#ifndef DEBUGMEMS
								if(numRefs!=1){ // multiple refs
									refId = GetSeqIdFromMergedSeqsPos(&refPos); // get ref id and pos inside that ref
									if( strcmp(allSequences[i]->name, allSequences[refId]->name) != 0 )
										fprintf(matchesOutputFile," %s\t",(allSequences[refId]->name));
									else continue;
								}

								struct mem my_mem;
								my_mem.ref_pos = 0;
								my_mem.query_pos = 0;
								my_mem.size = 0;

								int flag_mem = 0;
								for(size_t ite=0; ite<matchSize; ite++){
									if((i - numRefs + 1 == 1 && flag_rep==0) || marked[ite+(refPos)] == i - numRefs+second_lap || marked[ite+(refPos)] == i - numRefs + 1+second_lap){
										marked[ite+(refPos)]= i - numRefs + 1+second_lap;
										if(flag_mem == 1)
											my_mem.size++;
										if(flag_mem == 0){
											my_mem.ref_pos = (refPos+1+ite);
											my_mem.query_pos = (j+1+ite);
											my_mem.size = 1;
											flag_mem = 1;
										}

										if(ite == matchSize-1){
											if(flag_rep == 1)
												fprintf(matchesOutputFile,"%zu\t%zu\t%zu\n",(my_mem.ref_pos),(my_mem.query_pos),my_mem.size);
											//if(i == numSeqs -1 ){
												common_mems_size[i-numRefs+second_lap]++;
												common_mems[i-numRefs+second_lap] = (struct mem*)realloc(common_mems[i-numRefs+second_lap], (common_mems_size[i-numRefs+second_lap]) * sizeof(struct mem));
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].size = my_mem.size;
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].ref_pos = my_mem.ref_pos;
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].query_pos = my_mem.query_pos;
												flag_mem = 0;
												my_mem.ref_pos = 0;
												my_mem.query_pos = 0;
												my_mem.size = 0;
										}
										
									}else{
										marked[ite+(refPos)]= 0;
										if(flag_mem == 1 && my_mem.size >= minMatchSize){
											if(flag_rep == 1)
											fprintf(matchesOutputFile,"%zu\t%zu\t%zu\n",(my_mem.ref_pos),(my_mem.query_pos),my_mem.size);
											//if(i == numSeqs -1 ){
												common_mems_size[i-numRefs+second_lap]++;
												common_mems[i-numRefs+second_lap] = (struct mem*)realloc(common_mems[i-numRefs+second_lap], (common_mems_size[i-numRefs+second_lap]) * sizeof(struct mem));
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].size = my_mem.size;
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].ref_pos = my_mem.ref_pos;
												common_mems[i-numRefs+second_lap][common_mems_size[i-numRefs+second_lap]-1].query_pos = my_mem.query_pos;
											
											//}
										}
											
										flag_mem = 0;
										my_mem.ref_pos = 0;
										my_mem.query_pos = 0;
										my_mem.size = 0;	
									}	
								}
								#else
								fprintf(matchesOutputFile,"%u\t%d\t%d",(refPos+1),(j+1),matchSize);
								fputc('\t',matchesOutputFile);
								fputc((refPos==0)?('$'):(refText[refPos-1]+32),matchesOutputFile);
								fprintf(matchesOutputFile,"%.*s...%.*s",4,(char *)(refText+refPos),4,(char *)(refText+refPos+matchSize-4));
								fputc(((refPos+matchSize)==refSize)?('$'):(refText[refPos+matchSize]+32),matchesOutputFile);
								fputc('\t',matchesOutputFile);
								fputc((j==0)?('$'):(text[j-1]+32),matchesOutputFile);
								fprintf(matchesOutputFile,"%.*s...%.*s",4,(char *)(text+j),4,(char *)(text+j+matchSize-4));
								fputc(((j+matchSize)==textsize)?('$'):(text[j+matchSize]+32),matchesOutputFile);
								fputc('\n',matchesOutputFile);
								#endif
								numMatches++;
								sumMatchesSize += matchSize;
							}
						}
						prevTopPtr = topPtr;
						prevBottomPtr = bottomPtr;
						matchSize = GetEnclosingLCPInterval(&topPtr,&bottomPtr); // get parent interval and its depth
					}
					topPtr = savedTopPtr;
					bottomPtr = savedBottomPtr;
				}
				prevTopPtr=topPtr; // save pointer values in case there's no match on the next char, and they loose their values
				prevBottomPtr=bottomPtr;
			} // end of loop for all chars of seq
			totalNumMatches += numMatches;
			totalAvgMatchesSize += sumMatchesSize;
			matchSize=(int)((numMatches==0)?(0):(sumMatchesSize/(long long)numMatches));
			if(numMatches==0){
				printf("No matches with %s", (allSequences[i]->name));
				fprintf(filenon, "%s\n", (allSequences[i]->name));
			}
			printf(" (%d M%cMs ; avg size = %d bp)\n",numMatches,MATCH_TYPE_CHAR[matchType],matchSize);
			fflush(stdout);
		} // end of loop for both strands
		FreeSequenceChars(allSequences[i]);
		if(flag_rep == 0 && i == numSeqs -1){
			flag_rep = 1; 
			i=numRefs-1;
			second_lap = numSeqs-numRefs;
		}
	} // end of loop for all queries
	FMI_FreeIndex();
	FreeSampledSuffixArray();
	if((numSeqs-numRefs)!=1){ // if more than one query, print average stats for all queries
		printf(":: Average %d M%cMs found per query sequence (total = %lld, avg size = %d bp)\n",(int)(totalNumMatches/(numSeqs-numRefs)),MATCH_TYPE_CHAR[matchType],totalNumMatches,(int)(totalAvgMatchesSize/totalNumMatches));
	}
	fflush(stdout);
	printf("> Saving M%cMs to <%s> ... ",MATCH_TYPE_CHAR[matchType],outFilename);
	fclose(matchesOutputFile);
	fclose(filenon);
	printf("OK\n");
	fflush(stdout);
	printf("Marked Ref:\n");
	for(int iter=0;iter < FMI_GetBWTSize()-1; iter++){
		printf("%d\t",marked[iter]);
	}

	//common_mems[]
	// struct Node nodeList[common_mems_size[numSequences-numRefs-1]];
	// for(size_t iterador=0;iterador<common_mems_size;iterador++){
	// 	nodeList[iterador].index =  iterador;
	// 	nodeList[iterador].pointer = NULL;
	// 	nodeList[iterador].value = common_mems[numSeqs-numRefs-1][iterador].ref_pos;
	// 	nodeList[iterador].weight = common_mems[numSeqs-numRefs-1][iterador].size;
	// }
	// size_t numNodes = sizeof(nodeList) / sizeof(nodeList[0]);
	// his(nodeList, numNodes);

	// printf("\nEndMarked Ref:\n");

	// struct UniqueReferencesResult result = unique_queries(common_mems, common_mems_size, numSeqs, numRefs);
	
	// for(size_t iterador = 0; iterador < (numSeqs-numRefs); iterador++){
	// 	qsort(result.final_common_mems[iterador], result.common_mems_size[iterador], sizeof(struct mem), compare);
	// }

	// for(size_t iterador = 0; iterador < (numSeqs-numRefs)*2; iterador++){
	// 	qsort(common_mems[iterador], common_mems_size[iterador], sizeof(struct mem), compare);
	// }

	// //struct UniqueReferencesResult result2 = unique_references(result.final_common_mems, result.common_mems_size, numSeqs, numRefs);
	// struct UniqueReferencesResult result1 = unique_references(common_mems, common_mems_size, numSeqs, numRefs);

	// struct UniqueReferencesResult result2 = getCommonResult(result, result1, (numSeqs-numRefs));

	// int flag_common = 0;
	// size_t sum = 0;
	// size_t heaviest_sum = 0;
	// size_t heaviest_sum_pos = 0;
	// size_t flag_start_sum_pos = 0;
	// //TODO: falta ver o caso do ultimo MEM

	// for(size_t iterador = 0; iterador < result2.common_mems_size[0]; iterador++){
	// 	flag_common = 0;
	// 	for(size_t iterador2 = numSeqs-numRefs ;iterador2 < (numSeqs-numRefs)*2; iterador2++){
	// 		if(iterador == result2.common_mems_size[0]-1)
	// 		{
	// 			flag_common = 1;
	// 			break;
	// 		}
	// 		if(result2.final_common_mems[iterador2][iterador].query_pos > result2.final_common_mems[iterador2+1][iterador].query_pos){
	// 			flag_common = 1;
	// 			break;
	// 		}		
	// 	}
	// 	sum += result2.final_common_mems[0][iterador].size;

	// 	if(flag_common == 1){
	// 		if(sum > heaviest_sum){
	// 			heaviest_sum = sum;
	// 			if(flag_start_sum_pos == 0){
	// 				flag_start_sum_pos = iterador;
	// 			}
	// 			flag_start_sum_pos = 1;
	// 		}
	// 		flag_start_sum_pos = 0;
	// 		sum = 0;
	// 	}
	// }
	// printf("Tirei as queries repetidas!\n");
}

typedef struct _MEMInfo {
	char refName[65];
	int refPos;
	int queryPos;
	int size;
} MEMInfo;

int MEMInfoSortFunction(const void *a, const void *b){
	char *charsa, *charsb;
	int diff;
	diff = 0;
	charsa = (((MEMInfo *)a)->refName);
	charsb = (((MEMInfo *)b)->refName);
	while((diff=(int)((*charsa)-(*charsb)))==0 && (*charsa)!='\0'){ charsa++; charsb++; }
	if(diff==0){
		diff = ((((MEMInfo *)a)->refPos) - (((MEMInfo *)b)->refPos));
		if(diff==0) diff = ((((MEMInfo *)a)->queryPos) - (((MEMInfo *)b)->queryPos));
	}
	return diff;
}


// NOTE: the spacing of the MUMmer output format is in the form: "  <max_ref_name_size> <9_spaces_number> <9_spaces_number> <9_spaces_number>"
// NOTE: in multi-ref format (4 columns) the name of the ref is considered only up to the 1st space char
void SortMEMsFile(char *memsFilename){
	FILE *memsFile, *sortedMemsFile;
	char c, *sortedMemsFilename, seqname[256];
	int numMems, maxNumMems, numSeqs, refpos, querypos, memsize, formatNumFields, n;
	MEMInfo *memsArray;
	printf("> Sorting MEMs from <%s> ",memsFilename);
	fflush(stdout);
	if((memsFile=fopen(memsFilename, "r"))==NULL){
		printf("\n> ERROR: Cannot read input file\n");
		exit(-1);
	}
	c=fgetc(memsFile);
	if(c!='>'){
		printf("\n> ERROR: Invalid MEMs file\n");
		exit(-1);
	}
	while(c=='>'){
		c=fgetc(memsFile);
		while(c!='\n' && c!=EOF) c=fgetc(memsFile);
		c=fgetc(memsFile);
	}
	if(c==EOF){
		printf("\n> ERROR: No MEMs inside file\n");
		exit(-1);
	}
	formatNumFields=0;
	while(1){
		while(c==' ' || c=='\t') c=fgetc(memsFile);
		if(c!='\n'){
			formatNumFields++;
			while(c!=' ' && c!='\t' && c!='\n' && c!=EOF) c=fgetc(memsFile);
		}
		if(c=='\n' || c==EOF) break;
	}
	if(formatNumFields!=3 && formatNumFields!=4){
		printf("\n> ERROR: Invalid MEMs file format\n");
		exit(-1);
	}
	rewind(memsFile);
	if(formatNumFields==4) printf("(multiple references) ");
	printf("...\n");
	sortedMemsFilename=AppendToBasename(memsFilename,"-sorted.txt");
	if((sortedMemsFile=fopen(sortedMemsFilename,"w"))==NULL){
		printf("> ERROR: Cannot write output file\n");
		exit(-1);
	}
	seqname[0]='\0';
	refpos=-1;
	querypos=-1;
	memsize=-1;
	numSeqs=0;
	numMems=0;
	maxNumMems=0;
	memsArray=NULL;
	while(1){
		c=fgetc(memsFile);
		if(c=='>' || c==EOF){
			if(numSeqs!=0){
				printf("(%d MEMs)\n",numMems);
				fflush(stdout);
				qsort(memsArray,numMems,sizeof(MEMInfo),MEMInfoSortFunction);
				fprintf(sortedMemsFile,">%s\n",seqname);
				while(numMems!=0){
					numMems--;
					refpos=memsArray[numMems].refPos;
					querypos=memsArray[numMems].queryPos;
					memsize=memsArray[numMems].size;
					if(formatNumFields==4) fprintf(sortedMemsFile," %s\t",(memsArray[numMems].refName));
					fprintf(sortedMemsFile,"%d\t%d\t%d\n",refpos,querypos,memsize);
				}
			}
			if(c==EOF) break;
			numMems=0;
			n=fscanf(memsFile," %255[^\n]\n",seqname);
			printf(":: '%s' ... ",seqname);
			fflush(stdout);
			numSeqs++;
			continue;
		} else ungetc(c,memsFile);
		if(numMems==maxNumMems){
			maxNumMems+=1024;
			memsArray=(MEMInfo *)realloc(memsArray,(maxNumMems*sizeof(MEMInfo)));
		}
		if(formatNumFields==4) n=fscanf(memsFile," %64[^\t ]",(memsArray[numMems].refName));
		else memsArray[numMems].refName[0]='\0';
		n=fscanf(memsFile," %d %d %d ",&refpos,&querypos,&memsize);
		if(n!=3){
			printf("\n> ERROR: Invalid format\n");
			getchar();
			exit(-1);
		}
		memsArray[numMems].refPos=refpos;
		memsArray[numMems].queryPos=querypos;
		memsArray[numMems].size=memsize;
		numMems++;
	}
	printf("> Saving sorted MEMs to <%s> ...\n",sortedMemsFilename);
	fflush(stdout);
	fclose(memsFile);
	fclose(sortedMemsFile);
	free(sortedMemsFilename);
	free(memsArray);
	printf("> Done!\n");
	#ifdef PAUSE_AT_EXIT
	getchar();
	#endif
	exit(0);
}

// TODO: when supporting multiple references, also support MEMs file in proper format (4 fields per line, starting with ref name)
void CreateMemMapImage(char *memsFilename){
	FILE *memsFile;
	char c, seqname[256], **seqsNames, *imageFilename;
	int i, numSeqs, *seqsSizes, numMems, refPos, queryPos, memSize, strand, k;
	numSeqs=1;
	for(i=1;i<numSequences;i++){ // count number of references (seqs inside ref file)
		if((allSequences[i]->fileid)==(allSequences[0]->fileid)) numSeqs++;
	}
	if(numSeqs>1){
		printf("\n> ERROR: Support for visualizing multiple reference sequences is not implemented yet.\n");
		printf("\tIf you require this feature, please request it to the author.\n");
		exit(-1);
	}
	printf("> Processing MEMs from <%s> ...\n",memsFilename);
	fflush(stdout);
	if((memsFile=fopen(memsFilename,"r"))==NULL){
		printf("\n> ERROR: Cannot read file\n");
		exit(-1);
	}
	seqsSizes=(int *)malloc(numSequences*sizeof(int));
	seqsNames=(char **)malloc(numSequences*sizeof(char *));
	for(i=0;i<numSequences;i++){
		seqsSizes[i]=(allSequences[i]->size);
		seqsNames[i]=(allSequences[i]->name);
	}
	InitializeRefAlignmentImage(seqsSizes,numSequences);
	seqname[0]='\0';
	refPos=-1;
	queryPos=-1;
	memSize=-1;
	strand=0;
	numSeqs=0;
	numMems=0;
	while(1){
		c=fgetc(memsFile);
		if(c=='>' || c==EOF){
			if(numSeqs!=0){
				printf("(%d MEMs)\n",numMems);
				fflush(stdout);
			}
			if(c==EOF) break;
			numMems=0;
			k=fscanf(memsFile," %255[^\n]\n",seqname);
			printf(":: '%s' ... ",seqname);
			fflush(stdout);
			i=0; // check if seq name ends with string "Reverse"
			while(seqname[i]!='\0') i++;
			if(i>7) i-=7;
			for(k=0;k<8;k++,i++) if(("Reverse"[k])!=seqname[i]) break;
			if(k==8) strand=1; // reverse strand of the previous sequence
			else {
				strand=0;
				numSeqs++;
				if(numSeqs==numSequences){
					printf("\n> ERROR: MEMs file not generated from this query file (too many sequences)\n");
					exit(-1);
				}
				for(;numSeqs<numSequences;numSeqs++){ // get the sequence with the same name as in the MEMs file
					for(k=0;seqname[k]!='\0';k++) if((seqname[k])!=(seqsNames[numSeqs][k])) break;
					if((seqname[k])=='\0' && (seqsNames[numSeqs][k])=='\0') break;
				}
				if(numSeqs==numSequences){
					printf("\n> ERROR: Sequence name was not found in query file\n");
					exit(-1);
				}
			}
			continue;
		} else ungetc(c,memsFile);
		k=fscanf(memsFile," %d %d %d ",&refPos,&queryPos,&memSize);
		if(k!=3){
			printf("\n> ERROR: Invalid MEM format\n");
			exit(-1);
		}
		if( refPos==0 || queryPos==0 || memSize==0 ){
			printf("\n> ERROR: Invalid MEM values\n");
			exit(-1);
		}
		refPos--; // convert from 1-based position to 0-based position
		queryPos--;
		if(strand==0) DrawRefAlignmentBlock(queryPos,refPos,memSize,numSeqs);
		else { // the position is relative to the rev strand, so convert it to fwd strand
			queryPos=((seqsSizes[numSeqs])-(queryPos+memSize)); // pos is on right, add size to go to left, subtract from seq length
			DrawRefAlignmentBlock(queryPos,refPos,(-memSize),numSeqs);
		}
		numMems++;
	}
	fclose(memsFile);
	imageFilename=AppendToBasename(memsFilename,".bmp");
	FinalizeRefAlignmentImage(seqsNames,imageFilename);
	free(imageFilename);
	free(seqsSizes);
	free(seqsNames);
	DeleteAllSequences();
	printf("> Done!\n");
	#ifdef PAUSE_AT_EXIT
	getchar();
	#endif
	exit(0);
}

// Cleans all invalid characters from a FASTA file (leaving only ACGT) and joins multiple sequences if present
void CleanFasta(char *fastafilename){
	FILE *fastafile, *cleanfile;
	char *cleanfilename;
	unsigned int charcount,invalidcharcount,sequencecount,linesize;
	char c;
	printf("> Opening FASTA file <%s> ... ",fastafilename);
	fflush(stdout);
	if((fastafile=fopen(fastafilename,"r"))==NULL){
		printf("\n> ERROR: FASTA file not found\n");
		exit(-1);
	}
	c=fgetc(fastafile);
	if(c!='>'){
		printf("\n> ERROR: Invalid FASTA file\n");
		exit(-1);
	}
	printf("OK\n");
	cleanfilename=AppendToBasename(fastafilename,"-clean.fasta");
	printf("> Creating clean FASTA file <%s> ... ",cleanfilename);
	fflush(stdout);
	if((cleanfile=fopen(cleanfilename,"w"))==NULL){
		printf("\n> ERROR: Can't write clean FASTA file\n");
		exit(-1);
	}
	fprintf(cleanfile,">%s\n",fastafilename); // filename as label
	linesize=0;
	charcount=0;
	invalidcharcount=0;
	sequencecount=0;
	while(c!=EOF){
		if(c=='A' || c=='C' || c=='G' || c=='T'){
			fputc((int)c,cleanfile);
			charcount++;
			linesize++;
			if(linesize==100){ // split sequence by lines of size 100
				fputc('\n',cleanfile);
				linesize=0;
			}
		} else if(c=='a' || c=='c' || c=='g' || c=='t'){
			fputc((int)(c-32),cleanfile);
			charcount++;
			linesize++;
			if(linesize==100){ // split sequence by lines of size 100
				fputc('\n',cleanfile);
				linesize=0;
			}
		} else if(c=='>'){ // new sequence
			while(c!='\n' && c!=EOF) c=fgetc(fastafile); // skip description of sequence
			sequencecount++;
		} else if(c>32 && c<127) invalidcharcount++; // invalid alphanumeric character
		c=fgetc(fastafile);
	}
	if(linesize!=0) fputc('\n',cleanfile); // ending newline
	fclose(fastafile);
	fclose(cleanfile);
	free(cleanfilename);
	printf(" OK\n");
	printf(":: %u total chars",charcount);
	if(invalidcharcount!=0) printf(" (%u non ACGT chars removed)",invalidcharcount);
	if(sequencecount>1) printf(" ; %u sequences merged",sequencecount);
	printf("\n");
	printf("> Done!\n");
	#ifdef PAUSE_AT_EXIT
	getchar();
	#endif
	exit(0);
}

// TODO: if multiple refs exist draw all refs names bellow ref image, one name per row, and vertical lines on each split point extending to corresponding row (box around names or horizontal line above names right extended to the max of name length or split point; grey line when overlapping)
// TODO: enable option "-v" on normal mode to create image automatically after finding MEMs (directly draw each block inside MEMs finding function) (or if "-v" set, call function to check if any of the input files is a mems file, set new var "memsFile", and pass it to image function)
// TODO: draw gene annotations in image if GFF present (first detect type of all input files, set "memsFile"/"gffFile" variables, add all Fastas to a new array and use it as arg to loadSequences function)
// TODO: remove "baseBwtPos" field from SLCP structure to save memory and benchmark new running times
// TODO: output MUMs (only once in query) and Multi-MEMS (same number in ref and all queries)
int main(int argc, char *argv[]){
	int i, j, n, numFiles, numSeqsInFirstFile, refFileArgNum, memsFileArgNum, refNameSearchArgNum;
	int argMatchType, argBothStrands, argNoNs, argMinMemSize, argMinSeqLen, argCircular, argItself;
	char *outFilename, *isArgFastaFile, *refNameSearch, optionChar;
	printf("[ slaMEM v%s ]\n\n",VERSION);
	if(argc<3){
		printf("Usage:\n");
		printf("\t%s (<options>) <reference_file> <query_file(s)>\n",argv[0]);
		printf("Options:\n");
		printf("\t-mem\tfind MEMs: any number of occurrences in both ref and query (default)\n");
		printf("\t-mam\tfind MAMs: unique in ref but any number in query\n");
		//printf("\t-mum\tfind MUMs: unique both in ref and query\n");
		printf("\t-q\tcalculates circular mems\n");
		printf("\t-e\tcompares a fasta file with itself\n");
		printf("\t-l\tminimum match length (default=20)\n");
		printf("\t-o\toutput file name (default=\"*-mems.txt\")\n");
		printf("\t-b\tprocess both forward and reverse strands\n");
		printf("\t-n\tdiscard 'N' characters in the sequences\n");
		printf("\t-m\tminimum sequence size (e.g. to ignore small scaffolds)\n");
		printf("\t-r\tload only the reference(s) whose name(s) contain(s) this string\n");
		printf("Extra:\n");
		printf("\t-v\tgenerate MEMs map image from this MEMs file\n");
		//printf("\t-s\tsort MEMs file\n");
		//printf("\t-c\tclean FASTA file\n");
		printf("Example:\n");
		printf("\t%s -b -l 10 ./ref.fna ./query.fna\n",argv[0]);
		printf("\t%s -v ./ref-mems.txt ./ref.fna ./query.fna\n",argv[0]);
		return (-1);
	}
	if( ParseArgument(argc,argv,"S",0) ){ // Sort MEMs
		if(argc!=3){
			printf("Usage: %s -s <mems_file>\n\n",argv[0]);
			return (-1);
		}
		SortMEMsFile(argv[2]);
		return 0;
	}
	if( ParseArgument(argc,argv,"C",0) ){ // Clean FASTA file
		if(argc!=3){
			printf("Usage: %s -c <fasta_file>\n\n",argv[0]);
			return (-1);
		}
		CleanFasta(argv[2]);
		return 0;
	}
	n=0;
	isArgFastaFile=(char *)calloc(argc,sizeof(char));
	numFiles=0;
	for(i=1;i<argc;i++){
		if(argv[i][0]=='-'){ // skip arguments for options
			optionChar=argv[i][1];
			if(optionChar>='A' && optionChar<='Z') optionChar=(char)('a' + (optionChar - 'A'));
			if(optionChar=='l' || optionChar=='o' || optionChar=='m' || optionChar=='v') i++; // skip value of option "-l", "-o", "-s", "-v"
			else if(optionChar=='r'){ // skip reference name string (can span through multiple args)
				i++;
				if(i==argc) break;
				j=0;
				optionChar=argv[i][0];
				if(optionChar=='\'' || optionChar=='\"') j=1; // string can be enclosed in quotes
				else optionChar='\0';
				n=0;
				while(argv[i][j]!=optionChar){ // count chars of string
					if(argv[i][j]=='\0'){
						i++; // go to next argument
						if(i==argc) break;
						j=0;
					} else j++;
					n++;
				}
			}
			continue;
		}
		isArgFastaFile[i]=1;
		numFiles++;
	}
	argItself=ParseArgument(argc,argv,"e",0);

	int itself_flag = 0;
	//e a opção esteja selecionada
	if(numFiles==1 && argItself == 1)
		itself_flag = 1;
	else{
		if(numFiles<2) exitMessage("Not enough input sequence files provided");
	}
		
	argNoNs=ParseArgument(argc,argv,"N",0);
	argMinSeqLen=ParseArgument(argc,argv,"M",1);
	if(argMinSeqLen==(-1)) argMinSeqLen=0;
	refNameSearchArgNum=ParseArgument(argc,argv,"R",2);
	if(refNameSearchArgNum==(-1)) refNameSearch=NULL;
	else { // get reference name string
		if(n==0) exitMessage("No reference name string provided");
		refNameSearch=(char *)calloc((n+1),sizeof(char));
		i=refNameSearchArgNum;
		j=0;
		optionChar=argv[i][0];
		if(optionChar=='\'' || optionChar=='\"') j=1;
		else optionChar='\0';
		n=0;
		while(argv[i][j]!=optionChar){
			if(argv[i][j]=='\0'){
				i++;
				if(i==argc) break;
				j=0;
				refNameSearch[n]=' ';
			} else{
				refNameSearch[n]=argv[i][j];
				j++;
			}
			n++;
		}
		refNameSearch[n]='\0';
	}
	memsFileArgNum=ParseArgument(argc,argv,"V",2);
	refFileArgNum=(-1);
	numSeqsInFirstFile=0;
	numFiles=0;
	numSequences=0; // initialize global variable needed by sequence functions
	for(i=1;i<argc;i++){
		if(!isArgFastaFile[i]) continue; // skip options and their arguments
		n=LoadSequencesFromFile(argv[i],((numFiles==0 && memsFileArgNum==(-1))?1:0),((numFiles==0)?1:0),argNoNs,(unsigned int)argMinSeqLen,(numFiles==0)?refNameSearch:NULL);
		if(n!=0) numFiles++;
		if(numFiles==0) exitMessage("No valid sequences found in reference file");
		if(numFiles==1){ // reference file
			if(itself_flag == 1){
				n=LoadSequencesFromFile(argv[i],((numFiles==0 && memsFileArgNum==(-1))?1:0),((numFiles==0)?1:0),argNoNs,(unsigned int)argMinSeqLen,(numFiles==0)?refNameSearch:NULL);
			}
			refFileArgNum=i;
			numSeqsInFirstFile=n;
		}
	}
	free(isArgFastaFile);
	if(refNameSearch!=NULL) free(refNameSearch);
	//if(numFiles==0) exitMessage("No reference or query files provided");
	if(numFiles==1 && itself_flag != 1) exitMessage("No query files provided");
	n=(numSequences-numSeqsInFirstFile);
	if(n==0) exitMessage("No valid query sequences found");
	printf("> %d reference%s and %d quer%s successfully loaded\n",numSeqsInFirstFile,((numSeqsInFirstFile==1)?"":"s"),n,((n==1)?"y":"ies"));
	if(memsFileArgNum!=(-1)){ // Create MEMs image
		CreateMemMapImage(argv[memsFileArgNum]);
		return 0;
	}
	argMatchType=0; // MEMs mode
	if( ParseArgument(argc,argv,"MA",0) ) argMatchType=1; // MAMs mode
	argCircular=ParseArgument(argc,argv,"Q",0);
	argBothStrands=ParseArgument(argc,argv,"B",0);
	argMinMemSize=ParseArgument(argc,argv,"L",1);
	if(argMinMemSize==(-1)) argMinMemSize=20; // default minimum MEM length is 20
	n=ParseArgument(argc,argv,"O",2);
	if(n==(-1)) outFilename=AppendToBasename(argv[refFileArgNum],"-mems.txt"); // default output base filename is the ref filename
	else outFilename=argv[n];
	GetMatches(numSeqsInFirstFile,numSequences,argMatchType,argMinMemSize,argBothStrands,outFilename, argCircular);
	if(n==(-1)) free(outFilename);
	DeleteAllSequences();
	printf("> Done!\n");
	#ifdef PAUSE_AT_EXIT
	getchar();
	#endif
	return 0;
}
