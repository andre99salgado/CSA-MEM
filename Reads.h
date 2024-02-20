#ifndef INC_READS_H // Inclusion guard.
#define INC_READS_H // Inclusion guard.

struct MEM
{
    int posRef;
    int posQuery;
    int len;
};

struct AuxStruct
{
    int maxNumberOfMems;
    int idOfRead;
};

struct READ
{
    char name[256];
    int len;
    struct MEM **mems;
};

void init_mem(struct MEM *mem, int posRef, int posQuery, int memLen)
{
    mem->posRef = posRef;
    mem->posQuery = posQuery;
    mem->len = memLen;
}

void init_read(struct READ *read, char *name, int len)
{
    strcpy(read->name, name);
    read->len = len;
    read->mems = NULL;
}

void add_mem_to_read(struct READ *read, struct MEM *mem, int max_mems)
{
    if (read->len == 0)
    {
        read->mems = (struct MEM **)malloc(sizeof(struct MEM *) * 1);
    }
    else
    {
        read->mems = (struct MEM **)realloc(read->mems, sizeof(struct MEM *) * max_mems);
    }

    read->mems[read->len] = mem;
    read->len++;
}

void print_organism(struct READ *organism, int currentRead)
{
    printf("\nReads:\n");

    for (int i = 0; i < currentRead; i++)
    {
        printf("\nRead name: %s", organism[i].name);
        printf("Number of MEMs: %d\n", organism[i].len);

        for (int j = 0; j < organism[i].len; j++)
        {
            printf("MEM %d: posRef=%d, posQuery=%d, len=%d\n", j, organism[i].mems[j]->posRef, organism[i].mems[j]->posQuery, organism[i].mems[j]->len);
        }
    }
}

void print_organism_to_file(struct READ *organism, int currentRead, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        printf("Error: could not open file %s\n", filename);
        return;
    }

    fprintf(file, "\nReads:\n");

    for (int i = 0; i < currentRead; i++)
    {
        fprintf(file, "\nRead name: %s", organism[i].name);
        fprintf(file, "Number of MEMs: %d\n", organism[i].len);

        for (int j = 0; j < organism[i].len; j++)
        {
            fprintf(file, "MEM %d: posRef=%d, posQuery=%d, len=%d\n", j, organism[i].mems[j]->posRef, organism[i].mems[j]->posQuery, organism[i].mems[j]->len);
        }
    }

    fclose(file);
}

void free_organism(struct READ *organism, int currentRead)
{
    for (int i = 0; i < currentRead; i++)
    {
        for (int j = 0; j < organism[i].len; j++)
        {
            if (organism[i].mems[j] != NULL)
                free(organism[i].mems[j]);
        }
        if (organism[i].mems != NULL)
            free(organism[i].mems);
    }
}

#endif // Inclusion guard.