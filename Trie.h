#ifndef INC_TRIE_H // Inclusion guard.
#define INC_TRIE_H // Inclusion guard.

// Alphabet size (# of symbols)
#define ALPHABET_SIZE (30000)

// trie node
struct TrieNode
{
    struct TrieNode *children[ALPHABET_SIZE];
    int countToCompare;
    struct MEM *mem;
    // isEndOfWord is true if the node represents
    // end of a sequence
    bool isEndOfWord;
};

struct MaxPath
{
    struct TrieNode *node;
    int sum;
};

// Returns new trie node (initialized to NULLs)
struct TrieNode *getNode(void)
{
    struct TrieNode *pNode = NULL;

    pNode = (struct TrieNode *)malloc(sizeof(struct TrieNode));

    if (pNode)
    {
        int i;

        pNode->isEndOfWord = false;
        pNode->countToCompare = 0;
        pNode->mem = NULL;

        for (i = 0; i < ALPHABET_SIZE; i++)
            pNode->children[i] = NULL;
    }

    return pNode;
}

bool isEmpty(struct TrieNode *root)
{
    for (int i = 0; i < ALPHABET_SIZE; i++)
        if (root->children[i])
            return false;
    return true;
}

int hash(int value)
{
    const int prime = 10000003;
    long long hash_code = value % ALPHABET_SIZE;
    if (hash_code < 0)
    {
        hash_code += ALPHABET_SIZE;
    }
    hash_code = (hash_code * prime) % ALPHABET_SIZE;
    return (int)hash_code;
}

// If not present, inserts key into trie
// If the key is prefix of trie node, just marks leaf node
void insert(struct TrieNode *root, struct MEM **key, int length)
{
    int level;
    int index;

    struct TrieNode *pCrawl = root;

    for (level = 0; level < length; level++)
    {
        index = hash(key[level]->posRef);
        if (!pCrawl->children[index])
        {
            pCrawl->children[index] = getNode();
        }
        pCrawl->countToCompare = 0;
        pCrawl = pCrawl->children[index];
        pCrawl->mem = key[level];
    }

    // mark last node as leaf
    pCrawl->isEndOfWord = true;
}

void printTrie(struct TrieNode *node, int *buffer, int depth)
{
    FILE *file = fopen("trie.txt", "a");
    if (file == NULL)
    {
        printf("Error: could not open file %s\n", "trie.txt");
        return;
    }

    if (node->isEndOfWord)
    {
        fprintf(file, "[ ");
        int i;
        for (i = 0; i < depth; i++)
        {
            fprintf(file, "%d ", buffer[i]);
        }
        fprintf(file, "]\n");
    }
    fclose(file);
    int i;
    for (i = 0; i < ALPHABET_SIZE; i++)
    {
        if (node->children[i])
        {
            buffer[depth] = i;
            printTrie(node->children[i], buffer, depth + 1);
        }
    }
}

struct MaxPath findMaxPath(struct TrieNode *root)
{
    struct MaxPath maxPath = {NULL, 0};

    for (int i = 0; i < ALPHABET_SIZE; i++)
    {
        if (root->children[i] != NULL)
        {
            struct MaxPath childPath = findMaxPath(root->children[i]);
            int childSum = childPath.sum;
            if (childSum > maxPath.sum)
            {
                maxPath.sum = childSum;
                maxPath.node = childPath.node;
            }
        }
    }

    if (/*root->countToCompare != 0 &&*/ root->mem != NULL)
    {
        int curSum = maxPath.sum + root->mem->len;
        if (curSum > maxPath.sum)
        {
            maxPath.sum = curSum;
            maxPath.node = root;
        }
    }

    return maxPath;
}

void rem(struct TrieNode *root, struct TrieNode **parent_ptr, int aux, int readNumber)
{
    if (!root)
        return;

    // Recursively delete child nodes
    for (int i = 0; i < ALPHABET_SIZE; i++)
    {
        rem(root->children[i], &root->children[i], aux + 1, readNumber);
    }

    if (aux != 0 && root->countToCompare != readNumber)
    {
        if (root->isEndOfWord)
            root->isEndOfWord = false;

        if (isEmpty(root))
        {
            *parent_ptr = NULL;
            free(root);
        }
    }
}

// Returns true if key presents in trie, else false
bool search(struct TrieNode *root, const int *key, int length, int readNumber)
{
    int level;
    int index;
    struct TrieNode *pCrawl = root;

    for (level = 0; level < length; level++)
    {
        index = hash(key[level]);

        if (level != 0)
            pCrawl->countToCompare = readNumber;
        if (!pCrawl->children[index])
            return false;

        pCrawl = pCrawl->children[index];
    }

    pCrawl->countToCompare += readNumber;
    pCrawl->isEndOfWord = true;
    return (pCrawl->isEndOfWord);
}

void buildTrie(struct READ *organism, struct TrieNode *root)
{
    if (organism == NULL)
        return;
    if (organism[3].mems[0] == NULL)
        return;
    int i = 0;
    int aux = 0;
    struct MEM **auxArray = (struct MEM **)malloc((i + 1) * sizeof(struct MEM *));
    int auxSize = 0;
    while (organism[3].mems[i] != NULL)
    {
        aux = i;
        while ( (aux + 1) < organism[3].len && organism[3].mems[aux]->posRef > organism[3].mems[aux + 1]->posRef)
        {
            auxArray = (struct MEM **)realloc(auxArray, (auxSize + 1) * sizeof(struct MEM *));
            auxArray[auxSize] = organism[3].mems[aux];

            insert(root, auxArray, auxSize + 1);
            auxSize++;
            aux++;
        }
        auxArray = (struct MEM **)realloc(auxArray, (auxSize + 1) * sizeof(struct MEM *));
        auxArray[auxSize] = organism[3].mems[aux];
        if ((aux) < organism[3].len)
            insert(root, auxArray, auxSize + 1);
        auxArray = (struct MEM **)realloc(auxArray, (0) * sizeof(struct MEM *));
        auxSize = 0;
        i++;
    }
    // int buffer2[900001];
    // printTrie(root, buffer2, 0);
}

void compareWithTrie(struct READ *organism, int index, struct TrieNode *root)
{
    if (root == NULL)
        return;
    if (organism == NULL)
        return;
    if (organism[index].mems[0] == NULL)
        return;

    int i = 0;
    int aux = 0;
    int *auxArray = (int *)malloc((i + 1) * sizeof(int));
    int auxSize = 0;
    
    while (organism[index].mems[i] != NULL)
    {
        aux = i;
        while (organism[index].mems[aux + 1] != NULL && organism[index].mems[aux]->posRef < organism[index].mems[aux + 1]->posRef)
        {
            auxArray = (int *)realloc(auxArray, (auxSize + 1) * sizeof(int));
            auxArray[auxSize] = organism[index].mems[aux]->posRef;
            search(root, auxArray, auxSize + 1, index);
            auxSize++;
            aux++;
        }
        auxArray = (int *)realloc(auxArray, (auxSize + 1) * sizeof(int));
        auxArray[auxSize] = organism[index].mems[aux]->posRef;
        search(root, auxArray, auxSize + 1, index);
        auxArray = (int *)realloc(auxArray, (0) * sizeof(int));
        auxSize = 0;
        i++;
    }
    rem(root, NULL, 0, index);
    printf("Removi\n");
    int buffer2[ALPHABET_SIZE];
    printTrie(root, buffer2, 0);
}

#endif // Inclusion guard.
