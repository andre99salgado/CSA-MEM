int BuildSampledLCPArray(char *text, unsigned int textsize, unsigned char *lcparray, int minlcp, int verbose, unsigned int first_sa_position, int argCircular);
void FreeSampledSuffixArray();
int GetLCP(unsigned int bwtpos);
int GetEnclosingLCPInterval(unsigned int *topptr, unsigned int *bottomptr);
