#ifndef COMPARER_H
#define COMPARER_H

/**
 * @brief Class to comparison of data files
 */
class Comparer{

    /**
     * Count nodes by space and time
     */
    int nx,ny,time;

    /**
     * Step by space
     */
    double hx,hy;

    public:

    /**
     * @brief Constructor
     * @param nx num node by Ox
     * @param ny num node by Oy
     * @param hx step by Ox
     * @param hy step by Oy
     * @param time count steps by time
     */
    Comparer(int nx,int ny,double hx,double hy,int time);

    /**
     * @brief Function to comparison
     * @param file1 first file
     * @param file2 second file
     * @param result file name to output
     */
    void delta(char* file1,char* file2, char* result);
    ~Comparer();
};

#endif
