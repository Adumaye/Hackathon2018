#ifndef _defclass
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
template<typename vtype>
class myvector
{
  vtype* mydata;
  size_t mysize;
public:
  inline vtype & operator[]( int i ) const {return mydata[i];}
  inline size_t size(){ return mysize; }

  myvector( int size_ ){
    mysize = size_;
    mydata = new vtype[mysize];
    #pragma acc enter data copyin(this)
    #pragma acc enter data create(mydata[0:mysize])
  }
  ~myvector(){
    #pragma acc exit data delete(mydata[0:mysize],this)
    delete [] mydata;
  }
};
#define _defclass
#endif
