#ifndef HELPER_CPP
#define HELPER_CPP
 
using namespace std;

//*****************************************************************************
// Flip a particular bit in a bitstr 
//*****************************************************************************
int FlipBit(const int& BitStr, const int& index){
    return BitStr ^ (1 << index);
}

//*****************************************************************************
// Return a vector of indices that sorts a given vector in decreasing order 
//*****************************************************************************
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}


#endif
