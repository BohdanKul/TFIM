// lower_bound/upper_bound example
#include <iostream>     // std::cout
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <vector>       // std::vector
using namespace std;

int main () {
  float myints[] = {0.25,0.75,1.0};
  vector<float> v(myints,myints+3);           // 10 20 30 30 20 10 10 20
  
  for (auto elem=v.begin(); elem!=v.end(); elem++){
	cout << *elem << " ";
  }

  sort (v.begin(), v.end());                // 10 10 10 20 20 20 30 30

  std::vector<float>::iterator low,up;
  low=lower_bound (v.begin(), v.end(), 0.96); //      
  up= upper_bound (v.begin(), v.end(), 0.96); //     

  cout << endl;
  for (auto elem=v.begin(); elem!=v.end(); elem++){
	cout << *elem << " ";
  }
  cout << endl;
  cout << "lower_bound at position " << (low- v.begin()) << " " << *low << '\n';
  cout << "upper_bound at position " << (up - v.begin()) << " " << *up  << '\n';

  return 0;
}
