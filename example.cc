#include <iostream>

#include "include/rs/multi_map.h"

//using namespace std;

void RadixSplineExample() {
  std::cout << "Begin RadixSplineExample ...\n";
  // Create random keys.
  std::vector<uint64_t> keys(1e6);
  //std::vector<uint64_t> keys(32);
  std::generate(keys.begin(), keys.end(), rand);
  keys.push_back(6329);
  std::sort(keys.begin(), keys.end());

  // Build RadixSpline.
  uint64_t min = keys.front();
  uint64_t max = keys.back();
  
  rs::Builder<uint64_t> rsb(min, max);
  for (const auto& key : keys) rsb.AddKey(key);
  rs::RadixSpline<uint64_t> rs = rsb.Finalize();

  std::cout << "Keys.size()=" << keys.size() << "; " << std::endl;

  // Search using RadixSpline.
  rs::SearchBound bound = rs.GetSearchBound(6329);
  std::cout << "The search key is in the range: [" << bound.begin << ", "
       << bound.end << ")" << std::endl;
  for (size_t i = bound.begin; i <= bound.end; i++)
  {
    std::cout << "; " << keys[i] ;
  }
  std::cout << std::endl;
  
  auto start = std::begin(keys) + bound.begin, last = std::begin(keys) + bound.end;
  std::cout << "begin=" << *std::begin(keys) << "; start=" << *(std::begin(keys) + bound.begin)
            << "; last=" << *(std::begin(keys) + bound.end) << std::endl;
  std::cout << "The key is at position: "
       << std::lower_bound(start, last, 6329) - std::begin(keys) << std::endl;
}

void MultiMapExample() {
  std::cout << "----------------\n";
  std::vector<std::pair<uint64_t, char>> data = {{1ull, 'a'},
                                       {12ull, 'b'},
                                       {7ull, 'c'},  // Unsorted.
                                       {42ull, 'd'}};
  rs::MultiMap<uint64_t, char> map(std::begin(data), end(data));
  std::cout << "Begin MultiMapExample ...\n";
  std::cout << "find(7): '" << map.find(7)->second << "'" << std::endl;
  std::cout << "lower_bound(3): '" << map.lower_bound(3)->second << "'" << std::endl;
}

int main(int argc, char** argv) {
  std::cout << "Begin the radix spline text ...\n";
  RadixSplineExample();
  MultiMapExample();

  return 0;
}
