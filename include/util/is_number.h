#ifndef UTIL_IS_NUMBER_H_
#define UTIL_IS_NUMBER_H_


namespace util {

// from stack overflow
bool is_number(const std::string& s) {
  return !s.empty() && std::find_if(s.begin(), 
          s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

};

#endif // UTIL_IS_NUMBER_H_