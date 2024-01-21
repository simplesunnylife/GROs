#ifndef UTIL_LOG_H_
#define UTIL_LOG_H_

#include <string>
#include <iostream>

namespace util {

// runtime information, for monitoring the running situation
void Info(const std::string& info) {
  std::cout << "# INFO #  \t" << info << std::endl;
  return;
}

// for debug information
void Debug(const std::string& info) {
  #ifndef NDEBUG
  std::cout << "# DEBUG # \t" << info << std::endl;
  #endif // NDEBUG
  return;
}

// for error information
void Error(const std::string& info) {
  std::cout << "# ERROR # \t" << info << std::endl;
  return;
}

// essential output only
void Output(const std::string& info) {
  std::cout << "# Output #\t" << info << std::endl;
  return;
}
  
// constexpr std::string TEST_SYMBOL = "[TEST] ";
// constexpr std::string INFO_SYMBOL = "[INFO] ";
// constexpr std::string SYSTEM_SYMBOL = "[SYSM] ";
// constexpr std::string WARN_SYMBOL = "[WARN] ";
// constexpr std::string ERROR_SYMBOL = "[ERROR] ";

// /////////////// Log min-system //////////////
// enum LOG_GRADE {
//   LOG_T_GRADE = 0,  // Test Level
//   LOG_I_GRADE = 1,  // Info Level
//   LOG_S_GRADE = 2,  // System Level
//   LOG_W_GRADE = 3,  // Warning Level
//   LOG_E_GRADE = 4,  // Error Level
// };

// LOG_GRADE LOG_GRADE_SET = LOG_GRADE::LOG_S_GRADE;
// void set_log_grade(LOG_GRADE grade) { LOG_GRADE_SET = grade; }

// void LOG_T(const std::string str) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_T_GRADE) return;
//   std::cout << TEST_SYMBOL << str << "\n";
// }

// template <typename T>
// void LOG_T(const std::string str, T n) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_T_GRADE) return;
//   std::cout << TEST_SYMBOL << str << n << "\n";
// }

// void LOG(const std::string str) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_I_GRADE) return;
//   std::cout << INFO_SYMBOL << str << "\n";
// }

// template <typename T>
// void LOG(const std::string str, T n) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_I_GRADE) return;
//   std::cout << INFO_SYMBOL << str << n << "\n";
// }

// void LOG_S(const std::string str) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_S_GRADE) return;
//   std::cout << SYSTEM_SYMBOL << str << "\n";
// }

// template <typename T>
// void LOG_S(const std::string str, T n) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_S_GRADE) return;
//   std::cout << SYSTEM_SYMBOL << str << n << "\n";
// }

// void LOG_W(const std::string str) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_W_GRADE) return;
//   std::cout << WARN_SYMBOL << str << "\n";
// }

// template <typename T>
// void LOG_W(const std::string str, T n) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_W_GRADE) return;
//   std::cout << WARN_SYMBOL << str << n << "\n";
// }

// void LOG_E(const std::string str) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_E_GRADE) return;
//   std::cout << ERROR_SYMBOL << str << "\n";
// }

// template <typename T>
// void LOG_E(const std::string str, T n) {
//   if (LOG_GRADE_SET > LOG_GRADE::LOG_E_GRADE) return;
//   std::cout << ERROR_SYMBOL << str << n << "\n";
// }

};

#endif // UTIL_LOG_H_