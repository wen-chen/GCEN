#ifndef __STRIM_H_
#define __STRIM_H_

#include <algorithm>
#include <cctype>
#include <functional>
#include <locale>
#include <string>

// trim from start
static inline std::string &lstrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
          }));
  return s;
}

// trim from end
static inline std::string &rstrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); })
              .base(),
          s.end());
  return s;
}

// trim from both ends
static inline std::string &strim(std::string &s) { return lstrim(rstrim(s)); }

#endif