#include <sstream>
#include <string>
#include <TH1F.h>
namespace utils {
  TH1F* SmoothenHisto(TH1 *h1);
  //template <typename T>
  //std::string to_string_with_precision(const T a_value, const int n = 20);
  template <typename T>
  std::string to_string_with_precision(const T a_value, int n )
  {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
  }
}
