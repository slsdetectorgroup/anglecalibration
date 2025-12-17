#include <exception>
#include <stdexcept>

class NoBasePeakOverLapError : public std::runtime_error {
  public:
    NoBasePeakOverLapError()
        : std::runtime_error("There is no frame where base peak ROI overlaps "
                             "with module region."){};
    ~NoBasePeakOverLapError() = default;
};